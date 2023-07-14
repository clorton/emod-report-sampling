"""
Quick test of the dtk_to_recombination module.
"""

from argparse import ArgumentParser
from datetime import datetime
import os
from pathlib import Path
from typing import List
import numpy as np
import pandas as pd
import dtk_to_recombination as dtk_recom


def process_by_year(
        csv_filename: Path,
        n_samples_year: int,
        frac_pop: List[int],
        user_subset: bool,
        frac_poly: List[float],
        replicate: int,
        has_fever: bool,
        age_bin_weights: List[float],
        genome_filename: Path):

    """Process reports one year at a time rather than reading the entire CSV into memory at once."""

    fraction_counts, fraction_infections = [], []

    # for each year's worth of reports, sample requested number of infections with given criteria
    for year, reports in reports_by_year(csv_filename):

        print(f"{len(reports)} reports for year {year}")
        print(f"{len(set(reports['IndividualID']))} unique IndividualIDs in the reports")
        reports = normalize_population_indices(reports)
        reports = dtk_recom.filter_emod_infections(reports)     # one report per individual in the given timeframe
        print(f"{len(reports)} remaining rows after 'filter_emod_infections()'")
        n_per_pop = dtk_recom.n_samples_by_pop(reports, n_samples_year, frac_pop)
        print(f"{n_per_pop=}")

        if user_subset is None:
            subset_randomly_or_with_polygenomic_fraction(
                n_per_pop,
                reports,
                frac_poly,
                replicate,
                fraction_counts,
                fraction_infections
            )
        else:
            subset_considering_fever_and_age_bins(
                n_per_pop,
                reports,
                replicate,
                has_fever,
                age_bin_weights,
                fraction_counts,
                fraction_infections
            )

    counts_df = pd.concat(fraction_counts)
    sub_genome_df = pd.concat(fraction_infections)

    sub_genome_df, fpg_list_columns = dtk_recom.split_lists(sub_genome_df, genome_filename, fix_index=False)

    if genome_filename.exists():
        sub_ext_df = merge_samples_with_genomes(sub_genome_df, genome_filename)
    else:
        sub_ext_df = explode_fpg_list_columns(sub_genome_df, fpg_list_columns)

    sub_ext_df = sub_ext_df.sort_values(by=["day", "infIndex"])

    return counts_df, sub_ext_df


def reports_by_year(csv_filename: Path):

    """Return reports year by year rather than reading in the entire CSV at once."""

    current_year = -1   # invalid year
    current_rows = []   # empty list
    with pd.read_csv(csv_filename, chunksize=1024) as reader:               # 1024 rows at a time
        for chunk in reader:                                                # grab 1024 rows
            while any(chunk["year"] != current_year):                       # if any row in the current chunk != current_year
                current_rows.append(chunk[chunk["year"] == current_year])   # append any rows matching the current_year
                rows = pd.concat(current_rows)                              # concatenate collected rows into a dataframe
                if not rows.empty:                                          # if we collected any rows
                    print(f"yielding {len(rows)} rows for year {current_year}...")
                    yield current_year, rows                                # yield the result
                chunk = chunk[chunk["year"] != current_year]                # update the current chunk
                current_year = chunk.iloc[0]["year"]                        # update the current_year
                current_rows = []                                           # reset list
            current_rows.append(chunk[chunk["year"] == current_year])       # append any rows matching the current_year
        rows = pd.concat(current_rows)                                      # concatenate collected rows
        if not rows.empty:                                                  # if we collected any rows
            print(f"yielding {len(rows)} rows for year {current_year}...")
            yield current_year, rows                                        # yield the result

    return


def normalize_population_indices(dataframe: pd.DataFrame):

    """Fix up population indexing (sequential order starting with 0)."""

    population_dict = {j: i for i, j in enumerate(sorted(np.unique(dataframe.population.values)))}

    return dataframe.replace({"population": population_dict})


def subset_randomly_or_with_polygenomic_fraction(
        n_per_pop: List[int],
        recursive_df: pd.DataFrame,
        frac_poly: List[float],
        replicate: int,
        fraction_counts: List[int],
        fraction_infections: List[pd.DataFrame]):

    """Select number of samples per population optionally allocating between monogenomic and polygenomic infections."""

    for pop_id, pop_n_samp in enumerate(n_per_pop):
        print(
            "Subsetting population", pop_id, "with", pop_n_samp, "samples per year."
        )
        pop_df = recursive_df[recursive_df["population"] == pop_id]
        if frac_poly is None:
            frac_poly = []
        elif not isinstance(frac_poly, list):
            frac_poly = [frac_poly]

        if len(frac_poly) > 0:
            # convert single numbers into  list from model config is needed to avoid errors in for loop
            print("Subsetting samples based on polygenomic fraction.")
            for fraction in frac_poly:
                if isinstance(fraction, float):
                    print(
                        "Subset will contain up to ",
                        100 * fraction,
                        "% (",
                        pop_n_samp * fraction,
                        ") polygenomic infections per year. True numbers may be lower than this desired threshold.",
                    )
                    counts_df, sub_genome_df = dtk_recom.subset_by_coi_groups(
                        pop_df, pop_n_samp, fraction, replicate
                    )
        else:
            print("Subset infections will be chosen randomly.")
            counts_df, sub_genome_df = dtk_recom.subset_randomly(
                pop_df, pop_id, pop_n_samp, replicate
            )

        fraction_counts.append(counts_df)
        fraction_infections.append(sub_genome_df)

    return


def subset_considering_fever_and_age_bins(
        n_per_pop: List[int],
        recursive_df: pd.DataFrame,
        replicate: int,
        has_fever: bool,
        age_bin_weights: List[float],
        fraction_counts: List[int],
        fraction_infections: List[pd.DataFrame]
):

    """Choose desired number of samples per population considering fever status."""

    print("Subsetting samples based on epi parameters.")
    for pop_id, pop_n_samp in enumerate(n_per_pop):
        print(f"Subsetting population {pop_id} with {pop_n_samp} samples per year.")
        pop_df = recursive_df[recursive_df["population"] == pop_id]
        print(f"{len(pop_df)} samples for population {pop_id}")
        counts_df, sub_genome_df = dtk_recom.subset_by_epi(
            pop_df, pop_n_samp, replicate, has_fever, age_bin_weights
        )

        fraction_counts.append(counts_df)
        fraction_infections.append(sub_genome_df)

    return


def merge_samples_with_genomes(sub_genome_df: pd.DataFrame, genome_df_fn: Path):

    """I think this does a join of the report with a file with additional genome information."""

    # standard EMOD -> GenEpi inputs with genotype subset
    sub_ext_df = sub_genome_df.explode("recursive_nid").reset_index(drop=True)
    # give new column for the order of the indices
    sub_ext_df = dtk_recom.new_indices_df(sub_ext_df)
    sub_ext_df["recursive_nid"] = pd.to_numeric(sub_ext_df["recursive_nid"])
    # add parent infection column to track super and co-transmissions
    genome_df = pd.read_csv(genome_df_fn)[["infIndex", "nid"]]
    genome_df = genome_df.rename(
        columns={"infIndex": "infHep", "nid": "recursive_nid"}
    )
    sub_ext_df = pd.merge(sub_ext_df, genome_df, on="recursive_nid")

    return sub_ext_df


def explode_fpg_list_columns(sub_genome_df: pd.DataFrame, fpg_list_columns: List[str]):

    """Expands rows with a list of infection IDs (and matching bite and genome IDs) into a single row for each."""

    # standard FPG inputs, explode multiple columns
    sub_ext_df = sub_genome_df.explode(fpg_list_columns)
    sub_ext_df = sub_ext_df.rename(columns={"bite_ids": "infHep"})

    return sub_ext_df


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-r", "--report_filename", type=Path, default=Path("two_years_data.csv"))
    parser.add_argument("-n", "--n_samples_year", type=int, default=100)
    parser.add_argument("-f", "--frac_pop", type=float, action="append", default=[])
    parser.add_argument("-u", "--user_subset", action="store_true")
    parser.add_argument("-e", "--subset_fever", action="store_true")
    parser.add_argument("-a", "--age_bins", type=float, action="append", default=None)
    parser.add_argument("-g", "--genome_filename", type=Path, default=Path("dummy.txt"))
    parser.add_argument("-p", "--frac_poly", type=float, action="append", default=None)
    parser.add_argument("-s", "--replicate", type=int, default=1)

    args = parser.parse_args()

    t0 = datetime.now()
    counts, sub_indices = dtk_recom.subset_indices(
        genome_df_fn       = str(args.genome_filename),
        recursive_df_fn    = str(args.report_filename),
        n_samples_year     = args.n_samples_year,
        replicate          = args.replicate,
        frac_pop           = args.frac_pop,
        frac_poly          = args.frac_poly,
        user_subset        = args.user_subset,
        has_fever          = args.subset_fever,
        age_bin_weights    = args.age_bins
    )
    t1 = datetime.now()
    print(f"Time elapsed (dtk_recom.subset_indices()): {t1 - t0}")

    counts.to_csv("counts.csv", index=False)
    sub_indices.to_csv("sub_indices.csv", index=False)

    t0 = datetime.now()
    counts_test, sub_indices_test = process_by_year(
        args.report_filename,
        args.n_samples_year,
        args.frac_pop,
        args.user_subset,
        args.frac_poly,
        args.replicate,
        args.subset_fever,
        args.age_bins,
        args.genome_filename)
    t1 = datetime.now()
    print(f"Time elapsed (process_by_year()): {t1 - t0}")

    counts_test.to_csv("counts_test.csv", index=False)
    sub_indices_test.to_csv("sub_indices_test.csv", index=False)

    print("Done")
