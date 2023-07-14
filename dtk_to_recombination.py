#! /usr/bin/env python3

import itertools
import numpy as np
import pandas as pd
import os, sys
from os.path import join, dirname, basename, exists
from ast import literal_eval

# DanB
# from postparent.core.utils import read_json, save_json
from utils import read_json, save_json


#####################################################################################
# functions for matching genotypes to infections
#####################################################################################
def concur_index_expand(concur_list, index_dict):
    if isinstance(concur_list, str):
        concur_list = concur_list.strip("][").split(", ")

    coi_indices = []
    try:
        for inf in concur_list:
            coi_indices.append([a.item() for a in index_dict[int(inf)]])
    except:
        pass

    return list(coi_indices)


def recursive_genomes_infection(infection_df, genome_df):
    """
    Update 2022/11: Handling new input files from EMOD co-transmission that accounts for individuals at time steps.
    Differs from original structure that was centered around biting events.
    """
    genome_sub = genome_df[["infIndex", "infId", "nid"]].rename(
        columns={"infId": "InfectionID"}
    )
    inf_df = infection_df.merge(genome_sub, on="InfectionID")
    # remove infections that are undetectable
    # TODO: Update to a more biologically relevant number with MPG team
    inf_filter = (
        inf_df[inf_df["IRBCs"] > 1000].drop(labels=["nid"], axis=1).drop_duplicates()
    )

    # Can update in the future to contain more columns relevant to the observational model, i.e. clinical symptoms, sampling schemes
    grouping_columns = ["IndividualID", "Time"]
    if "AgeYears" in inf_filter.columns:
        inf_filter["age_day"] = np.round(inf_filter["AgeYears"] * 364.24, 2)
        grouping_columns = grouping_columns + ["age_day"]

    infection_temp = (
        inf_filter.groupby(grouping_columns)
        .agg(
            {
                "InfectionID": lambda x: x.tolist(),
                "infIndex": lambda x: x.tolist(),
                "IRBCs": lambda x: x.tolist(),
            }
        )
        .reset_index()
    )
    print(infection_temp.columns)

    recursive_temp = (
        inf_df[grouping_columns[0:2] + ["nid"]]
        .groupby(grouping_columns[0:2])
        .agg({"nid": lambda x: x.tolist()})
        .reset_index()
    )

    recursive_temp["recursive_count"] = recursive_temp["nid"].str.len()
    recursive_df = infection_temp.merge(
        recursive_temp, on=grouping_columns[0:2]
    ).rename(
        columns={"Time": "day", "nid": "recursive_nid", "infIndex": "genepi_infIndex"}
    )
    # Not ideal, but adding a new column for each infection index with the original Genepi 'infIndex' name to keep all downstream observational model functions in tact.
    recursive_df.reset_index(inplace=True)
    recursive_df = recursive_df.rename(columns={"index": "infIndex"})

    # TODO: Update for a multi-population scenario from NodeID
    if len(np.unique(genome_df.population)) > 1:
        population_dict = dict(
            (j, i) for i, j in enumerate(sorted(np.unique(genome_df.population.values)))
        )
        recursive_df = recursive_df.replace({"population": population_dict})
    else:
        recursive_df["population"] = 0

    # add time grouping variables
    recursive_df["year"] = np.floor(recursive_df["day"] / 365.24)
    if "month" in genome_df.columns:
        recursive_df["month"] = np.floor(recursive_df["day"] / 365.24 * 12)

    print("Last few lines of infIndexRecursive-genomes-df:")
    print(recursive_df.head())
    return recursive_df


def recursive_genomes(infection_df):
    """
    Update 2022/11: Will be deprecated soon; assumes that the concurrent infection ID column has all the information about infections in a person at time of biting. Have shown this is not true and misses infections.
    Identify all genotypes in a biting event by concurrent infections.
    """
    # genome_df = pd.read_csv(genome_df)
    _, concur_dict = recursive_infID_coi(genome_df)

    total_coi = []
    for i in genome_dict.keys():
        total_coi.append([i, genome_dict[i], len(genome_dict[i])])

    coi_df = pd.DataFrame(total_coi).rename(
        {0: "infIndex", 1: "recursive_nid", 2: "recursive_count"}, axis="columns"
    )

    # reset years for better wildcard matching/future plotting
    if min(genome_df["year"].values) > 0:
        genome_df["year"] = genome_df["year"] - min(
            i for i in genome_df["year"].values if i > 0
        )

    # check for epi columns - can expand list for for complex scenarios
    mandatory_cols = ["population", "year", "infIndex", "day"]
    xtra_epi_cols = [x for x in genome_df.columns if x in ["age_day", "fever_status"]]
    if len(xtra_epi_cols) > 0:
        mandatory_cols.extend(xtra_epi_cols)
    genome_df_simp = genome_df[mandatory_cols].drop_duplicates().reset_index(drop=True)
    recursive_df = genome_df_simp.merge(coi_df, how="left")

    print("Viewing last entries of new adjusted COI.")
    print(recursive_df.tail())

    return recursive_df


def recursive_infID_coi(genome_df):
    """
    Fxn compiles information for "true COI" concatenating the genome IDs assigned to the current infection ID and the genome IDs for the concurrent infections.
    """
    print("Checking the parents of concurrent infections.")
    # TODO: Incorporate clinic sampling time into the report to determine if the most recent infection ID would be detectable from the time the infection was acquired.
    # TODO: Incorporate strain proportions somehow to run with the proportional model.

    # remove genomes with a single hepatocyte per infection and no current infections
    infId_match = (
        genome_df.groupby("infIndex").apply(lambda x: x["nid"].unique()).to_dict()
    )
    # concur_df = genome_df[genome_df['concInfIds_count'] > 0]
    # concur_df = genome_df[['infIndex', 'concInfIdxs']].drop_duplicates()
    concur_dict = dict(zip(concur_df["infIndex"], concur_df["concInfIdxs"]))

    coi_gt = {}
    for inf_index in concur_df["infIndex"].unique():
        # inf_gen = infId_match[inf_index] # add itself to the the list before recursively counting
        inf_gen = []
        if inf_index in concur_dict:
            concur_gen = concur_index_expand(concur_dict[inf_index], infId_match)
            concur_gen = [inf_gen] + concur_gen
            inf_gen = list(itertools.chain(*concur_gen))
        coi_gt[inf_index] = inf_gen

    return infId_match, coi_gt


#####################################################################################
# functions for subsetting options
#####################################################################################
def new_indices_df(df):
    df["recursive_nid"] = pd.to_numeric(df["recursive_nid"])
    # give new column for the order of the indices
    new_indices = pd.DataFrame(
        list(enumerate(sorted(set(df.recursive_nid.values)))),
        columns=["recursive_nid", "original_nid"],
    )
    df = df.rename(columns={"recursive_nid": "original_nid"})
    df = pd.merge(df, new_indices, how="left", on="original_nid")
    return df


def filter_emod_infections(recursive_df, duplicate_window="year", duplicate_seed=123):

    """
    Filter report down, via random selection, to one report per individual per month or year.
    """

    # check that day is not sequential in output report
    times = np.unique(recursive_df.day.values)
    if sorted(times) == list(range(min(times), max(times) + 1)):
        print(
            "Warning: Population is output per day, which likely includes individuals being sampled multiple times for the same infections.\
            Will be choosing an infection per individual at random for each year in the report."
        )
        duplicate_window = "year"

    if duplicate_window == "month":
        print("Picking one infection per individual per month.")
        recursive_sub = recursive_df.groupby(["IndividualID", "month"]).sample(
            n=1, random_state=duplicate_seed
        )
    else:
        print("Picking one infection per individual per year.")
        recursive_sub = recursive_df.groupby(["IndividualID", "year"]).sample(
            n=1, random_state=duplicate_seed
        )

    return recursive_sub


def subset_indices(
    genome_df_fn,
    recursive_df_fn,
    n_samples_year,
    replicate,
    frac_pop=None,
    frac_poly=None,
    user_subset=None,
    has_fever=False,
    age_bin_weights=None,
):
    """
    Called by Snakemake to create infection subsets.
    """
    recursive_df = pd.read_csv(
#        recursive_df_fn, engine="python", error_bad_lines=False, encoding="utf-8"
        recursive_df_fn, engine="python", encoding="utf-8"
    )

    # convert strings of lists to actual lists, columns "recursive_nid", "infection_ids", "bite_ids", "genome_ids"
    # or just "recursive_nid" if genome_df_fn exists
    recursive_df, fpg_list_columns = split_lists(recursive_df, genome_df_fn)

    print("Filtering reports to one per individual per year...")
    recursive_df = filter_emod_infections(recursive_df)

    # Determine the number of samples per population per year
    n_per_pop = n_samples_by_pop(recursive_df, n_samples_year, frac_pop)

    print(f"Creating subset index columns per population and year for {replicate} replicate(s).")

    fraction_counts, fraction_infections = [], []
    if user_subset is None:
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
                        counts_df, sub_genome_df = subset_by_coi_groups(
                            pop_df, pop_n_samp, fraction, replicate
                        )
            else:
                print("Subset infections will be chosen randomly.")
                counts_df, sub_genome_df = subset_randomly(
                    pop_df, pop_id, pop_n_samp, replicate
                )

            fraction_counts.append(counts_df)
            fraction_infections.append(sub_genome_df)
    else:
        print("Subsetting samples based on epi parameters.")
        for pop_id, pop_n_samp in enumerate(n_per_pop):
            print(f"Subsetting population {pop_id} with {pop_n_samp} samples per year.")
            pop_df = recursive_df[recursive_df["population"] == pop_id]
            counts_df, sub_genome_df = subset_by_epi(
                pop_df, pop_n_samp, replicate, has_fever, age_bin_weights
            )

            fraction_counts.append(counts_df)
            fraction_infections.append(sub_genome_df)

    counts_df = pd.concat(fraction_counts)
    sub_genome_df = pd.concat(fraction_infections)

    if os.path.exists(genome_df_fn):
        # standard EMOD -> GenEpi inputs with genotype subset
        sub_ext_df = sub_genome_df.explode("recursive_nid").reset_index(drop=True)
        # give new column for the order of the indices
        sub_ext_df = new_indices_df(sub_ext_df)
        sub_ext_df["recursive_nid"] = pd.to_numeric(sub_ext_df["recursive_nid"])
        # add parent infection column to track super and co-transmissions
        genome_df = pd.read_csv(genome_df_fn)[["infIndex", "nid"]]
        genome_df = genome_df.rename(
            columns={"infIndex": "infHep", "nid": "recursive_nid"}
        )
        sub_ext_df = pd.merge(sub_ext_df, genome_df, on="recursive_nid")
    else:
        # standard FPG inputs, explode multiple columns
        sub_ext_df = sub_genome_df.explode(fpg_list_columns)
        sub_ext_df = sub_ext_df.rename(columns={"bite_ids": "infHep"})

    sub_ext_df = sub_ext_df.sort_values(by=["day", "infIndex"])
    return counts_df, sub_ext_df


def split_lists(recursive_df: pd.DataFrame, genome_df_fn: str, fix_index=True):

    """
    Split the lists in the recursive dataframe into multiple columns.
    This is needed for GenEpi inputs, but not for FPG inputs.
    """

    if not os.path.exists(genome_df_fn):
        # minor fixes to FPG input
        fpg_list_columns = ["recursive_nid", "infection_ids", "bite_ids", "genome_ids"]
        for list_column in fpg_list_columns:
            recursive_df[list_column] = [
                x.strip("[]").split(",") for x in recursive_df[list_column]
            ]
        # fix population indexing
        if fix_index:
            population_dict = dict(
                (j, i)
                for i, j in enumerate(sorted(np.unique(recursive_df.population.values)))
            )
            recursive_df = recursive_df.replace({"population": population_dict})
    else:
        recursive_df["recursive_nid"] = [
            x.strip("[]").split(", ") for x in recursive_df["recursive_nid"]
        ]
        fpg_list_columns = None

    return recursive_df, fpg_list_columns


def n_samples_by_pop(recursive_df, n_samples_year, frac_pop=None):
    """
    Returns the number of samples in each population based on the sample weights specified by the user.
    """
    n_pop = len(recursive_df["population"].unique())

    # There is only one population, take all the samples from that population
    if n_pop == 1:
        pop_samples = [n_samples_year]
    else:   # n_pop > 1
        print(f"Report contains {n_pop} populations.")
        if frac_pop is None or len(frac_pop) == 0:
            pop_samples = [np.floor(n_samples_year / n_pop)] * n_pop
            print("User did not specify per population sampling weights.")
            print(f"Allocating samples, {n_samples_year}, equally across populations = {pop_samples}, per year.")
        else:
            if len(frac_pop) == n_pop and np.floor(sum(frac_pop)) == 1:
                pop_samples = np.array(frac_pop) * n_samples_year
                print(f"User specified sampling weights have been applied. Subsetting {pop_samples} samples in respective population order by year.")
            else:
                sys.exit(
                    "User specified sampling weights do not equal 1 or differ in length to the number of samples in the population. \
                \n Update the model config appropriately for",
                    n_pop,
                    "populations.",
                )

    # convert to int
    pop_samples = list(map(int, pop_samples))

    return pop_samples


def subset_by_coi_groups(coi_genome_df, n_samples_year, frac_poly, replicates):

    mono_frac = int(n_samples_year * (1 - frac_poly))
    poly_frac = int(n_samples_year * frac_poly)
    adj_coi_counts = pd.DataFrame(
        {"coi_bool": [True, False], "n_target": [poly_frac, mono_frac]}
    )

    print("Checking for years with less than the specified number of samples per year.")
    coi_genome_df["coi_bool"] = coi_genome_df["recursive_count"].gt(1)

    uniq_infIndex = coi_genome_df.drop_duplicates("infIndex")

    group_counts = (
        uniq_infIndex.groupby(["population", "year", "coi_bool"])
        .size()
        .reset_index(name="n_infections")
    )
    group_counts = group_counts.merge(adj_coi_counts, how="left")
    group_counts["n_samples"] = group_counts[["n_infections", "n_target"]].min(axis=1)
    group_counts.insert(1, "frac_poly", frac_poly)

    # Exploiting the fact that pd.concat() takes a list of dataframes
    # Would prefer to create an empty dataframe with the correct schema and append the selected (sampled) rows at each iteration
    sample_dfs = []
    rep_list = np.arange(replicates)
    for rep in rep_list:
        for year in coi_genome_df["year"].unique():
            for polygenomicity, fraction in [
                (False, mono_frac),
                (True, poly_frac),
            ]:  # polygenomicity? - picking good names for things is hard
                take = min(
                    fraction,
                    int(
                        group_counts[
                            (group_counts["year"] == year)
                            & (group_counts["coi_bool"] == polygenomicity)
                        ]["n_samples"]
                    ),
                )
                rep_df = uniq_infIndex[
                    (uniq_infIndex["year"] == year)
                    & (uniq_infIndex["coi_bool"] == polygenomicity)
                ].sample(take, random_state=int(rep))
                rep_df.insert(0, "subset_replicate", rep)
                rep_df.insert(1, "frac_poly", frac_poly)
                sample_dfs.append(rep_df)
    inf_genome_df = pd.concat(sample_dfs)

    return group_counts, inf_genome_df


def subset_randomly(coi_genome_df, pop, n_samples_year, replicates):
    print(coi_genome_df.head())
    uniq_infIndex = coi_genome_df.drop_duplicates("infIndex")

    group_counts = (
        uniq_infIndex.value_counts(["population", "year"])
        .reset_index()
        .rename(columns={0: "count"})
    )

    sample_dfs = []
    rep_list = np.arange(replicates)
    for rep in rep_list:
        for year in coi_genome_df["year"].unique():
            take = int(
                group_counts[
                    (group_counts["year"] == year) & (group_counts["population"] == pop)
                ]["count"]
            )
            sample_n = take if take < n_samples_year else n_samples_year
            rep_df = uniq_infIndex[
                (uniq_infIndex["year"] == year) & (uniq_infIndex["population"] == pop)
            ].sample(sample_n, random_state=int(rep))
            rep_df.insert(0, "subset_replicate", rep)
            rep_df.insert(1, "frac_poly", "None")
            sample_dfs.append(rep_df)
    inf_genome_df = pd.concat(sample_dfs)

    # add column for COI
    inf_genome_df["coi_bool"] = inf_genome_df["recursive_count"].gt(1)
    # update count column to contain counts of the subsets
    group_counts = (
        inf_genome_df.value_counts(
            ["subset_replicate", "population", "year", "coi_bool"]
        )
        .reset_index()
        .rename(columns={0: "count"})
    )
    group_counts.insert(1, "frac_poly", "None")

    return group_counts, inf_genome_df


def subset_by_epi(recursive_df, n_samples_year, replicates, has_fever=False, age_bin_weights=None):

    """Subset the recursive_df by epi columns. If fever status is included, only subset infections that present fever for the observational model.
    If age bin weights are included, subset infections by age bin for the observational model.
    """

    mandatory_cols = [
        "population",
        "year",
        "infIndex",
        "day",
        "count",
        "recursive_nid",
        "recursive_count",
    ]

    # Find any "extra" columns in the data set that are not mandatory
    demo_cols = [col for col in recursive_df.columns if col not in mandatory_cols]

    # Check that there is at least one demographic column.
    if len(demo_cols) == 0:
        sys.exit(
            "There are no unique demography columns in the data set. Confirm infIndexRecursive-genomes-df.csv contains columns of interest.\
        \n If only subsetting N infections in each population per year, set `user_subset` to False in the run config and update subsetting parameters in the model config."
        )

    if "fever_status" in demo_cols and has_fever:
        print("Fever status filter has been applied. Will only subset infections that present fever for the observational model.")
        recursive_df = recursive_df[recursive_df["fever_status"] == 1]

    rep_list = np.arange(replicates)
    sample_dfs = []

    # Convert age to bins if included as an epi column
    if "age_day" in demo_cols and age_bin_weights is not None:

        days = 364.24
        recursive_df["age_bin"] = pd.cut(
            x=recursive_df["age_day"],
            bins=[
                0,
                int(np.floor(days * 5)),
                int(np.floor(days * 15)),
                int(max(recursive_df["age_day"]) * days),
            ],
            labels=["0-5", "5-15", "15+"]
        )
        print("Age bins have been created for subsetting.")
        print(recursive_df.head())

        # NOTE: This is where to update different subsetting proportions for age related bins
        print(f"{age_bin_weights=}")

        if (len(age_bin_weights) != 3) or (0.99 > sum(age_bin_weights)) or (sum(age_bin_weights) > 1):
            sys.exit(
                "Subsetting is currently set up for three age bins (0-5, 5-15, 15+), and the user did not provide three bins with a sum of one.\
            \n Update the model config or the `subset_by_epi` function in 'workflow/dtk_to_recombination.py' to match lengths."
            )

        n_samp_by_age = list(np.floor(np.array(age_bin_weights) * n_samples_year))
        print("Further subsetting by age bins with the following maximum samples per age bin per year.")
        print(list(zip(["0-5", "5-15", "15+"], n_samp_by_age)))
        n_samp_by_bin = {
            "age_bin": ["0-5", "5-15", "15+"],
            "n_subset_target": n_samp_by_age,
        }

        # calculate number of true and subset target number of infections per age bin
        group_counts = recursive_df.groupby(["population", "year", "age_bin"]).size().reset_index(name="n_infections")
        group_counts = group_counts.merge(pd.DataFrame(n_samp_by_bin), how="left")
        for rep in rep_list:
            for year in recursive_df["year"].unique():
                for bins, n_samp in zip(["0-5", "5-15", "15+"], n_samp_by_age):
                    take = int(min(n_samp, int(group_counts[(group_counts["year"] == year) & (group_counts["age_bin"] == bins)]["n_infections"])))
                    if take > 0:
                        bin_df = recursive_df[(recursive_df["year"] == year) & (recursive_df["age_bin"] == bins)].sample(take, random_state=rep)
                        bin_df.insert(0, "subset_replicate", rep)
                        sample_dfs.append(bin_df)
    else:
        print("No 'age_day' column or no age bin weights subsetting by year.")
        group_counts = recursive_df.groupby(["population", "year"]).size().reset_index(name="n_infections")
        group_counts["n_subset_target"] = n_samples_year
        for rep in rep_list:
            for year in recursive_df["year"].unique():
                take = min(n_samples_year, int(group_counts[group_counts["year"] == year]["n_infections"].iloc[0]))
                year_df = recursive_df[recursive_df["year"] == year].sample(take, random_state=int(rep))
                year_df.insert(0, "subset_replicate", rep)
                sample_dfs.append(year_df)

    inf_genome_df = pd.concat(sample_dfs)

    return group_counts, inf_genome_df


#####################################################################################
# functions for creating summary statistic groupings
#####################################################################################
def split_years(genome_df, recursive_df, user_subset=None):
    output_year = join(dirname(genome_df), "subsets")
    if not exists(output_year):
        os.makedirs(output_year)

    genome_df = pd.read_csv(genome_df)[["infIndex", "nid"]].rename(
        columns={"infIndex": "infHep", "nid": "recursive_nid"}
    )
    recursive_df = pd.read_csv(recursive_df).dropna()
    recursive_df = subset_reframe(recursive_df)
    df = pd.merge(recursive_df, genome_df)
    # df.to_csv(join(output_year, "..", "all-indices.csv"), index=False)

    for year in df["year"].unique():
        genomes_by_year = df[df["year"] == year]
        genomes_by_year = new_indices_df(genomes_by_year)
        # create subset long format recursive_nid dfs
        genomes_by_year.to_csv(
            join(output_year, f"year{int(year)}-indices.csv"), index=False
        )
        # create grouping json for summary stats
        grouping_json = create_summary_groupings(genomes_by_year, "year", user_subset)
        save_json(
            join(output_year, f"year{int(year)}-groups.json"), grouping_json, indent=4
        )
    return


def run_groupings(
    df, timescale, output_groups, output_indices=None, user_subset=None, coi_groups=None
):
    df = pd.read_csv(df)
    # run groupings
    groupings = create_summary_groupings(df, timescale, user_subset, coi_groups)
    save_json(output_groups, groupings, indent=4)

    if df["recursive_nid"].dtype != np.int64:
        df = subset_reframe(df)
        idx_hold = pd.DataFrame(
            columns=["year", "infIndex", "original_nid", "recursive_nid"]
        )
        for year in df["year"].unique():
            df_year = df[df["year"] == year]
            df_year = new_indices_df(df_year)
            idx_hold = pd.concat([idx_hold, df_year])
        idx_hold.to_csv(output_indices, index=False)


def check_epi_cols(df):
    """
    Confirms the inclusion of demographic columns for summary statistic subsetting.
    """
    mandatory_cols = [
        "population",
        "year",
        "month",
        "infIndex",
        "day",
        "count",
        "subset_replicate",
        "frac_poly",
        "original_nid",
        "recursive_nid",
        "recursive_count",
        "age_day",
        "fever_status",
        "age",
        "infHep",
    ]
    demo_cols = [col for col in df.columns if col not in mandatory_cols]
    if len(demo_cols) > 0:
        print(
            "The following column(s) were identified as potential demographic grouping",
            demo_cols,
        )
        # check the new epi columns offer a new subset - otherwise remove from analyses
        n_demo_groups = df[demo_cols].drop_duplicates().shape[0]
        if n_demo_groups < 2:
            print(
                "Demographic columns of",
                demo_cols,
                "did not further subset the population. Ignoring demographic subsets.",
            )
            demo_cols = []
    else:
        print(
            "Demographic groupings specified by the user, but no demographic columns found. \
        \nProceeding with full population summary statistics."
        )
    return demo_cols


def create_summary_groupings(df, timescale, user_subset=None, coi_groups=None):
    """
    Create new groups for running summary statistics on subsets on a yearly or monthly scale.
    Returns a nested dictionary with the grouping variable as the first key.
    The individual grouping keys contain:
        - Columns used in grouping for downstream reformatting
        - Indices of the infections belonging to that subset group for analyses ('infIndex')
    """
    demo_cols = []
    if user_subset is not None:
        demo_cols = check_epi_cols(df)

    if timescale == "month" and "month" not in df.columns:
        df = convert_month(df)

    if coi_groups is not None:
        print(
            "User has specified summary comparisons for",
            len(coi_groups),
            "cutoffs. \
        \n The following bounds are used for adding COI groups:",
            coi_groups,
        )
        df["coi_group"] = pd.cut(
            x=df["recursive_count"],
            bins=[0] + list(map(int, coi_groups)),
            labels=list(map(str, coi_groups)),
        )

    # identify first grouping keys
    subset_check, subset_cols = [timescale], ["subset_replicate", "frac_poly"]
    for col in subset_cols:
        if col in df.columns:
            if len(df[col].unique()) > 1:
                subset_check.append(col)

    # TODO: Figure out how to filter the empty dictionaries for the elegant solution
    # d_groups = ['sampling', 'population', 'demography', 'popxdemo', 'coi_group', 'popxcoi']
    # d = {group: {} for group in d_groups}
    d = {}
    # add to dictionary
    if len(subset_check) > 0:
        d["sampling"] = {}
        d["sampling"]["columns"] = subset_check
    if len(df["population"].unique()) > 1:
        d["population"] = {}
        d["population"]["columns"] = subset_check + ["population"]
    if len(demo_cols) > 0:
        d["demography"] = {}
        d["demography"]["columns"] = subset_check + demo_cols
    if coi_groups is not None:
        d["coi_group"] = {}
        d["coi_group"]["columns"] = subset_check + ["coi_group"]
    if len(df["population"].unique()) > 1:
        if len(demo_cols) > 0:
            d["popxdemo"] = {}
            d["popxdemo"]["columns"] = subset_check + ["population"] + demo_cols
        if coi_groups is not None and len(coi_groups) > 0:
            d["popxcoi"] = {}
            d["popxcoi"]["columns"] = subset_check + ["population", "coi_group"]

    # identify indices for all combinations
    for key in d.keys():
        for name, group in df.groupby(d[key]["columns"]):
            if not isinstance(name, tuple):
                nest_key = str(name)
            else:
                nest_key = "_".join(map(str, name))
            d[key][nest_key] = list(set(group.infIndex.values))
    return d


#####################################################################################
# functions for reformatting genomes dfs
#####################################################################################
def convert_month(df):
    df["month"] = df["day"].values / 365.24 * 12
    df["month"] = df["month"] - np.min(df["month"])
    df["month"] = df["month"].round(0).astype(int)
    return df


def subset_reframe(df):
    df["recursive_nid"] = df["recursive_nid"].apply(literal_eval)
    df = df.explode("recursive_nid").reset_index(drop=True)
    df["recursive_nid"] = pd.to_numeric(df["recursive_nid"])
    return df
