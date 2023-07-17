# Stream Process CSV with Pandas

The EMOD malaria model currently can generate a report, daily or monthly, with a line for each infected individual. The GenEpi project wants to use this report to get a random sample of infections each month or year. A five year report is not too large (three years, for example, is 145k lines in 32MB of text). However, a full run of 45 years including burn-in generates 2.8m lines in 459MB of text. Ten of the columns of the report are numeric data which, when stored in NumPy arrays (as Pandas does for columns), isn't much larger than the on disk size.

```python
import pandas as pd
df = pd.read_csv("infIndexRecursive-genomes-df.csv", usecols=["population", "year", "month", "infIndex", "day", "count", "age_day", "fever_status", "recursive_count", "IndividualID"])
```

The above command uses about **_280MB_** of RAM.

(The sample file comes from some of Dan Bridenbecker 's work - \\internal.idm.ctr\IDM2\home\dbridenbecker\output\ModelComparision-FPG-MaxInf=9-Allel_20230607_134936\f63\105\3e3\f631053e-3a05-ee11-aa07-b88303911bc1\output)

However, the remaining columns are lists of IDs which Pandas reads in as 1) lists and 2) of strings. The resulting columns in memory are columns of Python objects which are lists of Python objects (strings). Python objects are terribly expensive in memory with all their associated metadata.

```python
import pandas as pd
df = pd.read_csv("infIndexRecursive-genomes-df.csv")    # reads _all_ columns
```

The above command uses more than **_1GB_** of RAM.

The bottom line is that there are ~10k reports per year, but we only want 100 (maybe 1000) samples each year - 1% or 10% of the total data. Rather than reading the entire CSV into memory at once (full processing requires almost 6GB of RAM) we would like to process a "sampling interval's" worth of data at a time. For example, if we know we want 100 samples per year, let's only read one year's worth of data into memory at a time. How do we do this?

A naïve implementation uses the CSV reader and chunks the incoming lines into individual years. However, the CSV reader does not attempt to determine the type of the incoming data, only returning strings, and the Pandas command to create a DataFrame from a list of rows also skips the type inference which comes when using `Pandas.read_csv()`.

Fortunately, [Pandas provides a `chunksize` parameter in `read_csv()`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html) which makes the function into an iterator which returns `chunksize` rows of data at a time. Now we just need to collect a year's worth of rows even though we do not know a priori how many rows of data are in each year of the report. We solve this problem by using `read_csv()` in a generator function which `yield`s a year's worth of data on each iteration. Then we can choose the requisite number of samples from each year, concatenate the results, and Bob's your uncle. In case you are not familiar with generator functions, the resulting function, below, can be used in a `for var in generator_function():`  for loop. See the second "here" link at the bottom for the code using this code.

```python
def reports_by_year(csv_filename: Path):
 
    """Return reports year by year rather than reading in the entire CSV at once."""
 
    current_year = -1   # invalid year
    current_rows = []   # empty list
    with pd.read_csv(csv_filename, chunksize=1024) as reader:               # 1024 rows at a time
        for chunk in reader:                                                # grab 1024 rows
            while any(chunk["year"] != current_year):                       # if any row in the current chunk != current_year
                current_rows.append(chunk[chunk["year"] == current_year])   # append any rows matching the current_year
                rows = pd.concat(current_rows)                              # concatenate collected rows into a dataframe
                if not rows.empty:                                          # if we collected any rows for the current year...
                    print(f"yielding {len(rows)} rows for year {current_year}...")
                    yield current_year, rows                                # yield the result
                chunk = chunk[chunk["year"] != current_year]                # update the current chunk
                current_year = chunk.iloc[0]["year"]                        # update the current_year
                current_rows = []                                           # reset list
            current_rows.append(chunk[chunk["year"] == current_year])       # append any rows matching the current_year
        rows = pd.concat(current_rows)                                      # concatenate collected rows
        if not rows.empty:                                                  # if have any collected rows remaining for the current year...
            print(f"yielding {len(rows)} rows for year {current_year}...")
            yield current_year, rows                                        # yield the result
 
    return
```

The chunked version of the code uses ~1/20th the RAM of the monolithic version (~300MB vs. ~6GB) and runs approximately 3.3x faster for the test case with forty-five years of data (this is reading from the local filesystem - reading from a network drives needs additional testing).

The generator function to use `read_csv()` with the chunksize  parameter is [here](https://github.com/clorton/emod-report-sampling/blob/main/test.py#L76).

The outer function iterating over each year's data by using the generator function is [here](https://github.com/clorton/emod-report-sampling/blob/main/test.py#L31). It processes a year's worth of data, filtering and sampling, and appends each year to a list. [At the end it concatenates all the DataFrames for the individual years into a single DataFrame](https://github.com/clorton/emod-report-sampling/blob/main/test.py#L62) and performs a little more processing - splitting multiple infections for an individual into separate rows and sorting.

## Notes from Dan Bridenbecker

I had to have an environment with pandas 1.3.5.

The messaging in the output about age-bins is interesting, because we aren't doing anything with them.
