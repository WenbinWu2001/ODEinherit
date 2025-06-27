# Cell Lineage Time Series Dataset

*These data are also available in the R package and can be directly loaded with `data(cell_lineage_data)`.*

This repository contains metadata and time series data for the mother-daughter cell pairs analyzed in our study.

## Metadata (`metadata.csv`)

This file provides lineage and timing information for each cell:

- **`cell_id`**: Unique identifier for each cell.
- **`mother_id`**: Identifier of the mother cell. `"root"` indicates a founder cell (seeded at the start of the experiment).
- **`cell_birth_timepoint`**: Time index at which the cell is born. This corresponds to the column index in the time series matrices (excluding the first column).

## Time Series (`*.csv`)

The following files contain denoised time series data for six proteins:

- `Cdc10.csv`
- `Stb3.csv`
- `CLB5.csv`
- `Whi5.csv`
- `Xbp1.csv`
- `Tup1.csv`

Each file is a matrix where:

- Each **row** represents a cell.
- Each **column** (excluding the first) corresponds to a uniformly spaced time point.
- The **first column** is `cell_id`, which matches the IDs in `metadata.csv`.
- Values are set to `NA` before the cell's birth time.

### Preprocessing Notes

- Time series were denoised using functional PCA.
- Interpolated to 5× resolution via local polynomial smoothing.
- Dataset includes:
  - 25 mother cells
  - 60 daughter cells
  - 240 time points over the course of the experiment (originally 48 time points before interpolation)

## Loading the Data in R

You may use the following code snippet to load the data in R:

```
metadata <- read.csv("metadata.csv", row.names = NULL)

variable_names <- c("Cdc10","Stb3", "CLB5", "Whi5", "Xbp1", "Tup1")
mat_list <- lapply(variable_names, function(var_name) {read.csv(paste0(var_name, ".csv"), row.names = 1)})
names(mat_list) <- variable_names
```