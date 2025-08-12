#' Cell Lineage Time Series Dataset
#'
#' A dataset containing cell lineage metadata and time series of protein expression for a set of mother–daughter cell pairs.
#'
#' @format A named list with two elements:
#' \describe{
#'   \item{\code{metadata}}{A data frame providing cell lineage and timing information. Columns include:
#'     \itemize{
#'       \item \code{cell_id}: Unique identifier for each cell.
#'       \item \code{mother_id}: Identifier of the mother cell; \code{"root"} indicates a founder cell (seeded at the start of the experiment).
#'       \item \code{cell_birth_timepoint}: Time index at which the cell is born. This corresponds to the column index in the time series matrices (excluding the first column).
#'     }
#'   }
#'   \item{\code{time_series}}{A named list of six protein time series matrices:
#'     \itemize{
#'       \item \code{Cdc10}, \code{Stb3}, \code{CLB5}, \code{Whi5}, \code{Xbp1}, \code{Tup1}
#'     }
#'     Each matrix has rows representing cells and columns representing uniformly spaced time points.
#'     In other words, each row is the time series of a protein for a single cell.
#'     The first column is \code{cell_id}, matching those in \code{metadata}.
#'     Values before each cell's birth are \code{NA}.
#'   }
#' }
#'
#' @details
#' Time series data were denoised using functional principal component analysis and interpolated to 5× temporal resolution using local polynomial smoothing.
#' The dataset includes 25 mother cells and 60 daughter cells, each measured on a common grid of 240 time points (uniformly spaced) across the experiment (originally 48 time points before interpolation).
#' In the paper, these time series are referred to as "trajectories".
#'
#' @examples
#' data(cell_lineage_data)
#' names(cell_lineage_data)
#' head(cell_lineage_data$metadata)
#' dim(cell_lineage_data$time_series$CLB5)
#'
"cell_lineage_data"

