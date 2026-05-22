#' @keywords internal
"_PACKAGE"

#' @import Matrix
#' @importFrom grDevices chull
#' @importFrom stats kmeans cov var sd pt model.matrix pnorm
#' @importFrom igraph graph_from_adjacency_matrix components make_empty_graph add_edges
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom cli cli_alert_info
#' @importFrom FNN get.knn get.knnx
#' @importFrom data.table data.table rbindlist
#' @importFrom spatstat.geom nnwhich nndist
NULL
