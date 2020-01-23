#' Internally sort each element of a GenomicRangesList
#'
#' @param x
#' GRangesList object
#' @param ...
#' Comma separated list of unquoted variable names to sort by. Variables are names of metadata columns
#' found in GRangesList object. Use desc() to sort a variable in descending order.
#' Input can be `exonorder` to sort each element in exon order
#'
#'
#' @return
#' Sorted GRangesList object
#' @export
#'
#' @author Fursham Hamid
sorteach <- function(x, ...) {
  # define global variables
  strand <- group <- NULL
  stopifnot(is(x, "GRangesList"))
  expr <- rlang::quos(...)
  if ("~exonorder" %in% as.character(expr)) {
    index <- which("~exonorder" %in% as.character(expr))
    expr[[index]] <- rlang::quo(ifelse(strand == "-", dplyr::desc(start), start))
  }
  return(x %>% as.data.frame() %>%
    dplyr::arrange(!!!expr) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(
      split.field = "group_name",
      keep.extra.columns = T
    ))
}

#' Internally filter each element of a GenomicRangesList
#'
#' @param x
#' GRangesList object
#' @param ...
#' Logical conditions to filter each element in the GRanges by. Multiple conditions
#' can be provided as comma-delimited inputs
#'
#' @return
#' Filtered GRangesList object
#' @export
#'
#' @author Fursham Hamid
filtereach <- function(x, ...) {
  group <- NULL
  stopifnot(is(x, "GRangesList"))
  return(x %>% as.data.frame() %>%
    dplyr::filter(...) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(
      split.field = "group_name",
      keep.extra.columns = T
    ))
}

#' Internally create or transform metadata of a GenomicRangesList
#'
#' @param x
#' GRangesList object
#' @param ...
#' Name-value pairs of expressions. The name of each argument will be the name
#' of a new metadata column, and the value will be its corresponding value.
#'
#' @return
#' Transformed GRangesList object
#' @export
#'
#' @author Fursham Hamid
mutateeach <- function(x, ...) {
  group <- NULL
  stopifnot(is(x, "GRangesList"))
  return(x %>% as.data.frame() %>%
    dplyr::mutate(...) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(
      split.field = "group_name",
      keep.extra.columns = T
    ))
}
