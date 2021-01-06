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
#' @examples
#' # sort elements in each GRangesList in descending coordinate order
#' query_exons_desc <- sorteach(query_exons, dplyr::desc(start))
#'
#' # sort elements in each GRangesList in its order in transcript
#' query_exons_exonorder <- sorteach(query_exons_desc, exonorder)
#'
#' # test similarity of query_exons and query_exons_exonorder
#' identical(query_exons, query_exons_exonorder)
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
                   keep.extra.columns = TRUE
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
#' @examples
#' # select first element of each GRangesList item
#' filtereach(query_exons, dplyr::row_number() == 1)
#' @author Fursham Hamid
filtereach <- function(x, ...) {
    group <- NULL
    stopifnot(is(x, "GRangesList"))
    return(x %>% as.data.frame() %>%
               dplyr::group_by(group) %>%
               dplyr::filter(...) %>%
               dplyr::ungroup() %>%
               dplyr::select(-group) %>%
               GenomicRanges::makeGRangesListFromDataFrame(
                   split.field = "group_name",
                   keep.extra.columns = TRUE
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
#' @examples
#' # Create chr:start-end id for each entry
#' mutateeach(query_exons, id = paste0(seqnames, ":", start, "-", end))
#' @author Fursham Hamid
mutateeach <- function(x, ...) {
    group <- NULL
    stopifnot(is(x, "GRangesList"))
    return(x %>% as.data.frame() %>%
               dplyr::mutate(...) %>%
               dplyr::select(-group) %>%
               GenomicRanges::makeGRangesListFromDataFrame(
                   split.field = "group_name",
                   keep.extra.columns = TRUE
               ))
}
