#' Resize 5' and 3' ends of a transcript GenomicRanges
#'
#' @param x GenomicRanges object containing exon coordinates from a transcript
#' @param start Length to append from the start of transcript
#' @param end Length to append from the end of transcript
#'
#' @return Appended GenomicRanges object
#' @export
#' @author Fursham Hamid
#'
#' @examples
#' library(GenomicRanges)
#' gr1 <- GRanges(
#'   seqnames = "chr1", strand = c("+", "+", "+"),
#'   ranges = IRanges(
#'     start = c(1, 500, 1000),
#'     end = c(100, 600, 1100)
#'   )
#' )
#'
#' resizeTranscript(gr1, 20, 80)
#' resizeTranscript(gr1, 110, 150)
resizeTranscript <- function(x, start = 0, end = 0) {

  # define global variables
  width <- tmp.fwdcumsum <- tmp.revcumsum <- tmp.end <- tmp.start <- NULL
  
  # check inputs
  if (is(x, "GRanges")) {
    type <- "GR"
    x <- GenomicRanges::GRangesList(x)
  } else if (is(x, "GRangesList")) {
    type <- "GRlist"
  } else {
    rlang::abort("x is not of class GRanges or GRangesList")
  }
  if (length(start) ==1) {
    start <- rep(start, length(x))
  } else if(length(start) != length(x)){
    startlen <- length(start)
    xlen <- length(x)
    rlang::abort(sprintf("Length of `start` (%s) is not equal to `x` (%s)", startlen, xlen))
  }
  if (length(end) ==1) {
    end <- rep(end, length(x))
  } else if(length(end) != length(x)){
    endlen <- length(end)
    xlen <- length(x)
    rlang::abort(sprintf("Length of `end` (%s) is not equal to `x` (%s)", endlen, xlen))
  }
  
  # return if appending length is longer than transcript
  if (any(sum(BiocGenerics::width(x)) < (start + end))) {
    rlang::abort("Appending length is larger than size of transcript")
  }
  

  # retrieve strand information
  strandList <- as.character(S4Vectors::runValue(BiocGenerics::strand(x)))

  # subsequently,we will treat granges as forward stranded
  # so if granges is initially on rev strand, we will swap
  trim_df <- tibble::tibble(
    strand = rep(strandList, lengths(x)),
    start = rep(start, lengths(x)),
    end = rep(end, lengths(x))
  )
  trim_df <- trim_df %>%
    dplyr::mutate(tmp.headlength = ifelse(strand == "-", end, start)) %>%
    dplyr::mutate(tmp.taillength = ifelse(strand == "-", start, end)) %>%
    dplyr::select(-strand,-start, -end)
  
  
  x <- x %>%
    as.data.frame() %>%
    dplyr::bind_cols(trim_df) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(
      tmp.fwdcumsum = cumsum(width) - tmp.headlength,
      tmp.revcumsum = rev(cumsum(rev(width))) - tmp.taillength,
      tmp.start = start,
      tmp.end = end
    ) %>%
    dplyr::filter(tmp.fwdcumsum > 0 & tmp.revcumsum > 0) %>%
    dplyr::mutate(
      start = ifelse(dplyr::row_number() == 1, tmp.end - tmp.fwdcumsum + 1, start),
      end = ifelse(dplyr::row_number() == dplyr::n(), tmp.start + tmp.revcumsum - 1, end)
    ) %>%
    dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
    dplyr::select(-dplyr::starts_with("tmp.")) %>%
    dplyr::ungroup()
  
  # format output
  if (type == "GR") {
    x <- x %>% dplyr::select(-group, -group_name) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  } else if (type == "GRlist") {
    x <- x %>% dplyr::select(-group) %>%
      GenomicRanges::makeGRangesListFromDataFrame(keep.extra.columns = T, 
                                                  split.field = "group_name")
  }
  
  return(x)
}
