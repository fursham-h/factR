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

  # return if appending length is longer than transcript
  if (sum(BiocGenerics::width(x)) < (start + end)) {
    rlang::abort("Appending length is larger than size of transcript")
  }

  # define global variable
  width <- tmp.fwdcumsum <- tmp.revcumsum <- tmp.end <- tmp.start <- NULL

  # retrieve strand information
  strand <- as.character(BiocGenerics::strand(x))[1]

  # subsequently,we will treat granges as forward stranded
  # so if granges is initially on rev strand, we will swap
  # the length to append from both ends
  headlength <- ifelse(strand == "-", end, start)
  taillength <- ifelse(strand == "-", start, end)



  x <- x %>%
    as.data.frame() %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(
      tmp.fwdcumsum = cumsum(width) - headlength,
      tmp.revcumsum = rev(cumsum(rev(width))) - taillength,
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
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

  return(GRanges)
}
