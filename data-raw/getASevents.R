
getASevents <- function(transcript1, transcript2, testedNMD, orf, is_NMD) {

  # prepare list to be returned
  ASlist <- list(
    CE = "NA", MX = "NA", SD = "NA", SA = "NA", FE = "NA", TS = "NA", LE = "NA", PA = "NA", RI = "NA",
    ce = "NA", mx = "NA", sd = "NA", sa = "NA", fe = "NA", ts = "NA", le = "NA", pa = "NA", ri = "NA"
  )
  if (testedNMD == TRUE) {
    ASlist <- utils::modifyList(ASlist, list(NMDcausing = as.character("NA"), NMDcausing.coord = as.character("NA")))
  }

  ASreport <- classifyAS(transcript1, transcript2)
  ASreport_out <- ASreport %>%
    as.data.frame() %>%
    dplyr::group_by(AS) %>%
    dplyr::mutate(coord = paste0(start, "-", end)) %>%
    dplyr::summarise(combcoord = as.character(paste(coord, collapse = ";"))) %>%
    dplyr::select(AS, combcoord)
  ASlist <- utils::modifyList(ASlist, split(ASreport_out$combcoord, ASreport_out$AS))


  if (testedNMD == T & is_NMD == T & length(ASreport) > 0) {
    NMDexonreport <- identifyNMDcausing(ASreport, orf)
    ASlist <- utils::modifyList(ASlist, NMDexonreport)
  }

  return(ASlist)
}
