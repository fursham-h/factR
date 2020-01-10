#' Compare and classify alternative spliced segments
#'
#' @param exons 
#' In pair-wise mode: GRanges object containing exons for a particular transcript.
#' 
#' In intra-list mode: GRangesList bject containing exons for each transcripts.  
#' Transcripts will be paired and compared based on its gene family or groups. 
#' In order to do so, object has to contain gene_id attribute. Alternatively, user 
#' may provide a dataframe with a list of transcripts and its groupings as a `groupings`
#' argument. See `groupings`.
#' @param ... 
#' In pair-wise mode, argument is a GRanges object containing exons for a 
#' particular transcript as a comparison to `exons`
#' @param groupings 
#' Dataframe describing the groupings of the transcripts in `exons`. Ideally, Transcripts should be 
#' grouped by gene families. Therefore, first column in the dataframe is a list of gene_id or
#' gene_names and the second column is the names of transcripts in `exons` which fall into
#' the gene groupings. This argument is optional if gene_id metadata is present in `exons`
#' 
#'
#' @return
#' In pairwise mode: a GRangesList object with included and skipped segments in `exons` 
#' as compared to ...
#' 
#' In intralist mode: a dataframe with coordinates and information of alternative segments
#' between every transcripts in its groupings.
#' @export
#'
#' @examples
#' # pair-wise comparison
#' compareAS(query_exons[[1]], query_exons[[3]])
#' 
#' #intra-list comparison
#' compareAS(query_exons)
compareAS <- function(exons, ..., groupings = NULL) {
  
  # catch missing args
  mandargs <- c("exons")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variables
  tx.id <- index.y <- index.x <- tx.id.x <- tx.id.y <- NULL
  Gene <- seqnames <- compare.to <- coord <- strand <- NULL
  . <- AS.type <- AS.direction <- NULL
  
  argnames <- as.character(match.call())[-1]

  # if a second GRanges object is given, carry out pairwise comparison
  if (length(list(...)) > 0) {
    # multiple comparisons not supported. can be a feature in the future
    if (length(list(...)) > 1) {
      # idea for multiple comparison. convert grangeslist to df for comparison
      rlang::warn("Multiple comparisons is not yet supported in pair-wise mode. First item in ... used")
    }
    tx2 <- list(...)[[1]]
    if (!is(tx2, "GRanges")){
      # escape out and carry out intra-list comparisons
      rlang::warn(sprintf("`%s` is not a GRanges object. Intra-list comparison mode was executed",
                      argnames[[2]]))
    } else {
      # check object types and return warnings if not GRanges
      if (!is(exons, "GRanges")){
        if(is(exons, "GRangesList")){
          exons <- exons[[1]]
          rlang::warn(sprintf("GRangesList objects are not supported in pair-wise mode. First item in `%s` was used",
                          argnames[[1]]))
        } else {
          #error out
          rlang::abort(sprintf("`%s` is not of class GRanges or GRangesList",
                       argnames[[1]]))
        }
      } 
      
      return(getAS_(exons, tx2))
    }
  }
  
  # run intraList comparison
  
  # check for exons object type
  if (is(exons, "GRanges")){
    rlang::abort(sprintf("%s has to be GRangesList in intraList mode", 
                 argnames[[1]]))
  }
  
  # check for >1 items in list
  if (length(exons) == 1){
    rlang::abort(sprintf("%s has to contain >1 items intraList mode", 
                 argnames[[1]]))
  }
  
  # prepare comparison df based on groupings/metadata
  if (!is.null(groupings)) {
    # get header names and prepare comparisons
    groupname <- names(groupings)[1]
    names(groupings) <- c(groupname, 'tx.id')
    
    if (any(!groupings[[2]] %in% names(exons))){
      missing <- sum(!groupings[[2]] %in% names(exons))
      groupings <- groupings %>%
        dplyr::filter(tx.id %in% names(exons))
      rlang::warn(sprintf("%s transcript(s) in `%s` have missing GRanges in `%s`. These were not analyzed",
                      missing, tail(argnames,1), argnames[[1]]))
    }
    if (nrow(groupings) > nrow(dplyr::distinct(groupings))){
      groupings <- dplyr::distinct(groupings)
      rlang::warn(sprintf("Duplicate ids in `%s` were removed", tail(argnames,1)))
    }
    
    groupings <- groupings %>%
      dplyr::group_by(!!as.symbol(groupname)) %>% 
      dplyr::mutate(index = dplyr::row_number()) %>% 
      dplyr::ungroup() %>%
      dplyr::left_join(., ., by = groupname) %>%
      dplyr::filter(index.y > index.x) %>%
      dplyr::select(-dplyr::starts_with('index')) %>%
      dplyr::rename(tx.id = tx.id.x, compare.to = tx.id.y)
  } else {
    # check for metadata. error out if missing gene_id metadata
    if (is.null(unlist(exons)$gene_id)) {
      msg <- 'Unable to create comparison list as `gene_id` attribute is missing from `%s`. 
      Please provide `groupings` df'
      
      rlang::abort(sprintf(msg, argnames[[1]]))
    } else {
      tmp.exons <- unlist(exons)
      groupings <- data.frame(Gene = tmp.exons$gene_id, tx.id = names(tmp.exons)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(Gene) %>% 
        dplyr::mutate(index = dplyr::row_number()) %>% 
        dplyr::ungroup() %>%
        dplyr::left_join(., ., by = 'Gene') %>%
        dplyr::filter(index.y > index.x) %>%
        dplyr::select(-dplyr::starts_with('index')) %>%
        dplyr::rename(tx.id = tx.id.x, compare.to = tx.id.y)
    }
  }
  # run getAS analysis based on comparisons
  
  # create CDS list for all remaining tx
  out <- BiocParallel::bpmapply(function(x, y) {
    ASreport <- getAS_(exons[[x]], exons[[y]]) %>%
       unlist()
    if (is.null(ASreport)) {
      return(NULL)
    }
    ASreport$AS.direction <- names(ASreport)
    names(ASreport) <- NULL
    ASreport <- ASreport %>%
      as.data.frame() %>%
      dplyr::mutate(tx.id = x, compare.to = y,
                    coord = paste0(seqnames,':',start,'-',end)) %>%
      dplyr::select(tx.id, compare.to, coord, strand, AS.type, AS.direction)
    return(ASreport)
  }, groupings[[2]], groupings[[3]],
  BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
  ) %>%
    dplyr::bind_rows()
  out <- groupings %>%
    dplyr::select(tx.id, compare.to) %>%
    dplyr::left_join(out, by = c('tx.id', 'compare.to'))
  
  # add mirror analysis
  out <- out %>%
    dplyr::select(tx.id = compare.to, compare.to = tx.id, coord:AS.direction) %>%
    dplyr::mutate(AS.direction = ifelse(AS.direction == 'included',
                                        'skipped', 'included')) %>%
    dplyr::mutate(AS.type = ifelse(AS.direction == 'included',
                                        toupper(AS.type), tolower(AS.type))) %>%
    dplyr::bind_rows(out,.) %>%
    dplyr::arrange(tx.id, compare.to, 
                   ifelse(strand == '-', dplyr::desc(coord), coord))
  
  return(out)
}





getAS_ <- function(tx1, tx2) {
  
  # define global variables
  revmap <- type <- upcoord <- downcoord <- downdiff <- NULL
  AS.type <- updiff <- countersource <- sourcedown <- NULL
  sourceup <- seqnames <- NULL
  
  # get information on transcripts
  tx1index <- c(1:length(tx1))
  tx2index <- c((length(tx1) + 1):(length(tx1) + length(tx2)))
  strand <- BiocGenerics::strand(tx1)[1] %>% as.character()

  # combine transcripts and disjoin
  disjoint <- suppressWarnings(BiocGenerics::append(tx1, tx2) %>%
    GenomicRanges::disjoin(with.revmap = T))
  
  # return if tx do not have common segments
  if (length(disjoint) == tail(tx2index,1)) {
    return(NULL)
  }
  
  # annotate position of alt segments
  disjoint <- disjoint %>%
    as.data.frame() %>%
    dplyr::mutate(type = ifelse(lengths(revmap) == 2, "cons", "alt"))

  # return NULL GRanges if there is no alternative segments
  if (!"alt" %in% disjoint$type) {
    return(NULL)
  }
  constartindex <- min(which(disjoint$type == "cons"))
  if (constartindex > 1) {
    disjoint[1:constartindex - 1, ]$type <- "up"
  }

  conslastindex <- max(which(disjoint$type == "cons"))
  if (conslastindex < nrow(disjoint)) {
    disjoint[(conslastindex + 1):nrow(disjoint), ]$type <- "down"
  }

  # prepare dataframe for classifcation
  disjoint <- disjoint %>%
    dplyr::mutate(source = ifelse(type != "cons" & revmap %in% tx1index, 1, 0)) %>%
    dplyr::mutate(source = ifelse(type != "cons" & revmap %in% tx2index, 2, source)) %>%
    dplyr::mutate(countersource = chartr("12", "21", source)) %>%
    dplyr::mutate(sourcedown = dplyr::lag(source, default = 0), sourceup = dplyr::lead(source, default = 0)) %>%
    dplyr::mutate(upcoord = dplyr::lag(end, default = .data$start[1]), downcoord = dplyr::lead(start, default = .data$end[dplyr::n()])) %>%
    dplyr::mutate(updiff = start - upcoord, downdiff = downcoord - end)

  # classify AS based on features
  disjoint <- disjoint %>%
    dplyr::filter(type != "cons") %>%
    dplyr::mutate(AS.type = as.character(NA)) %>%
    dplyr::mutate(AS.type = ifelse(type == "up" & downdiff == 1, "ts", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "up" & downdiff > 1, "FE", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "down" & updiff == 1, "pa", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "down" & updiff > 1, "LE", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "alt" & updiff == 1 & downdiff == 1, "RI", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "alt" & updiff == 1 & downdiff > 1, "SD", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "alt" & updiff > 1 & downdiff == 1, "SA", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "alt" & updiff > 1 & downdiff > 1 & (countersource %in% c(sourcedown, sourceup)), "ME", AS.type)) %>%
    dplyr::mutate(AS.type = ifelse(type == "alt" & updiff > 1 & downdiff > 1 & (!countersource %in% c(sourcedown, sourceup)), "CE", AS.type))

  # convert from alternative first exon to alternative (e.g.) if Ranges is on the neg strand
  if (strand == "-") {
    disjoint <- disjoint %>%
      dplyr::mutate(AS.type = chartr("DAFLtspa", "ADLFpats", AS.type)) %>%
      dplyr::arrange(dplyr::desc(start))
  }

  # convert cases depending on whether segment is found in tx1 or tx2
  # return as GRanges
  disjoint <- disjoint %>%
    dplyr::mutate(AS.type = ifelse(source == 1, toupper(AS.type), tolower(AS.type))) %>%
    dplyr::mutate(source = ifelse(source == 1, 'included', 'skipped')) %>%
    dplyr::select(seqnames:strand, source, AS.type) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'source',keep.extra.columns = T)
  return(disjoint)
}
