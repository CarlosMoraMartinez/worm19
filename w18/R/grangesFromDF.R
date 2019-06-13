
#' grangesFromDF
#' 
#' Creates a GenomicRanges object from a data frame containing seqID, 
#' chrom_start, chrom_end and strand variables.
#'
#' @param dataframe input dataframe. Must contain variables named:
#' seqID, chrom_start, chrom_end, strand.
#' @param seqnames names for sequences (defaults to "seqID")
#' @return A GenomicRanges object with the same intervals as the dataframe
#' @keywords granges
#' @import GenomicRanges
#' @export 

grangesFromDF <- function(dataframe, seqnames = "seqID"){
  gr <- makeGRangesFromDataFrame(dataframe,
                               keep.extra.columns = FALSE,
                               ignore.strand = FALSE,
                               seqinfo = NULL,
                               seqnames.field = seqnames,
                               start.field = "chrom_start",
                               end.field = "chrom_end",
                               strand.field = "strand",
                               starts.in.df.are.0based = FALSE)
  elementMetadata(gr) <- dataframe
  return(gr)
}