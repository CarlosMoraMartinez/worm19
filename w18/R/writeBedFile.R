
#' writeBedFile
#' 
#' Creates a dataframe with bed-like file structure and writes it into 
#' a .bed file without the row names and the column names.
#' @param df input dataframe with "chrom", "chrom_start", "chrom_end"
#' and "strand" variables
#' @param IDs character vector with the same length as the number of rows
#' of the input dataframe, containing the interval names
#' @param fname output file name (.bed will be pasted to it)
#' @keywords bed
#' @export

writeBedFile <- function(df, IDs, fname){
  bed <- data.frame(chrom = df$chrom, chromStart = df$chrom_start, 
                    chromEnd = df$chrom_end, name = IDs, 
                    score = 0, strand = df$strand)
  write.table(bed, file = paste(fname, ".bed", sep="", collapse=""), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
