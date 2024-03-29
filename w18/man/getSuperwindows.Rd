% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSuperwindows.R
\name{getSuperwindows}
\alias{getSuperwindows}
\title{getSuperwindows}
\usage{
getSuperwindows(win)
}
\arguments{
\item{win}{dataframe with window data (must have at least "seqID", 
"chrom_start", "chrom_end", "strand" columns)}
}
\value{
A dataframe with the superwindow coordinates.
}
\description{
Given a data frame object with window data, it creates a GenomicRanges 
object and uses the "reduce" function to get non-redundant intervals
(i.e., superwindows). Then, it converts these data to dataframe again.
}
\seealso{
grangesFromDF
}
\keyword{superwindows}
\keyword{windows,}
