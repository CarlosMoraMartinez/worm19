
#' filterSWbyDist
#' 
#' Filters out the superwindows which distances to the first gene exon
#' are not within a given interval (parameter dist).
#' 
#' @param sw dataframe with superwindows (intervals). 
#' Must have a "dist_to_first_exon" column.
#' @param dist numeric vector of length 2. The first element is the maximum 
#' upstream distance; the second element is the maximum downstream distance.
#' @param include if TRUE, includes partially overlapping sw and trims "width" 
#' to the overlapping region. If FALSE, omits overlapping windows (defaults
#' to TRUE)
#' @return A dataframe with the same columns as the superwindows dataframe.
#' @keywords filters, distance, superwindows
#' @export

filterSWbyDist <- function(sw, dist, include = TRUE){
  if(include == FALSE){
    sw2 <- sw[sw$dist_to_first_exon > dist[1] & 
                sw$dist_to_first_exon < dist[2], ]
  } else {
    ## In case that superwindows that partially overlap with the dist interval 
    ## need to be included
    sw2 <- sw[(((sw$dist_to_first_exon + sw$width) > dist[1]) & 
                 (sw$dist_to_first_exon < 0)) | 
                (((sw$dist_to_first_exon - sw$width) < dist[2]) & 
                   (sw$dist_to_first_exon >= 0)), ]
    exced <- which(sw2$dist_to_first_exon < dist[1] & 
                     sw2$dist_to_first_exon > dist[2])
    
    ## Trim sw length so that it matches with the interval given by dist
    if(length(exced > 0)){
      for(i in 1:length(exced)){
        if(sw[exced[i], "dist_to_first_exon"] > 0){
          sw2[exced[i], "width"] <- sw2[exced[i], "width"] - (abs(sw2[exced[i], "dist_to_first_exon"]) - abs(dist[2]))
        } else {
          sw2[exced[i], "width"] <- sw2[exced[i], "width"] - (abs(sw2[exced[i], "dist_to_first_exon"]) - abs(dist[1]))
        }
      }
    }
  }
  return(sw2)
}