#' Calculate values in given rotation matrix
#'
#' @param a A number in the first entry
#' @param b A number in the second entry
#'
#' @return The values in given rotation matrix
#' @export
#'
#' @examples given_rotation(3, 4)
given_rotation <- function(a, b)
{  if(a == 0){
  c <- 0
  s <- -sign(b)}
  else if(abs(a) > abs(b)){
    t <- b / a
    u <- sign(a) * sqrt(1 + t * t)
    c <- 1 / u
    s <- -c * t}
  else{
    t <- a / b
    u <- sign(b) * sqrt(1 + t * t)
    s <- -1 / u
    c <- t / u}
  return(list(c = c, s = s))
}
