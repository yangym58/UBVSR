#' Calculate values in given rotation matrix
#'
#' @param a A number(scalar) in the first entry
#' @param b A number(scalar) in the second entry
#'
#' @return A list with the elements
#' \item{c}{The value of c in th given rotation matrix.}
#' \item{s}{The value of s in th given rotation matrix.}
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
