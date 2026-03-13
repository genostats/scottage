#' Score test for Andersen-Gill model
#'
#' @param null.model an object of class nullModel
#' @param X a vector 
#'
#' @export
score.test <- function(null.model, X) {
  if(null.model@float) 
    score_test_float(null.model@ptr, X)
  else
    score_test_double(null.model@ptr, X)
}

