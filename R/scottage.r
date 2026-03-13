#' Scottage options
#' 
#' @description Set or get options values for package 'scottage'
#' 
#' @param ... options to be defined, using 'name = value', or name(s) of 
#' option(s) to get.
#'
#' @details Currently, only one option can be defined, \code{use.float}
#' The default value is set to 'FALSE'
#' 
#' @return A list with the defined options, or a single element when \code{scottage(name)}
#' is used.
#' @examples
#' old.value <- scottage("use.float")
#' scottage(use.float = TRUE)
#' scottage("use.float")
#' scottage(use.float = old.value)
#' @export
scottage <- function(...) {
  L <- list(...)
  if(length(L) == 0) return(as.list(scottage.opts)) 
  naL <- names(L)
  R <- list()
  if(is.null(naL)) {
    for(x in unlist(L)) R[[x]] <- scottage.opts[[x]]
    if(length(R) == 1) return(unlist(R)) else return(R)
  }
  naOpts <- ls(envir = scottage.opts)
  for(i in seq_along(naL)) {
    if(naL[i] != "") {
      if(!(naL[i] %in% naOpts)) {
        warning(naL[i], " is not a known option")
        next
      }
      assign(naL[i], L[[i]], envir = scottage.opts)
      R[[ naL[i] ]] <- L[[i]]
    } else {
      R[[ L[[i]] ]] <- scottage.opts[[ L[[i]] ]]
    }
  }
  R
}

scottage.opts <- list2env(list(use.float = FALSE))

