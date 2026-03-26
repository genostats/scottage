#' Null object of Andersen-Gill score test
#'
#' @description S4 class etc
#'
#' @exportClass nullModel
setClass("nullModel", slots = c(formula = "formula", data = "data.frame", float = "logical", ids = "character", ptr = "externalptr"))

setMethod("show", "nullModel", function(object) {
            cat("An object of class nullModel, corresponding to\n")
            print(formula)
})

# avant de tout passer à la fonction Cpp on décale tous les indices de 1 -> utilisables tels quels en C++ 
nullModel <- function(...) {
  L <- list(...)
  ids <- as.character(L$id)
  L$id <- shift1(L$id)
  L$ev.index <- shift1(L$ev.index)
  L$events <- shift1(L$events)
  L$L.times <- lapply(L$L.times, shift1)
  L$L.events <- lapply(L$L.events, shift1)
  float <- scottage("use.float")
  if(float) 
    ptr <- mk_nullObject_float(L)
  else
    ptr <- mk_nullObject_double(L)
  new("nullModel", formula = L$formula, data = L$data, float = float, ids = ids, ptr = ptr)
}

shift1 <- function(x) as.integer(x) - 1L
