numpy <- NULL
tensorly <- NULL
sklearn <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
  tensorly <<- reticulate::import("tensorly", delay_load = TRUE)
  sklearn <<- reticulate::import("sklearn", delay_load = TRUE)
  py <<- reticulate::import_builtins()
}
