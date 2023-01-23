numpy <- NULL
tensorly <- NULL
sklearn <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference
  tensorly <<- reticulate::import("tensorly", delay_load = TRUE)
  sklearn <<- reticulate::import("sklearn", delay_load = TRUE)
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
  reticulate::source_python(system.file("source_py.py", package = "spaco"))
}
