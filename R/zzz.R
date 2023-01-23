.onLoad <- function(libname, pkgname) {
  reticulate::virtualenv_create(envname = "spaco",
                                packages = c("numpy", "tensorly",
                                             "scikit-learn"))
  reticulate::use_virtualenv("spaco")
  reticulate::source_python(system.file("source_py.py", package = "spaco"))
}
