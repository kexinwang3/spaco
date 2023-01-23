.onLoad <- function(libname, pkgname) {
  reticulate::install_miniconda(force = TRUE)
  if (!("spaco" %in% reticulate::conda_list()$name)) {
    reticulate::conda_install(envname = "spaco",
                              packages = c("numpy", "tensorly", "scikit-learn"))
  }
  reticulate::use_condaenv("spaco", required = TRUE)
  reticulate::source_python(system.file("source_py.py", package = "spaco"))
}
