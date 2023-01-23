.onLoad <- function(libname, pkgname) {
  user_permission <- utils::askYesNo("Install miniconda?")
  if (isTRUE(user_permission)) {
    reticulate::install_miniconda(force = TRUE)
  } else {
    warning("You should run `reticulate::install_miniconda()` before using this package.")
  }
  if (!("spaco" %in% reticulate::conda_list()$name)) {
    reticulate::conda_install(envname = "spaco",
                              packages = c("numpy", "tensorly", "scikit-learn"))
  }
  reticulate::use_condaenv("spaco", required = TRUE)
  reticulate::source_python(system.file("source_py.py", package = "spaco"))
}
