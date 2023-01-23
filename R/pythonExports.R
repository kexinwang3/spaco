reticulate::use_condaenv("spaco", required = TRUE)
reticulate::source_python(system.file("source_py.py", package = "spaco"))
