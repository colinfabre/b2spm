.onAttach <- function(libname, pkgname) {
  packageStartupMessage("+--------------------------------------------------------------------------------------------------+")
  packageStartupMessage("|---------------------------------------------- B2SPM ---------------------------------------------|")
  packageStartupMessage("|---------------------------------------------- v1.1 ----------------------------------------------|")
  packageStartupMessage("|##################################################################################################|")
  packageStartupMessage("|-------------------------------------- Dependencies checking -------------------------------------|")

  required_pkgs <- c("utils", "stats", "terra", "sf", "gstat")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) == 0) {
    packageStartupMessage("|-------------------------------------- Dependencies -- ready -------------------------------------|")
  } else {
  packageStartupMessage("|---------------------------------- Dependencies -- /!\ MISSING /!\ ---------------------------------|")
  packageStartupMessage(paste0("|----- Please run first: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\")) -----|"))
  }

  packageStartupMessage("+--------------------------------------------------------------------------------------------------+")
}