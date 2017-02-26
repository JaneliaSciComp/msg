options <- commandArgs()
MSG_directory <- dirname(substring(options[grep("^--file=", options)], 8)) #As done in fit-hmm.R

MSG_packages <- c("R.methodsS3", "R.oo", "zoo", "HiddenMarkov")

check_package <- function(pkg_name, MSG_directory) {
   if (!require(pkg_name, lib.loc=MSG_directory, character.only=TRUE)) {
      print(paste0("R package ", pkg_name, " not found or incorrectly installed.  Re-installing.\n"))
      install.packages(pkg_name, lib=MSG_directory, dependencies=TRUE, type="source", repos="http://cran.us.r-project.org")
   }
   if (!require(pkg_name, lib.loc=MSG_directory, character.only=TRUE)) {
      print(paste0("Unable to install R package ", pkg_name, ".\n"))
      return(FALSE)
   } else {
      return(TRUE)
   }
}

quit(save="no", status=!all(sapply(MSG_packages, function(x, MSG_directory) {check_package(x, MSG_directory)}, MSG_directory)))
