.onAttach <- function(libname, pkgname) {
  semPowerVersion <- tryCatch(utils::packageDescription("semPower", fields = "Version"),
                            warning = function(w) return(""),
                            error = function(e) return(""))
  txt <- paste0("\n### Welcome to semPower ", semPowerVersion, " ###",
                "\n",
                "\nSee https://github.com/moshagen/semPower for quick examples.",
                "\nSee https://moshagen.github.io/semPower/ for a detailed manual.",
                "\n",
                "\nPlease cite as:",
                "\nMoshagen, M., & Bader, M. (2023). semPower: General Power Analysis for Structural Equation Models.",
                "\nBehavior Research Methods, 56, 2901-2922. https://doi.org/10.3758/s13428-023-02254-7",
                "\n")
  packageStartupMessage(txt)
}


