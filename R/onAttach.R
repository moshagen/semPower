.onAttach <- function(libname, pkgname) {
  semPowerVersion <- tryCatch(utils::packageDescription("semPower", fields = "Version"),
                            warning=function(w) return(""),
                            error=function(e) return(""))
  txt <- paste0("\n### Welcome to semPower ", semPowerVersion, " ###",
                "\n",
                "\nSee https://github.com/moshagen/semPower for quick examples and a detailed manual.",
                "\n",
                "\nPlease cite as:",
                "\nMoshagen, M., & Erdfelder, E. (2016). A new strategy for testing structural equation models.",
                "\nStructural Equation Modeling, 23, 54-60. doi: 10.1080/10705511.2014.950896",
                "\n")
  packageStartupMessage(txt)
}