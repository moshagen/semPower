#' sempower
#'
#' Perform a power analysis. This is a wrapper function for a-priori, post-hoc, and compromise power analyses.
#'
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param ... other parameters related to the specific type of power analysis requested
#' @return list
#' @examples
#' \dontrun{
#' 
#' ap <- semPower(type = 'a-priori', 
#'                effect = .08, effect.measure = "RMSEA", 
#'                alpha = .05, beta = .05, df = 200)
#' summary(ph)
#' 
#' ph <- semPower(type = 'post-hoc', 
#'                effect = .08, effect.measure = "RMSEA", 
#'                alpha = .05, N = 250, df = 200)
#' summary(ph)
#' 
#' cp <- semPower(type = 'compromise', 
#'                effect = .08, effect.measure = "RMSEA", 
#'                abratio = 1, N = 250, df = 200)
#' summary(ph)
#' }
#' @seealso [semPower.aPriori()], [semPower.postHoc()], [semPower.compromise()]
#' @export
semPower  <- function(type, ...){
  type <- checkPowerTypes(type)
  switch(type, 
         'a-priori'= {
           semPower.aPriori(...)
         },
         'post-hoc' = {
           semPower.postHoc(...)
         },
         'compromise' = {
           semPower.compromise(...)
         }
         )
}
