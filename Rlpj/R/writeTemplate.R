#' @title A writing template function for LPJ
#' @description  This function reads a template, and replaces parameters with the
#' provides parameters list.If any parameters values is not provided, the function will
#' set it to the default values. The function assumes that a copy of the template
#'  is already placed in the run folder.
#' @param template1 a character string providing the general model template,
#'  e.g, global.ins. Provide only the file name, not the path
#' @param parameterList  a named list containing the parameters to be calibrated
#' @param runDir a character string indicating path to the run directory
#' @return none
#' @export
#' @keywords Rlpj
#' @author Florian Hartig, Ramiro Silveyra Gonzalez, Maurizio Bagnara
#' @note Based an older code of Istem Fer, Uni Potsdam
#' @examples \dontrun{
#' writeTemplate("global.ins", list(run_lamda_max = 0.5, run_emax= 5),
#'               "/home/lpjRun/runDirectory1")
#' }
writeTemplate <- function(template1 = NULL, parameterList = NULL, runDir = NULL){

  # Checking provided parameters
  if (is.null(runDir) || !file.exists(runDir) ){
    stop("Please provide a valid run directory")
  }
  if (is.null(template1) || !file.exists(file.path(runDir, template1))){
    stop("The provided template  does not exits. Please provide a template name")
  }
  if (is.null(parameterList) || !class(parameterList) =="list"){
    stop("Please provide a parameterList")
  }
  # call the function
  if (grepl("global",template1)){
      # Checking paramaterList and if null, setting to default values
      # 256 physiological parameters can be calibrated right now. Exceptions are the fine roots distributions for all PFTs
      parameterList <- checkParameters(scale= "global", parameterList)
  }else if (grepl("europe", template1)){
      # Checking paramaterList and if null, setting to default values
      # 585 physiological parameters, both general and species-specific, can be calibrated right now.
      # Exceptions are the fine roots distributions for all species and PFTs
      parameterList <- checkParameters(scale = "europe", parameterList)
  }else{
    stop("Cannot recognize the template: neither global nor europe")
  }
  # getting parameters names
  parameterNames <- names(parameterList)
  # looping over parameters ## Faster
  template <- readLines(file.path(runDir, template1))
  for(i in 1:length(parameterNames))  {
    template <- sub(parameterNames[i], parameterList[i], template)
  }
  writeLines(template, file.path(runDir, template1))
}


