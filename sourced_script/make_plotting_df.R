#' Create a plotting dataframe from an omicsData object
#'
#' @param omicsData An object of the class 'pepData', 'proData', 'lipidData', or 'metabData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param sample_colname The name of the column which contains the sample ID's
#' @param value_colname The name of the column which contains peak measurement values for a particular sample.
#'
#' @return A dataframe consisting of the joined information in e_data f_data and e_meta within the omicsData object
#' @export
#'
#' @examples
#' library(pmartR)
#' library(pmartRdata)
#' library(datadr)
#' 
#' #make a ddf to pass to trelliscope
#' 
#' data(pep_object)
#' omicsData_long <- make_plotting_df(pep_object)
#' omicsData_ddf <- divide(omicsData_long, by = "Protein")

make_plotting_df <- function(omicsData, sample_colname = "Sample", value_colname = "Value"){
  edata_cname = get_edata_cname(omicsData)

  res <- omicsData$e_data %>% gather(!!rlang::sym(sample_colname), !!rlang::sym(value_colname), -edata_cname)

  if(!is.null(omicsData$f_data)){
    fdata_cname = get_fdata_cname(omicsData)
    res <- res %>% left_join(omicsData$f_data, by = setNames(fdata_cname, sample_colname))
  }

  if(!is.null(omicsData$e_meta)){
    res <- res %>% left_join(omicsData$e_meta, by = setNames(edata_cname, edata_cname))
  }

  res
}
