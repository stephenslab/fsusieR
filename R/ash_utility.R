######## Adaptation of ash function compute local False Sign Rate
#Utility function from ash
get_exclusions=function(data){
  return((data$s==0 | data$s == Inf | is.na(data$x) | is.na(data$s)))
}

# local False Discovery Rate
#
#' @importFrom ashr cdf_post
#' @importFrom ashr mixcdf
calc_np = function(g,data){
  exclude  =  get_exclusions(data)
  NegativeProb = rep(0,length = length(data$x))
  NegativeProb[!exclude] = cdf_post(g, 0, data)[!exclude] - calc_lfdr(g,data)[!exclude]
  NegativeProb[exclude] = mixcdf(g,0)
  ifelse(NegativeProb<0,0,NegativeProb) #deal with numerical issues that lead to numbers <0
}

# local False Discovery Rate
#
#' @importFrom ashr comp_postprob
#' @importFrom ashr mixprop
#' @importFrom ashr pm_on_zero
#'
calc_lfdr = function(g,data){
  exclude  = get_exclusions(data)
  ZeroProb = rep(0,length =  length(data$x))
  ZeroProb[!exclude] = colSums(comp_postprob(g,data)[pm_on_zero(g),,drop = FALSE])[!exclude]
  ZeroProb[exclude] = sum(mixprop(g)[pm_on_zero(g)])
  ifelse(ZeroProb<0,0,ZeroProb) #deal with numerical issues that lead to numbers <0
}

#' @importFrom ashr compute_lfsr
# local False Sign Rate
calc_lfsr <- function(g, data)
{

    return(compute_lfsr(calc_np(ashr::get_fitted_g(g),
                                data),
                        calc_lfdr(ashr::get_fitted_g(g),
                                  data)
                        )
           )


}
