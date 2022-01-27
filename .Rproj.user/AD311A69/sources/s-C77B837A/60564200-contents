
#'@title Get mixture proportion for mixture normal prior
#'
#'@description
#'@param G_prior mixture normal prior
#'@return vector of mixture proportion

get_pi_G_prior.mixture_normal <- function(G_prior)
{
  out <- G_prior[[1]]$fitted_g$pi
  return(out)
}



#'@title Get mixture proportion for mixture normal prior per scale
#'
#'@description
#'@param G_prior mixture normal prior
#'@return list of vector of mixture proportion

get_pi_G_prior.mixture_normal_per_scale <- function(G_prior)
{
  out <- lapply(G_prior, function(x) x$fitted_g$pi)
  return(out)
}

#'@title Get mixture standard deviations for mixture normal prior
#'
#'@description
#'@param G_prior mixture normal prior
#'@return vector of standard deviations

get_sd_G_prior.mixture_normal <- function(G_prior)
{
  out <- G_prior[[1]]$fitted_g$sd
  return(out)
}



#'@title Get mixture standard deviations mixture normal prior per scale
#'
#'@description
#'@param G_prior mixture normal prior
#'@return list of vectors of standard deviations

get_sd_G_prior.mixture_normal_per_scale <- function(G_prior)
{
  out <- lapply(G_prior, function(x) x$fitted_g$pi)
  return(out)
}
