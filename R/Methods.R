

log_BF    <- function( G_prior , ...) UseMethod("get_log_BF")
post_mat_mean <- function( G_prior , ...) UseMethod("post_mat_mean")
post_mat_sd   <- function( G_prior , ...) UseMethod("post_mat_sd")
EM_pi         <- function( G_prior , ...) UseMethod("EM_pi")
get_pi_G_prior <- function(G_prior, ...) UseMethod("get_pi_G_prior")
get_sd_G_prior <- function(G_prior, ...) UseMethod("get_sd_G_prior")
L_mixsq  <- function(G_prior, ...) UseMethod("L_mixsq")
