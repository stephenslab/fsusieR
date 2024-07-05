plot_colors <- c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                 "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
                 "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                 "deeppink1", "blue1", "steelblue4", "darkturquoise",
                 "green1", "yellow4", "yellow3", "darkorange4", "brown")

#' @rdname fsusie_plots
#'
#' @title fSuSiE Plots
#' 
#' @description Various visualizations of fSuSiE results.
#
#' @param obj Output of the susiF function.
#' 
#' @param point_size numeric, size of the point
#' 
#' @param pos_SNP vector, containing the base pair of the SNPs
#' 
#' @param point_shape vector, containing the shape of dots
#' 
#' @param pip_only logical, if TRUE only ouput the PIP plot
#' 
#' @param font_size Passed as the \dQuote{ont_size} argument to
#'   \code{\link[cowplot]{theme_cowplot}}.
#' 
#' @param title The title of the plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 scale_fill_manual
#' 
#' @export
#
plot_susiF <- function (obj,
                        title = "",
                        effect = "all",
                        cred_band = TRUE,
                        lfsr_curve = TRUE,
                        linewidth = 0.5,
                        point_size = 2,
                        pos_SNP,
                        point_shape,
                        pip_only = FALSE) {
  
  if(missing(pos_SNP)){
    pos_SNP<-  1:length(obj$pip)
  }
  if( missing(point_shape)){
    point_shape <- rep( 19, length(pos_SNP))
  }
  
  P1  <- plot_susiF_pip  (obj         = obj  ,
                           point_size  = point_size,
                           pos_SNP     = pos_SNP,, 
                           point_shape = point_shape)
  
  P2 <- plot_susiF_effect(obj = obj ,
                          effect = effect,
                          cred_band = cred_band,
                          lfsr_curve = lfsr_curve,
                          linewidth = linewidth)
  
  if (pip_only) {
    return(P1)
  }
  return(out <- gridExtra::grid.arrange(P1,P2,ncol=2,top =title))
}

#' @rdname fsusie_plots
#'
#' @export
#'
plot.susiF <- function(x, ...) {
  return(plot_susiF(x,...))
}
 
#' @rdname fsusie_plots
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab 
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
plot_susiF_pip <- function (obj,
                            title = "", 
                            point_size = 2,
                            pos_SNP, 
                            point_shape,
                            font_size = 10 ) {
  if (missing(pos_SNP)) {
    pos_SNP <- 1:length(obj$pip)
  }
  if (missing(point_shape)) {
    point_shape <- rep(19,length(pos_SNP))
  }
  L <- obj$L
  y <- obj$pip
  col_y <- rep(0,length(y))
  for (l in 1:L) {
    col_y[which(1:length(y) %in% obj$cs[[l]])] <- l
  }
  CS <- factor(col_y,
               levels = 0:L,
               labels = c("none",1:L))
  df <- data.frame(y = y,CS = CS,pos_SNP)
  return(ggplot(df,aes(y = y,x = pos_SNP,color = CS)) +
         geom_point(size = point_size,
                    shape = point_shape) +
         scale_color_manual("CS",values = plot_colors) +
         labs(x = "SNP",y = "PIP",title = title) +
         theme_cowplot(font_size = font_size))
}

#' @rdname fsusie_plots
#' 
#' @param effect The indices of the effects to be plotted, or use
#'   \code{effect = "all"} to plot all effects.
#' 
#' @param cred_band logical. If \code{TRUE}, plot credible bands if
#'   the fSuSiE model was fitted with wavelet regression.
#'
#' @param show_affected_region. If \code{show_affected_region = TRUE},
#'   the regions in which the credible bands cross zero are also shown.
#' 
#' @param lfsr_curve Logical. If \code{TRUE}, plot estimated lfsr of the
#'   effect at each base pair if obj fitted with HMM regression. This
#'   has no effect unless the \code{\link{susiF}} option
#'   \code{post_processing = "HMM"} was used.
#' 
#' @param linewidth Numeric. Width of the plotted lines.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#
plot_susiF_effect <- function (obj,
                               effect = "all",
                               title = "",
                               cred_band = TRUE,
                               show_affected_region = TRUE,
                               lfsr_curve = TRUE,
                               linewidth = 0.5,
                               font_size = 10) {
  L     <- obj$L
  n_wac <- obj$n_wac
  y     <- obj$pip
  col_y <- rep(0,length(y))
 
  if (is.character(effect)) {
    if (effect == "all") {
      indx_effect <- 1:obj$L
    } else {
      stop(paste("Entry not supported"))
    }
  } else{
    if (is.numeric(effect) == 1 & length(which(effect > L)) > 0) {
      stop(paste("The specified effect should be lower or equal to",L))
    }
    indx_effect <- effect
    indx_effect <- indx_effect[order(indx_effect)]
  }

  L        <- length(indx_effect)
  fun_plot <- do.call(c,obj$fitted_func[indx_effect])
  n_wac    <- obj$n_wac
  fun_plot <- c(rep(0,n_wac),fun_plot)

  if (cred_band & is.null(obj$lfsr_func)) {

    # Include the credible bands.
    cred_band_dat <- data.frame(t(do.call(cbind,obj$cred_band[indx_effect])))
    cred_band_dat <- rbind(data.frame(up = rep(0,n_wac),low = rep(0,n_wac)),
                           cred_band_dat)
    x  <- rep(obj$outing_grid,length(indx_effect) + 1)
    CS <- rep(c(0,indx_effect),each = n_wac)
    df <- data.frame(fun_plot = fun_plot,CS = factor(CS),x = x,
                     upr = cred_band_dat$up,lwr = cred_band_dat$low)
    df  <- df[-which(df$CS == 0),]
    out <- ggplot(df,aes(y = fun_plot,x = x,color = CS)) +
      geom_line(linewidth = linewidth) +
      geom_ribbon(aes(ymin = lwr,ymax = upr,fill = CS,color = CS),
                  linewidth = 0,alpha = 0.3)
    facet_scale <- "fixed"
  } else {
    if (lfsr_curve & !is.null(obj$lfsr_func)) {

      # Add lfsr information.
      lfsr_curve_dat <- do.call(c,obj$lfsr_func[indx_effect])
      x  <- rep(obj$outing_grid,length(indx_effect) + 1)
      CS <- rep(c(0,indx_effect), each = n_wac)
      df <- data.frame(fun_plot = fun_plot,CS = factor(CS),x = x)
      df <- df[-which(df$CS == 0),]
      df$lfsr_curve <- lfsr_curve_dat
      out <- ggplot(df,aes(y = fun_plot,x = x,color = CS)) +
        geom_line(linewidth = linewidth) +
        geom_line(aes(y = lfsr_curve,x = x),color = "black",
                  show.legend = FALSE) +
        geom_hline(yintercept = 0.05,color = "black",linetype = "dashed")
      facet_scale <- "free"
    } else {
        
      # Do not include the credible bands.
      x  <- rep(obj$outing_grid,length(indx_effect) + 1)
      CS <- rep(c(0,indx_effect),each = n_wac)
      df <- data.frame(fun_plot = fun_plot,CS = factor(CS),x = x)
      df <- df[-which(df$CS == 0),]
      out <- ggplot(df,aes(y = fun_plot,x = x,color = CS)) +
        geom_line(linewidth = linewidth)
      facet_scale <- "free"
    }
  }

  # Add a horizontal line.
  out <- out + geom_hline(yintercept = 0,color = "black",linetype = "dotted")

  # Add the "affected region".
  if (show_affected_region) {
    affected_region_dat <- cbind(affected_reg(obj),
                                 data.frame(ystart = 0,yend = 0))
    affected_region_dat$CS <- factor(affected_region_dat$CS)
    out <- out + geom_segment(aes(x = Start,xend = End,y = ystart,yend = yend),
                              data = affected_region_dat,linewidth = 0.75,
                              lineend = "square",color = "black")
  }

  # Finish up the plot.
  plot_colors <- plot_colors[-1][indx_effect]
  return(out +
         facet_grid(CS~.,scales = facet_scale) +
         scale_color_manual("Credible set",values = plot_colors) +
         scale_fill_manual("Credible set",values = plot_colors) +
         labs(x = "position",y = "estimated effect",title = title) +
         theme_cowplot(font_size = font_size))
}
