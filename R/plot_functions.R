
#' @title Plot posterior inclusion probabilities from susiF object
#' 
#' @description Plot posterior inclusion probabilities from susiF object
#' 
#' 
#' @param point_size numeric, size of the point
#' 
#' @param pos_SNP vector, containing the base pair of the SNPs
#' 
#' @param point_shape vector, containing the shape of dots
#' 
#' @param font_size numeric, the size of the font
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab 
#' @importFrom ggplot2 aes
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 labs
#' @export
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
  color <- c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
             "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
             "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
             "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
             "yellow3", "darkorange4", "brown")
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
  return(ggplot(df, aes_string(y = "y",x = "pos_SNP",color = "CS")) +
         geom_point(size = point_size,
                    shape = point_shape) +
         scale_color_manual("CS",values = color) +
         labs(x = "SNP",y = "PIP",title = title) +
         theme_cowplot(font_size = font_size))
 
  
}

#' @title Plot specific effect from susiF object
#' 
#' @description Plot specific effect from susiF object
#' 
#' @param obj Output of the susiF function.
#' 
#' @param effect The indices of the effects to be plotted, or use
#'   "all" to plot all effects.
#' 
#' @param cred.band logical, if TRUE, plot credible bands if obj
#'   fitted with wavelets regression.
#' 
#' @param lfsr.curve logical, if TRUE, plot estimated lfsr of the
#'   effect at each base pair if obj fitted with HMM regression.
#' 
#' @param linewidth numeric, width of the plotted lines
#'
#' @title font_size Passed as the \dQuote{ont_size} argument to
#'   \code{\link[cowplot]{theme_cowplot}}.
#' 
#' @param title The title of the plot.
#'  
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
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
#' @importFrom ggplot2 aes
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 labs
#' @export
#
plot_susiF_effect  <- function (obj,
                                effect = "all",
                                title = "",
                                cred.band = TRUE,
                                lfsr.curve = TRUE,
                                linewidth = 0.5,
                                font_size = 10 ) {
  color <- c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
             "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
             "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
             "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
             "yellow3", "darkorange4", "brown")
  L     <- obj$L
  n_wac <- obj$n_wac
  y     <- obj$pip
  col_y <- rep(0,length(y))
 
  if (is.character(effect)) {
    if (effect == "all") {
      indx_effect <- 1:obj$L
    } else {
      stop(paste("entry not supported"))
    }
  } else{
    if (is.numeric(effect) == 1 & length(which(effect > L)) > 0)  {
      stop(paste("the specified effect should be lower or equal to ",L))
    }
    indx_effect <- effect
    indx_effect <- indx_effect[order(indx_effect)]
  }

  L <- length(indx_effect)
  fun_plot <- do.call(c, obj$fitted_func[indx_effect])
  n_wac    <- obj$n_wac
  fun_plot <- c(rep(0, n_wac), fun_plot)

  if (cred.band & is.null(obj$lfsr_func)) {
    cred_band <- data.frame(t(do.call(cbind,
                                      obj$cred_band[indx_effect])))
    cred_band <- rbind(data.frame(up = rep(0, n_wac),
                                  low = rep(0, n_wac)),
                       cred_band)
    
    x  <- rep(obj$outing_grid,length(indx_effect) + 1)
    CS <- rep(c(0,indx_effect), each = n_wac)
    df <- data.frame(fun_plot = fun_plot,CS = as.factor(CS),x = x,
                     upr = cred_band$up,
                     lwr = cred_band$low)
    df  <- df[-which(df$CS == 0),]
    out <- ggplot(df, aes(y = fun_plot, x = x, color = as.factor(CS))) +
           geom_line(linewidth = linewidth) +
           geom_ribbon(aes(ymin = lwr, ymax = upr, fill = as.factor(CS), color = as.factor(CS)), linewidth = 0, alpha = 0.3) +
           scale_color_manual("Credible set", values = color[-1][indx_effect]) +
           scale_fill_manual("Credible set", values = color[-1][indx_effect]) +
           facet_grid(as.factor(CS) ~ .)
    
  } else {
    if (lfsr.curve & !is.null(obj$lfsr_func)){
      lfsr_curve <- do.call(c,obj$lfsr_func[indx_effect])
      x  <- rep(obj$outing_grid,length(indx_effect) + 1)
      CS <- rep(c(0,indx_effect), each = n_wac)
      df <- data.frame(fun_plot = fun_plot,
                       CS = as.factor(CS),
                       x = x)
      df <- df[-which(df$CS==0),]
      df$lfsr_curve <- lfsr_curve
      out <-   ggplot(df, aes(y = fun_plot, x = x, color = as.factor(CS))) +
               geom_line(linewidth = linewidth, show.legend = TRUE) +
               geom_line(aes(y = lfsr_curve, x = x), color = "black", show.legend = FALSE) +
               geom_hline(yintercept = 0.05, color = "black", linetype = "dashed") +
               geom_hline(yintercept = 0.0, color = "black") +
               scale_color_manual("Credible set", values = color[-1][indx_effect]) +
               facet_grid(CS ~ ., scales = "free")
      

    } else{
      x  <- rep(obj$outing_grid,length(indx_effect) + 1)
      CS <- rep(c(0,indx_effect), each = n_wac)
      df <- data.frame(fun_plot = fun_plot,
                       CS = as.factor(CS),
                       x = x)
      df <- df[-which(df$CS==0),]
      out <- ggplot(df, aes(y = fun_plot, x = x, color = as.factor(CS))) +
             geom_line(linewidth = linewidth, show.legend = TRUE) +
             scale_color_manual("Credible set", values = color[-1][indx_effect]) +
             facet_grid(CS ~ ., scales = "free")
    }
  }
  out <- out +
    labs(x = "position",y = "estimated effect",title = title) +
    theme_cowplot(font_size = font_size)
  return(out)
}

#' @title Plot susiF object
#
#' @param obj output of the susiF function
#' 
#' @param effect  the index of the effect to be plotted or use "all" to plot all effects
#' 
#' @param cred.band logical, if TRUE, plot credible bands if obj fitted with wavelets regression. Set as TRUE by default
#' 
#' @param  lfsr.curve logical, if TRUE, plot estimated lfsr of the effect at each base pair  if obj fitted with HMM regression. Set as TRUE by default
#' 
#' @param linewidth numeric, width of the plotted lines
#' 
#' @param point_size numeric, size of the point
#' 
#' @param pos_SNP vector, containing the base pair of the SNPs
#' 
#' @param point_shape vector, containing the shape of dots
#' 
#' @param pip_only logical, if TRUE only ouput the PIP plot
#' 
#' @param title character 
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
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
#' @importFrom ggplot2 aes 
#' @export
#
plot_susiF<- function (obj,
                        title="",
                        effect= "all",
                        cred.band = TRUE,
                        lfsr.curve=TRUE,
                        linewidth=1.5,
                        point_size=2,
                        pos_SNP,
                        point_shape,
                        pip_only =FALSE)
{
  
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
  
  P2 <- plot_susiF_effect ( obj = obj ,
                      effect= effect,
                      cred.band = cred.band,
                      lfsr.curve= lfsr.curve,
                      linewidth=linewidth  )
  
  if(pip_only){
   
    return(P1)
  }
  
  return(out <- gridExtra::grid.arrange(P1,P2,ncol=2,top =title))
}



#'
#' @export
#'
plot.susiF <- function(x,
                      ... ) {
  
                      return(plot_susiF (obj=x  ))
  }

 
