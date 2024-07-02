

plot_susiF_pips <- function (obj, title="", 
                             size_point=4,
                             pos_SNP, 
                             point_shape, ...){
  
  if(missing(pos_SNP)){
    pos_SNP<-  1:length(obj$pip)
  }
  if( missing(point_shape)){
    point_shape <- rep( 19, length(pos_SNP))
  }
  
  
  color = c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
            "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
            "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
            "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
            "yellow3", "darkorange4", "brown")
  L <- obj$L
  
  y <- obj$pip
  col_y <- rep(0, length(y))
  
  for (l in 1:L) {
    col_y[which(1:length(y) %in% obj$cs[[l]])] <- l
  }
  CS <-  factor(col_y,
                levels=0:L,
                labels = c("Not in CS", 1:L))
  
  
  df <-  data.frame(y = y, CS = CS)
  
  
  P1 <- ggplot(df, aes_string(y = "y",
                              x ="pos_SNP",
                              col = "CS")) +
    geom_point(size = size_point,
               shape=point_shape) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_color_manual("Credible set",
                       values = color) +
    xlab("SNP index") +
    ylab("Posterior Inclusion Probability (PIP)")
  
  return(P1)
  
}





#' @title Plot specific effect from susiF object
#' @description  Plot specific effect from susiF object
#' @param obj output of the susiF function
#' @param effect  the index of the effect to be plotted or use "all" to plot all effects
#' @param cred.band logical, if TRUE, plot credible bands if obj fitted with wavelets regression. Set as TRUE by default
#' @param  lfsr.curve logical, if TRUE, plot estimated lfsr of the effect at each base pair  if obj fitted with HMM regression. Set as TRUE by default
#' @param size_line numeric, width of the plotted lines 
#' @param title character
#' @param \dots Other arguments..
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
#
#' @export
#
plot_effect  <- function( obj,
                          effect=1,
                          title="",
                          cred.band = TRUE,
                          lfsr.curve=TRUE,
                          size_line=2, ...){
  
  
  
  
  
   
  color = c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
            "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
            "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
            "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
            "yellow3", "darkorange4", "brown")
  L <- obj$L
  n_wac <- obj$n_wac
  y <- obj$pip
  col_y <- rep(0, length(y))
  
  
  
  if ( is.character(effect)   ){
    if(effect =="all"){
      indx_effect <- 1:obj$L
    }else{
      stop(paste("entry not supported" ))
    }
    
    
    
  }else{
    if (is.numeric(effect)==1 & length(which(effect>L))>0)  {
      stop(paste("the specified effect should be lower or equal to ",
                 L))
    }
    indx_effect <- effect
  }
  
    

 
  
   
  L= length(indx_effect)
  fun_plot <- do.call(c, obj$fitted_func[indx_effect])
  n_wac <- obj$n_wac
  fun_plot <- c(rep(0, n_wac), fun_plot)
  
  
  if (cred.band& is.null(obj$lfsr_func)) {
    cred_band <- data.frame(t(do.call(cbind,
                                      obj$cred_band[indx_effect] )
    )
    )
    cred_band <- rbind(data.frame(up = rep(0, n_wac),
                                  low = rep(0, n_wac)),
                       cred_band)
    
    x  <- rep(obj$outing_grid, (length(indx_effect) + 1))
    CS <- rep(0:L, each = n_wac)
    df <- data.frame(fun_plot = fun_plot,
                     CS = as.factor(CS),
                     x = x,
                     upr = cred_band$up,
                     lwr = cred_band$low)
    df <- df[-which(df$CS==0),]
    P2 <- ggplot(df, aes_string(y = "fun_plot",
                                x = "x",
                                col = "CS")) +
      geom_line(linewidth = size_line) +
      geom_ribbon(aes_string(ymin = "lwr",ymax = "upr",fill = "CS",
                             col = "CS"),alpha = 0.3) +
      scale_color_manual("Credible set", values = color[-1][indx_effect]) +
      scale_fill_manual("Credible set", values = color[-1][indx_effect]) +
      facet_grid(CS~.) +
      xlab("postion") + ylab("Estimated effect")+
      theme(legend.position = "none", 
            strip.background = element_blank(), 
            strip.text = element_blank())
     
  }else {
    
    
    if (lfsr.curve& !is.null(obj$lfsr_func)){
      lfsr_curve <- do.call(c , obj$lfsr_func)
      
      
      
      x <- rep(obj$outing_grid, (indx_effect + 1))
      CS <- rep(0:L, each = n_wac)
      df <- data.frame(fun_plot = fun_plot,
                       CS = as.factor(CS),
                       x = x
      )
      df <- df[-which(df$CS==0),]
      df$   lfsr_curve <-    lfsr_curve
      P2 <- ggplot(df, aes_string(y = "fun_plot",
                                  x = "x",
                                  col = "CS")) +
        geom_line(linewidth = size_line) +
        geom_line(aes(y=lfsr_curve, x=x, col="black"))+
        geom_hline(yintercept = 0.05)+
        scale_color_manual("Credible set",
                           values = color[-1][indx_effect]) +
        facet_grid(CS~., scales = "free")+
        xlab("postion") + ylab("Estimated effect")+
        theme(legend.position = "none", 
              strip.background = element_blank(), 
              strip.text = element_blank())
    }else{
      
      
      x <- rep(obj$outing_grid, (indx_effect + 1))
      CS <- rep(0:L, each = n_wac)
      df <- data.frame(fun_plot = fun_plot,
                       CS = as.factor(CS),
                       x = x)
      df <- df[-which(df$CS==0),]
      P2 <- ggplot(df, aes_string(y = "fun_plot",
                                  x = "x",
                                  col = "CS")) +
        geom_line(linewidth = size_line) + scale_color_manual("Credible set",
                                                              values = color[-1][indx_effect]) +
        geom_hline(yintercept=0, linetype='dashed', col = 'grey', linewidth = 1.5)+
        facet_grid(CS~., scales = "free")+
        xlab("postion") + ylab("Estimated effect")+
        theme(legend.position = "none", 
              strip.background = element_blank(), 
              strip.text = element_blank())
    }
    
    
    
  }
  return(P2)
  
}




#' @title Plot susiF object
#
#' @param obj output of the susiF function
#' @param effect  the index of the effect to be plotted or use "all" to plot all effects
#' @param cred.band logical, if TRUE, plot credible bands if obj fitted with wavelets regression. Set as TRUE by default
#' @param  lfsr.curve logical, if TRUE, plot estimated lfsr of the effect at each base pair  if obj fitted with HMM regression. Set as TRUE by default
#' @param size_line numeric, width of the plotted lines
#' @param size_point numeric, size of the point
#' @param pos_SNP vector, containing the base pair of the SNPs
#' @param point_shape vector, containing the shape of dots
#' @param pip_only logical, if TRUE only ouput the PIP plot
#' @param title character
#' @param \dots Other arguments..
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
#
#' @export
#
plot_susiF  = function (obj, title="",
                        effect= "all",
                        cred.band = TRUE,
                        lfsr.curve=TRUE,
                        size_line=2,
                        size_point=4,
                        pos_SNP,
                        pip_only=FALSE,
                        point_shape, ...)
{
  
  if(missing(pos_SNP)){
    pos_SNP<-  1:length(obj$pip)
  }
  if( missing(point_shape)){
    point_shape <- rep( 19, length(pos_SNP))
  }
  
  
  P1  <- plot_susiF_pips  (obj         = obj  ,
                           size_point  = size_point,
                           pos_SNP     = pos_SNP,, 
                           point_shape = point_shape)
  
  P2 <- plot_effect ( obj = obj ,
                      effect= effect,
                      cred.band = cred.band,
                      lfsr.curve= lfsr.curve,
                      size_line=size_line,
                      size_point= size_point,
                      pos_SNP = pos_SNP )
  
  if(pip_only){
   
    return(P1)
  }
  
   
  
  return(out <- gridExtra::grid.arrange(P1,P2,ncol=2,top =title))
}





 
