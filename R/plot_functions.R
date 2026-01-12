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
#' @param which_plot Which plots to return; a PIP plot, effect plot,
#'   or both (in which case the return value is a list containing the
#'   two plots.
#'
#' @param point_size numeric, size of the points.
#'
#' @param pos_SNP vector, containing the base pair of the SNPs
#'
#' @param point_shape vector, containing the shape of dots
#'
#'
#' @param show_outing_grid logical, if TRUE show grid
#'
#'
#' @param show_affected_region logical, if TRUE show affected regions
#'
#' @param font_size Passed as the \dQuote{ont_size} argument to
#'   \code{\link[cowplot]{theme_cowplot}}.
#'
#' @param title The title of the plot.
#'
#' @param \dots additional arguments
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
                        which_plot = c("both","pip","effect"),
                        effect = "all",
                        cred_band = TRUE,
                        show_affected_region = TRUE,
                        show_outing_grid = ifelse(diff(range(diff(obj$outing_grid))) < 1e-6,FALSE,TRUE),
                        lfsr_curve = TRUE,
                        line_width = 0.35,
                        point_size = 1.25,
                        dot_size = 0.5,
                        pos_SNP,
                        point_shape,
                        font_size = 10,
                        title = "",
                         ...) {

  if (missing(pos_SNP)) {
    pos_SNP <- 1:length(obj$pip)
  }
  if (missing(point_shape)) {
    point_shape <- rep(19,length(pos_SNP))
  }
  which_plot <- match.arg(which_plot)

  p1 <- plot_susiF_pip(obj,title,pos_SNP,point_shape,point_size,font_size)

  p2 <- plot_susiF_effect(obj,effect,title,cred_band,show_affected_region,
                          show_outing_grid,lfsr_curve,line_width,dot_size,
                          font_size)

  if (which_plot == "both")
    return(list(pip = p1,effect = p2))
  else if (which_plot == "pip")
    return(p1)
  else
    return(p2)
}

#' @rdname fsusie_plots
#'
#' @param x Output of the susiF function.
#'
#' @importFrom graphics plot
#'
#' @method plot susiF
#'
#' @export
#'
plot.susiF <- function (x, ...)
  plot_susiF(x,...)

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
                            pos_SNP,
                            point_shape,
                            point_size = 1.25,
                            font_size = 10) {
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
#' @param show_affected_region If \code{show_affected_region = TRUE},
#'   the regions in which the credible bands cross zero are also shown.
#'
#' @param show_outing_grid If \code{show_outing_grid = TRUE}, show
#'   the grid positions at which the effects were estimated. By default,
#'   this option is set to \code{TRUE} only when the grid poisitions are
#'   uneven.
#'
#' @param lfsr_curve Logical. If \code{TRUE}, plot estimated lfsr of the
#'   effect at each base pair if obj fitted with HMM regression. This
#'   has no effect unless the \code{\link{susiF}} option
#'   \code{post_processing = "HMM"} was used.
#'
#' @param line_width Numeric. Width of the plotted lines.
#'
#' @param dot_size numeric, size of the points in the effect plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_point
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
                               show_outing_grid = ifelse(diff(range(diff(obj$outing_grid))) < 1e-6,FALSE,TRUE),
                               lfsr_curve = TRUE,
                               line_width = 0.35,
                               dot_size = 0.5,
                               font_size = 10) {
  # Declare variables to avoid R CMD check notes
  lwr <- upr <- Start <- End <- ystart <- yend <- NULL
  L     <- obj$L
  n_wac <- obj$n_wac
  y     <- obj$pip
  col_y <- rep(0,length(y))
  #this is to handle the case where no postprocessing is used
  if(sum( obj$cred_band[[1]]==0)== prod(dim( obj$cred_band[[1]]))){
    cred_band=FALSE
    show_affected_region =FALSE
  }
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
      geom_line(linewidth = line_width) +
      geom_ribbon(aes(ymin = lwr,ymax = upr,fill = CS,color = CS),
                  linewidth = 0,alpha = 0.3)
    if (show_outing_grid)
      out <- out + geom_point(size = dot_size)
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
        geom_line(linewidth = line_width) +
        geom_line(aes(y = lfsr_curve,x = x),color = "black",
                  show.legend = FALSE) +
        geom_hline(yintercept = 0.05,color = "black",linetype = "dashed")
      if (show_outing_grid)
        out <- out + geom_point(size = dot_size,shape = 21,fill = "white")
      facet_scale <- "free"
    } else {

      # Do not include the credible bands.
      x  <- rep(obj$outing_grid,length(indx_effect) + 1)
      CS <- rep(c(0,indx_effect),each = n_wac)
      df <- data.frame(fun_plot = fun_plot,CS = factor(CS),x = x)
      df <- df[-which(df$CS == 0),]
      out <- ggplot(df,aes(y = fun_plot,x = x,color = CS)) +
        geom_line(linewidth = line_width)
      if (show_outing_grid)
        out <- out + geom_point(size = dot_size,shape = 21,fill = "white")
      facet_scale <- "free"
    }
  }

  # Add a horizontal line.
  out <- out + geom_hline(yintercept = 0,color = "black",linetype = "dotted")

  # Add the "affected region".
  if (show_affected_region) {
    affected_region_dat <- cbind(affected_reg(obj),
                                 data.frame(ystart = 0,yend = 0))
    rows <- which(is.element(affected_region_dat$CS,indx_effect))
    affected_region_dat <- affected_region_dat[rows,]
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



#' @title fSuSiE Plots using Gviz
#'
#' @param obj is the fsusie object
# (i.e., the result of running fsusie)
#' @param chr the chromosome number
#' @param pos0 the start of the RNAseq count (genomic position of the first column of
#' the Y matrix used to fit the fsusie object) description
#' @param pos1 the end of the RNAseq count (genomic position of the last column of
#' the Y matrix used to fit the fsusie object) description
#'
#' @param X the X matrix used to fit fsusie
#' @param Y  the Y matrix used to fit fsusie
#' @param snp_info 'optional) a matrix containing the information of the genotype matrix see vignette on RNaseq
#' @param cs the cs number to be plotted
#' @param log1p_count logical (default to FALSE) show the observe count conditional on the leads SNP in the
#' using log1p_count
#' @param effect_log logical (set to TRUE) , the plot assume that you fitted the fsusie object
#' on log +1 count and so if you set this parameter to FALSE the displayed effect will be the expected
#' difference in count instead of the fitted curve
#' @param thresh_lfsr if the susiF object is fitted using HMM postprocessing you can use this argument
#' to set to 0 the estimated effect that have an local false sign rate higher than a given threshold (e.g., 0.05) description
#' @param  type_data  set to "p" change the type of plot for the observed count conditional of the lead SNP
#' @param data_splice position of splicing sites
#' see GViz documentation
# (scaled) read count matrices, respectively.
fsusie_log_plot <- function (obj, chr, pos0, pos1, X, Y, snp_info, cs = 1,
                             log1p_count=FALSE,
                             effect_log=TRUE,
                             thresh_lfsr=NULL,
                             data_splice=NULL,
                             type_data="p") {

  # Extract the relevant genes and exons in the specified region
  region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

  # Subset the genes and exons to the region of interest.
  region_genes <- subsetByOverlaps(region_genes,
                                   GRanges(seqnames = chr,
                                           ranges = IRanges(pos0,pos1)))

  # Generate a sequence of positions with a length of 1,024.
  positions <- seq(pos0,pos1,length.out = 1024)

  markers <- obj$cs[[cs]]
  j       <- which.max(obj$pip[markers])
  marker  <- markers[j]
  x       <- X[,marker]




  if(log1p_count){
    if (! length(which(x==2))>1){
      read_counts <- rbind(colMeans(log1p(Y[x == 0,])),
                           colMeans(log1p(Y[x == 1,])) )
    }else{
      read_counts <- rbind(colMeans(log1p(Y[x == 0,])),
                           colMeans(log1p(Y[x == 1,])),
                           colMeans(log1p(Y[x == 2,])))
    }

  }else{
    if (! length(which(x==2))>1){

      read_counts <- rbind(colMeans(Y[x == 0,]),
                           colMeans(Y[x == 1,]) )
    }else{
      read_counts <- rbind(colMeans(Y[x == 0,]),
                           colMeans(Y[x == 1,]),
                           colMeans(Y[x == 2,]))
    }

  }

  if(effect_log){
    effect=obj$fitted_func[[cs]]

    if(!is.null(thresh_lfsr)){
      effect=(obj$fitted_func[[cs]]) * ifelse(obj$lfsr_func[[cs]]< thresh_lfsr,1,0 )


    }


  }else{



    Y_mean= colMeans(log1p(Y[x == 0,]))
    effect <-  exp(Y_mean)* exp(obj$fitted_func[[cs]] )-exp(Y_mean)
    if(!is.null(thresh_lfsr)){
      effect=  exp(Y_mean)* exp(obj$fitted_func[[cs]] * ifelse(obj$lfsr_func[[cs]]< thresh_lfsr,1,0 ))-exp(Y_mean)
    }



  }


  # Create a "data track" to show the CS effect.
  cex <- 0.6

  effect_track <-
    DataTrack(range = GRanges(seqnames = chr,
                              ranges = IRanges(start = positions,
                                               end = positions + 1)),
              data = effect, genome = "hg38",
              name = paste("CS",cs),type = "l",col = "royalblue",
              track.margin = 0.05,cex.title = cex,cex.axis = cex,
              col.axis = "black",col.title = "black",
              fontface = "plain",background.title = "white",
              fontface.title = 1)

  # Create another "data track" to show the read counts.

  n0  <- sum(x == 0)
  n1  <- sum(x == 1)
  n2  <- sum(x == 2)
  id  <- snp_info[marker,"ID"]
  ref <- snp_info[marker,"REF"]
  alt <- snp_info[marker,"ALT"]

  if (! length(which(x==2))>1){
    groups <- c(sprintf("%s %s%s (n = %d)",id,ref,ref,n0),
                sprintf("%s %s%s (n = %d)",id,ref,alt,n1) )
    geno_colors <- c("navyblue","turquoise" )
  }else{
    groups <- c(sprintf("%s %s%s (n = %d)",id,ref,ref,n0),
                sprintf("%s %s%s (n = %d)",id,ref,alt,n1),
                sprintf("%s %s%s (n = %d)",id,alt,alt,n2))
    geno_colors <- c("navyblue","turquoise","darkorange")
  }


  if (mean(effect) > 0) {
    groups <- factor(groups,rev(groups))
    geno_colors <- rev(geno_colors)
  } else {
    groups <- factor(groups,groups)
  }


  lab_y =ifelse(log1p_count, "avg. log1p count","avg. count")
  data_track <- DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions + 1)),
                          data = read_counts,genome = "hg38",
                          groups = groups,
                          name = lab_y  , type = type_data, #"p",#type = "l",
                          col = geno_colors  ,
                          track.margin = 0.05,cex.title = cex,cex.axis = cex,
                          col.axis = "black",col.title = "black",
                          fontface = "plain",background.title = "white",
                          fontface.title = 1,cex.legend = cex, cex=0.2)

  # Create an "ideogram" track.
  ideo_track <- IdeogramTrack(genome = "hg38",chromosome = chr)

  # Create a "genome axis" track.
  genome_track <- GenomeAxisTrack(col.axis = "black",col.title = "black")

  # Create a "gene region" track.
  gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                                pos0 = pos0,pos1 = pos1,name = "",
                                showId = TRUE,geneSymbol = TRUE,
                                col.axis = "black",col.title = "black",
                                transcriptAnnotation = "symbol",
                                rotation.title = 0,cex.title = cex,
                                col = "salmon",fill = "salmon",
                                background.title = "white")

  # Map gene IDs to gene symbols.
  gene_ids <- unique(unlist(region_genes$gene_id))

  # Map to gene symbols using org.Hs.eg.db
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,keys = gene_ids,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID")
  n <- nrow(gene_symbols)
  if (n > 0) {
    for (i in 1:n) {
      j <- which(gene_track@range@elementMetadata@listData$gene ==
                   gene_symbols$ENTREZID[i])
      gene_track@range@elementMetadata@listData$id[j]     <- gene_symbols$SYMBOL[i]
      gene_track@range@elementMetadata@listData$symbol[j] <- gene_symbols$SYMBOL[i]
    }
  }


  if( !is.null(data_splice)){
    junction_ranges <- GRanges(
      seqnames = chr,
      ranges = IRanges(
        start = as.numeric(sub(".*_(\\d+)_.*", "\\1", data_splice$Name)),
        end = as.numeric(sub(".*_(\\d+)$", "\\1", data_splice$Name))
      ),
      names = data_splice$Description
    )

    # Offset overlapping junctions for better visualization
    # Create AnnotationTrack for splicing junctions with adjusted thickness
    junction_track <- AnnotationTrack(
      range = junction_ranges,
      genome = "hg38",
      chromosome = chr,
      name = "Splicing Junctions",
      stacking = "squish",
      col = "darkgreen",
      fill = "lightgreen",
      track.margin = 0.05,
      cex.title = cex,
      cex.axis = cex,
      col.axis = "black",
      col.title = "black",
      fontface = "plain",
      background.title = "white",
      height = 0.3  # Adjust height here to make rectangles thinner
    )
    # Combine all tracks into a single plot
    tracks <- c(
      ideo_track,
      genome_track,
      effect_track,
      data_track,
      gene_track,
      junction_track
    )

    # Plot the tracks
    return(plotTracks(tracks, from = pos0, to = pos1, sizes = c(1, 1.75, 2, 4, 5, 1)))
  }else{
    # Combine all tracks into a single plot.
    tracks <- c(ideo_track,
                genome_track,
                effect_track,
                data_track,
                gene_track)
    return(plotTracks(tracks,from = pos0,to = pos1,sizes = c(1,1.75,2,4,5)))
  }


}

