## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#",collapse = TRUE,results = "hold",
                      fig.align = "center",dpi = 120)

## -----------------------------------------------------------------------------
library(fsusieR)
data("RNA_seq_example")
Y=RNAseqData$Y
X=RNAseqData$X
base_pair=RNAseqData$pos

## -----------------------------------------------------------------------------
print(dim(Y))
plot( base_pair, Y[1,], main=" Observed base pair count for individual 1")

## -----------------------------------------------------------------------------
Y_binned= bin_Y(Y,base_pair = base_pair)
Y=Y_binned$binned_data
print( dim( Y))
pos = Y_binned$pos
print(pos[1:10])
bin_size=Y_binned$bin_size
print( bin_size)

## -----------------------------------------------------------------------------
plot(pos, Y[1,], main=" Observed binned count for individual 1")

## -----------------------------------------------------------------------------
size_factor=rowSums(Y)/ mean(rowSums(Y))
Y_cor =   (  Y /as.vector(size_factor))

plot(pos, log1p(Y[1,] ))

## -----------------------------------------------------------------------------
uni_res= univariate_functional_regression(Y=log1p(Y_cor),
                                          X=as.matrix(X[,261],
                                                      ncol=1),
                                          method="smash",
                                          alpha = 0.1
                                          )

## -----------------------------------------------------------------------------
plot(pos, uni_res$effect_estimate, 
     ylab = "estimated effect",
     xlab="base pair",
     type="l",
     lwd=2,
     col="royalblue",
     ylim=c (min(uni_res$cred_band),max(uni_res$cred_band)))
lines( pos,uni_res$cred_band[1,], 
       lty=2,
       lwd=1.4,
       col="royalblue")
lines(pos, uni_res$cred_band[2,], 
        lty=2,
       lwd=1.4,
       col="royalblue")
abline(h=0)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
pos0 = pos[1]
pos1= pos[length(pos)]+bin_size
positions=pos
chr=RNAseqData$chr

## -----------------------------------------------------------------------------
effect=effect=rbind(uni_res$effect_estimate,
             uni_res$cred_band,
             rep(0,length(uni_res$effect_estimate)))
group_colors <- c("black" ,"royalblue","royalblue","royalblue" )
group_lwd= c(1,2,1,1)
x=X[,261]
read_counts <-  rbind( colMeans(log1p(Y_cor[x == 1,])),
                       colMeans(log1p(Y_cor[x == 0,])) 
                           )
obs_effect=rbind (read_counts[2,]-read_counts[1,],
                      rep(0, length(read_counts[2,] )) ) 

ylim= c( min( c(effect,obs_effect)), max(c(  effect,obs_effect)))

group_cred= c(1:3,0)  
effect_track <-
    DataTrack(range = GRanges(seqnames = chr,
                              ranges = IRanges(start = positions,
                                                 end = positions + 1)),
              ylim=ylim,
              data = effect,
              genome = "hg38",
               groups= group_cred,
              name = "Effect SNP 261",
              type = "l",   
              lwd= group_lwd,  
              col = group_colors,
              track.margin = 0.05,
              cex.title = 0.6,
              cex.axis = 0.6,
              col.axis = "black",
              col.title = "black",
              fontface = "plain",
              background.title = "white",
              fontface.title = 1, 
              legend=FALSE)

## -----------------------------------------------------------------------------

  # Create a "gene region" track.
  gene_track <- GeneRegionTrack(txdb,
                                genome = "hg38",
                                chromosome = chr,
                                pos0 = pos0,
                                pos1 = pos1,
                                name = "",
                                showId = TRUE,
                                geneSymbol = TRUE,
                                col.axis = "black",
                                col.title = "black",
                                transcriptAnnotation = "symbol",
                                rotation.title = 0,
                                cex.title = 0.6,
                                col = "salmon",
                                fill = "salmon",
                                background.title = "white")

## -----------------------------------------------------------------------------
x       <- X[,261]

read_counts <-  rbind( colMeans(log1p(Y_cor[x == 1,])),
                       colMeans(log1p(Y_cor[x == 0,])) 
                           )
 
groups <- c("SNP 261 = 0",
                "SNP 261 = 1" )
geno_colors <- c("turquoise", "navyblue" )
lab_y = "avg. log1p count" 
  data_track <- DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions + 1)),
                          data = read_counts,
                          genome = "hg38",
                          groups = groups,
                          name = lab_y  , 
                          type = "p" , #"p",#type = "l",
                          col = geno_colors  ,
                          track.margin = 0.05,
                          cex.title = 0.6,
                          cex.axis = 0.6,
                          col.axis = "black",
                          col.title = "black",
                          fontface = "plain",
                          background.title = "white",
                          fontface.title = 1,
                          cex.legend = 0.6,
                          cex=0.2)

## -----------------------------------------------------------------------------
obs_diff_track <-
    DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions + 1)),
                          data =  read_counts[1,]- read_counts[2,],
              genome = "hg38",
              ylim=ylim,             
              name = "observed difference \n between genotype",
              type = "p", 
              lwd= 2,
              col = "blue",
              track.margin = 0.05,
              cex.title = 0.6,
              cex.axis = 0.6,
              col.axis = "black",
              col.title = "black",
              fontface = "plain",
              background.title = "white",
              fontface.title = 1, 
              legend=FALSE)
obs_diff_track2 <-
    DataTrack(range = GRanges(seqnames = chr,
                              ranges = IRanges(start = positions,
                                               end = positions + 1)),
              data =  rep(0, 1024), 
              genome = "hg38",
              ylim= ylim, 
              name ="Observed difference ",
              type = c("l"  ),
              col = c("black" ),
              track.margin = 0.05,
              cex.title =  0.6,
              cex.axis =  0.6,
              col.axis = "black",
              col.title = "black",
              fontface = "plain",
              background.title = "white",
              fontface.title = 1,
              legend = FALSE )
  obs_diff_track =OverlayTrack(trackList = list(obs_diff_track,
                                                  obs_diff_track2),background.title = "white")
  
  

## -----------------------------------------------------------------------------
 tracks <- c( 
      effect_track,
      obs_diff_track,
      data_track,
      gene_track 
    )

plotTracks(tracks, from = pos0, to = pos1, sizes = c(  3,3, 3, 1 ))

## ----running_fsusie-----------------------------------------------------------
library(fsusieR)
result= susiF(Y=log1p(Y_cor),X=as.matrix(X), L=20 , post_processing = "smash"             )
#result= RNAseqData$results

## -----------------------------------------------------------------------------
result$cs

## -----------------------------------------------------------------------------
plot_susiF(result
           )

## -----------------------------------------------------------------------------
res_HMM = change_fit   ( result,
                      log1p(Y_cor),
                      X, 
                      to="HMM"
)

## -----------------------------------------------------------------------------
cs=2 
pos=  result$outing_grid
plot(x= pos ,
     y= result$fitted_func[[cs
                         ]], 
     ylab = paste("estimated effect for CS  ",cs),
     xlab="base pair",
     type="l",
     lwd=2,
     col="royalblue",
     ylim=c (min(result$cred_band[[cs]]),max(result$cred_band[[cs]])))
lines( x= pos ,
       y= result$cred_band[[cs]][1,], 
       lty=2,
       lwd=1.4,
       col="royalblue")
lines(  x= pos ,
        y= result$cred_band[[cs]][2,], 
        lty=2,
       lwd=1.4,
       col="royalblue")
abline(h=0)

## -----------------------------------------------------------------------------
tt= do.call( rbind , result$lBF)

plot(apply(tt,2,max), col =ifelse(result$pip>0.2,2,3), pch=20)


