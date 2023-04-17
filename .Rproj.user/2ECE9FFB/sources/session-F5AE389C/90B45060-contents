library(susiF.alpha)
library(sim1000G)


set.seed(1)

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 200000 , min_maf = 0.02 , max_maf = NA )

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4)

startSimulation( vcf )



'%!in%' <- function(x,y)!('%in%'(x,y))
id = c()
for(i in 1:100) id[i] = SIM$addUnrelatedIndividual()





genotypes = SIM$gt1[1:100,] + SIM$gt2[1:100,]


print(dim(genotypes))
N <- nrow(genotypes)



L <- sample(1:20, size=1)#actual number of effect
lf <-  list()
for(l in 1:L){
  lf[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
}


G <- genotypes

if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}

true_pos <- sample( 1:ncol(G), L)# Actual position of the causal effect

Y <- matrix(0 , ncol=  2^7 , nrow = N)
for ( i in 1:N){
  for ( l in 1:L){
    Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]


  }

}


#### Running susiF
m1 <- susiF   (Y=Y, X=G,L=20,L_start=11   )
m1$cs
