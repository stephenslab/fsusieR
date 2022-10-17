library(susiF.alpha)
library(sim1000G)
set.seed(1)
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 200 , min_maf = 0.02 , max_maf = NA )

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4 )

startSimulation( vcf )



id = c()
for(i in 1:100) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:100,] + SIM$gt2[1:100,]

print(dim(genotypes))

str(genotypes)

library(gplots)

  res <- list()
  for (o  in (length(res)+1):100) {
       L <- sample(4:10, size=1)
       print(L)
       lf <-  list()
       for(l in 1:L){
         lf[[l]] <- simu_IBSS_per_level(lev_res=5)$sim_func
       }


       tt <- sample(1:7,1)
       G <- genotypes#(tt*101):((tt+3)*102)]

       if( length(which(apply(G,2,var)==0))>0){
         G <- G[,-which(apply(G,2,var)==0)]
       }
       # G <- matrix( rnorm(100*300), nrow = 100)
       true_pos <- sample( 1:ncol(G), L)

       Y <- matrix(rnorm((2^5)*100 ,sd=1), nrow = 100)
       for ( i in 1:100){
         for ( l in 1:L){
           Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
         }
       }

       m1 <- susiF(Y=Y, X=G,L=20 ,L_start=4 ,nullweight=10,  prior="mixture_normal", cal_obj =FALSE, tol = 1e-6, maxit=10)
       m1$cs
       m1$est_pi
       true_pos[order(true_pos)]
       m2 <- susiF(Y=Y, X=G,L=20,L_start=4 ,nullweight=10 , maxit=10)
       m2$cs
       out <- c( length(m1$cs), length(which(true_pos%in% do.call(c, m1$cs))),
                 length(m2$cs),
                 length(which(true_pos%in% do.call(c, m2$cs))),
                 L)
       res[[o]] <- out
       print(res)
       save(res, file="check_L_accuracy.RData")
  }
