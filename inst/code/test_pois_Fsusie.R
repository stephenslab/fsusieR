library(susiF.alpha)
library(sim1000G)
library(mvPoisVA)
lev_res=7
N=20
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 200000 , min_maf = 0.02 , max_maf = NA )

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4)

startSimulation( vcf )



'%!in%' <- function(x,y)!('%in%'(x,y))
id = c()
for(i in 1:100) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:N,] + SIM$gt2[1:N,]


caca <- genotypes
print(dim(genotypes))

str(genotypes)
load("~/inst/check_pois_perf_weak.RData")
library(gplots)
res <-list()
for (o  in (length(res)+1):2000) {

  L <- sample(1:20, size=1)#actual number of effect
  print(L)
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- abs(simu_IBSS_per_level(lev_res=lev_res)$sim_func) #functional effect for effect l
  }


  tt <- sample(0:4,1)
  G <- genotypes

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(0,ncol=2^lev_res, nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
    Y[i, ]<- rpois (n=length(Y[i,]),lambda=Y[i,])

  }

  m1 <- susiF(Y=Y, X=G,L=20,L_start=5 ,nullweight=10 , maxit=10)
  m2 <- HF_susiF(Y=Y, X=G,L=20,L_start=5 ,nullweight=10 , maxit=10)
  m1$cs
  m1$est_pi


  cal_purity <- function(l_cs,X){
    tt <- list()
    for (k in 1:length(l_cs)){
      if(length(unlist(l_cs[[k]]))==1 ){
        tt[[k]] <- 1
      }else{
        x <-abs( cor(X[,unlist(l_cs[[k]]   ) ]))


        tt[[k]] <-  min( x[col(x) != row(x)])
      }
    }
    return( tt )
  }

  out <- c( length(m1$cs), #number of CS
            length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
            Reduce("+",sapply(1:length(m1$cs), function(k)
              ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m1$cs, X=as.matrix(G)),#mean purity
            length(m2$cs), #number of CS
            length(which(true_pos%in% do.call(c, m2$cs))), #number of effect found
            Reduce("+",sapply(1:length(m2$cs), function(k)
              ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m2$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m2$cs, length)),  #CS size
            L,tt)
  res[[o]] <- out
  print(res)
  save(res, file="~/inst/check_pois_perf_weak.RData")
}
