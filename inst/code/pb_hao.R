# install.packages("remotes")
# remotes::install_github("stephenslab/susiF.alpha@v0.1.210")
# v1.0.12

rm(list=ls())
library(susiF.alpha)
sessionInfo()
a =  readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/pb_hao/fsusie.yuqi_mqtl.tad808.uni_Fsusie.mixture_normal_per_scale (3).rds")
X = a$input_data$X[1:100,]
Y <- a$input_data[[2]][1:100,]
result1=  susiF(Y= Y  ,X = X ,pos =a$input_data$pos,   L=20 ,prior = "mixture_normal_per_scale",max_SNP_EM = 100)
save(result1, file='D:/Document/Serieux/Travail/Data_analysis_and_papers/pb_hao/res808.RData' )




rm(list=ls())
library(susiF.alpha)
a =  readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/pb_hao/fsusie.yuqi_mqtl.tad1182.uni_Fsusie.mixture_normal_per_scale (3).rds")
X = a$input_data$X
Y =a$input_data[[2]]
result2=  susiF(X = a$input_data$X , a$input_data[[2]]  ,pos =a$input_data$pos,   L=20 ,prior = "mixture_normal_per_scale",max_SNP_EM = 1000)
save(result2, file='D:/Document/Serieux/Travail/Data_analysis_and_papers/pb_hao/res1182.RData' )
plot_susiF(result2
)

