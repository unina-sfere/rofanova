
# setwd("~/CloudStation/ROFANOVA/Case-study")
library(mvtnorm)

library(fda.usc)
library(fda)
library(pbmcapply)
library(fdANOVA)
library(Rcpp)
library(parallel)
library(robustbase)
library(abind)
library(car)
source("fun_RoFanova.R")
source("plot_functions.R")
sourceCpp('fun.cpp')

load(file = "mat_list0.RData")
mat_list<-mat_list
lab_lev<-label_lev
label_lay<-label_lay
eff=0.95
B=1000
# Get fdata ---------------------------------------------------------------
length_grid_s<-dim(mat_list[[1]])[1]
length_grid_t<-dim(mat_list[[1]])[2]
s_grid<-seq(0,1,length.out = length_grid_s)
t_grid<-seq(0,1,length.out = length_grid_t)
zero_coor<-c(0.7710,0.8347)
n<-length(lab_lev)
z<-array(NA,dim=c(n,length_grid_s,length_grid_t))
for (ii in 1:n) z[ii,,]<-mat_list[[ii]]
X_fdata<-fdata(z,argvals = list(s_grid,t_grid))
label_1<-label_lev
label_2<-label_lay
# ALIGN -------------------------------------------------------------------
#


load(file="X_fdatas0.RData")
#
# mod_bi<-RoFanova_twoway_sur(X_fdata_new,label_1,label_2,eff=eff,family="bisquare",maxit=200)
#
# k_1=k_2=6
# kkk=1
# sum_full_list<-list()
# for (ii in 1:k_1) {
#   for (jj in 1:k_2) {
#     sum_full_list[[kkk]]<-standardize_sur(ex_fdata(X_fdata_new,label_1==ii&label_2==jj),mod_bi$group_mean_ij[[ii]][[jj]])
#     kkk=kkk+1
#   }
# }
# sum_full_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
# res_fdata<-fdata(sum_full_data,argvals =X_fdata_new$argvals )
# norm_res<-norm_fdata1_c_sur(sum_full_i)
# save(res_fdata,norm_res, mod_bi, file = "./res_analy.RData")

# ALL data ----------------------------------------------------------------

X_fdata_i<-X_fdata_new
B=10
cores=5
maxit=5
# RoFanova two-way
per_list_median<-rofanova(X_fdata_i,label_1,label_2,B = B,eff=eff,family="median",maxit=maxit,cores=cores)
pvalue_median_vec<-per_list_median$pval_vec
per_list_huber<-rofanova(X_fdata_i,label_1,label_2,B = B,eff=eff,family="huber",maxit=maxit,cores=cores)
pvalue_huber_vec<-per_list_huber$pval_vec
per_list_bisquare<-rofanova(X_fdata_i,label_1,label_2,B = B,eff=eff,family="bisquare",maxit=maxit,cores=cores)
pvalue_bisquare_vec<-per_list_bisquare$pval_vec
per_list_hampel<-rofanova(X_fdata_i,label_1,label_2,B = B,eff=eff,family="hampel",maxit=maxit,cores=cores)
pvalue_hampel_vec<-per_list_hampel$pval_vec
per_list_optimal<-rofanova(X_fdata_i,label_1,label_2,B = B,eff=eff,family="optimal",maxit=maxit,cores=cores)
pvalue_optimal_vec<-per_list_optimal$pval_vec


# RoFanova one-way
per_list_median<-rofanova(X_fdata_i,label_1,B = B,eff=eff,family="median",maxit=maxit,cores=cores)
pvalue_median_vec<-per_list_median$pval_vec
per_list_huber<-rofanova(X_fdata_i,label_1,B = B,eff=eff,family="huber",maxit=maxit,cores=cores)
pvalue_huber_vec<-per_list_huber$pval_vec
per_list_bisquare<-rofanova(X_fdata_i,label_1,B = B,eff=eff,family="bisquare",maxit=maxit,cores=cores)
pvalue_bisquare_vec<-per_list_bisquare$pval_vec
per_list_hampel<-rofanova(X_fdata_i,label_1,B = B,eff=eff,family="hampel",maxit=maxit,cores=cores)
pvalue_hampel_vec<-per_list_hampel$pval_vec
per_list_optimal<-rofanova(X_fdata_i,label_1,B = B,eff=eff,family="optimal",maxit=maxit,cores=cores)
pvalue_optimal_vec<-per_list_optimal$pval_vec



# FANOVA
per_list_fanova<-fanova_twoway_perm_sur(X_fdata_i,label_1,label_2,B=B)
pvalue_numden_vec<-per_list_fanova$pval_mat[,1]
pvalue_full_vec<-per_list_fanova$pval_mat[,2]

p_value<-cbind(pvalue_numden_vec,pvalue_full_vec,pvalue_median_vec,pvalue_huber_vec,pvalue_bisquare_vec,pvalue_hampel_vec,pvalue_optimal_vec)
colnames(p_value)<-c("NUMDEN","FULL","MED","HUB","BIS","HAM","OPT")

# Interaction_plot

# mod_bisquare_i<-RoFanova_twoway_sur(X_fdata_i,label_1,label_2,eff=eff,family="bisquare",maxit=200)
# x11()
# interaction_plot_sur(mod_bisquare_i)
# x11()
# mean_plot_sur(mod_bisquare_i,ylim=c(0,6))

#post hoc test
# mod<-pairwise_comparisons(X_fdata_i ,label_1,label_2,B=B,eff=eff,family="bisquare",maxit=200)
# mat_1<-get_matrix_pairwise(mod,label_1,label_2)
# mat_2<-get_matrix_pairwise(mod,label_1,label_2,fac=2)

save(per_list_median,per_list_huber,per_list_bisquare,per_list_hampel,per_list_optimal,per_list_fanova,p_value,file = "Case-study_new0.RData")
