
# Test oneway -------------------------------------------------------------


media<-"M1"
contaminazione<-'C0'
ll_cont<-1
ll_m<-1
ll_sd<-4
n_i=20

alpha=0.05
N=500
B=1000
k=3
n=n_i*k
length_grid<-25
grid=seq(0,1,length.out = length_grid)
eff<-0.95
tol=1e-20
p_cont_vec<-c(0.05,0.1)
M_vec<-c(1,5,10)
sd_vec<-c(0.2,1,1.8,2.6,3.4,4.2,5)/length_grid



alpha_vec<- c(0,.05,.10,.25,.50)#c(0,.025,.05,.075,.10,.15,.25,.35,.50)
beta_vec<-c(0,.05,.10,.25,.50)#c(0,.025,.05,.075,.10,.15,.25,.35,.50)

p_cont<-p_cont_vec[ll_cont]
M_given<-M_vec[ll_m]
sd_given<-sd_vec[ll_sd]
data_out<-simulate_data_oneway(k=k,mean = media,con=contaminazione,p = p_cont,n_i = n_i,M = M_given,sd = sd_given,grid =grid)
data=data_out$data
label=data_out$label
grid<-data_out$grid
X_fdata<-fdata(t(data),argvals = grid)
# RoFanova ----------------------------------------------------------------
mu0=func.trim.FM(X_fdata,trim=0.2)
per_list_median<-RoFanova_oneway_perm(X_fdata,label,B = B,eff=eff,family="median",mu0_g=mu0)
pvalue_median<-per_list_median$pval
per_list_huber<-RoFanova_oneway_perm(X_fdata,label,B = B,eff=eff,family="huber",mu0_g=mu0)
pvalue_huber<-per_list_huber$pval
per_list_bisquare<-RoFanova_oneway_perm(X_fdata,label,B = B,eff=eff,family="bisquare",mu0_g=mu0)
pvalue_bisquare<-per_list_bisquare$pval
per_list_hampel<-RoFanova_oneway_perm(X_fdata,label,B = B,eff=eff,family="hampel",mu0_g=mu0)
pvalue_hampel<-per_list_hampel$pval
per_list_optimal<-RoFanova_oneway_perm(X_fdata,label,B = B,eff=eff,family="optimal",mu0_g=mu0)
pvalue_optimal<-per_list_optimal$pval



# Test two way ------------------------------------------------------------

library(mvtnorm)

library(fda.usc)
library(fda)
library(pbmcapply)
library(fdANOVA)
library(Rcpp)
library(parallel)
library(robustbase)
source("../fun_RoFanova.R")
sourceCpp('../fun.cpp')



media<-"M1"
contaminazione<-'C0'
n_i=20
N=500
B=1000
alpha_test=0.05
ll_cont<-1
ll_m<-1
ll_sd<-2
ll_alpha<-1
ll_beta<-1
k_1=k_2=2
n=n_i*k_1*k_2
length_grid<-25
grid=seq(0,1,length.out = length_grid)
eff<-0.95
tol=1e-20
p_cont_vec<-c(0.05,0.1)
M_vec<-c(1,5,25)
sd_vec<-c(.15,.30,.45)#seq(0.5,5.3,by=0.8)/10
alpha_vec<- c(0,.05,.10,.25,.50)#c(0,.025,.05,.075,.10,.15,.25,.35,.50)
beta_vec<-c(0,.05,.10,.25,.50)#c(0,.025,.05,.075,.10,.15,.25,.35,.50)

alpha<-alpha_vec[ll_alpha]
beta<-beta_vec[ll_beta]
p_cont<-p_cont_vec[ll_cont]
M_given<-M_vec[ll_m]
sd_given<-sd_vec[ll_sd]

data_out<-simulate_data_twoway(con=contaminazione,k_1=k_1,k_2=k_2,alpha = alpha,beta=beta,
                               p = p_cont,n_i = n_i,M = M_given,sd = sd_given,grid =grid)
data=data_out$data
label_1=data_out$label_1
label_2=data_out$label_2
grid<-data_out$grid
X_fdata<-fdata(t(data),argvals = grid)

# RoFanova ----------------------------------------------------------------
mu0=func.trim.FM(X_fdata,trim=0.2)
scale_glob<-rofanova:::scale_fun_pw(X_fdata,eff = eff,tol=tol, maxit = 50,mu0_g=mu0)
per_list_median<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="median",mu0_g=mu0,scale = scale_glob)
pvalue_median_vec<-per_list_median$pval_vec
per_list_huber<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="huber",mu0_g=mu0,scale=scale_glob)
pvalue_huber_vec<-per_list_huber$pval_vec
per_list_bisquare<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="bisquare",mu0_g=mu0,scale=scale_glob)
pvalue_bisquare_vec<-per_list_bisquare$pval_vec
per_list_hampel<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="hampel",mu0_g=mu0,scale=scale_glob)
pvalue_hampel_vec<-per_list_hampel$pval_vec
per_list_optimal<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="optimal",mu0_g=mu0,scale=scale_glob)
pvalue_optimal_vec<-per_list_optimal$pval_vec



