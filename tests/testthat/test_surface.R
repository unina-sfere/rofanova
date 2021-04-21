

# Test one-way surface ----------------------------------------------------



media<-"M2"
contaminazione<-'C0'
ll_cont<-1
ll_m<-1
ll_sd<-4
n_i=20

alpha=0.05
N=500
B=20
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
data_out<-simulate_data(scenario="one-way surface",k_1=k,mean = media,con=contaminazione,p = p_cont,n_i = n_i,M = M_given,sd = sd_given,grid =grid)
label_1=data_out$label
X_fdata<-data_out$X_fdata
slice=20
plot(X_fdata$data[1,slice,],type="l",ylim=c(0,0.2))
sapply(1:20,function(ii)lines((X_fdata$data[ii,slice,]),type="l"))
sapply(21:40,function(ii)lines((X_fdata$data[ii,slice,]),type="l",col=2))
sapply(41:60,function(ii)lines((X_fdata$data[ii,slice,]),type="l",col=3))

cores=1
# RoFanova one-way
per_list_median<-rofanova(X_fdata,label_1,B = B,eff=eff,family="median",maxit=maxit,cores=cores)
pvalue_median_vec<-per_list_median$pval_vec
per_list_huber<-rofanova(X_fdata,label_1,B = B,eff=eff,family="huber",maxit=maxit,cores=cores)
pvalue_huber_vec<-per_list_huber$pval_vec
per_list_bisquare<-rofanova(X_fdata,label_1,B = B,eff=eff,family="bisquare",maxit=maxit,cores=cores)
pvalue_bisquare_vec<-per_list_bisquare$pval_vec
per_list_hampel<-rofanova(X_fdata,label_1,B = B,eff=eff,family="hampel",maxit=maxit,cores=cores)
pvalue_hampel_vec<-per_list_hampel$pval_vec
per_list_optimal<-rofanova(X_fdata,label_1,B = B,eff=eff,family="optimal",maxit=maxit,cores=cores)
pvalue_optimal_vec<-per_list_optimal$pval_vec



# Test two-way surface ----------------------------------------------------


media<-"M2"
contaminazione<-'C0'
n_i=20
N=500
B=100
alpha_test=0.05
ll_cont<-1
ll_m<-1
ll_sd<-2
ll_alpha<-2
ll_beta<-2
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
data_out<-simulate_data(scenario="two-way surface",con=contaminazione,k_1=k_1,k_2=k_2,alpha = alpha,beta=beta,
                        p = p_cont,n_i = n_i,M = M_given,sd = sd_given,grid =grid)
label_1=data_out$label_1
label_2=data_out$label_2
X_fdata<-data_out$X_fdata
plot(X_fdata$data[1,slice,],type="l",ylim=c(0,1))
sapply(1:20,function(ii)lines((X_fdata$data[ii,slice,]),type="l"))
sapply(21:40,function(ii)lines((X_fdata$data[ii,slice,]),type="l",col=2))
sapply(41:60,function(ii)lines((X_fdata$data[ii,slice,]),type="l",col=3))


cores=1
# RoFanova two-way
per_list_median<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="median",maxit=maxit,cores=cores)
pvalue_median_vec<-per_list_median$pval_vec
per_list_huber<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="huber",maxit=maxit,cores=cores)
pvalue_huber_vec<-per_list_huber$pval_vec
per_list_bisquare<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="bisquare",maxit=maxit,cores=cores)
pvalue_bisquare_vec<-per_list_bisquare$pval_vec
per_list_hampel<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="hampel",maxit=maxit,cores=cores)
pvalue_hampel_vec<-per_list_hampel$pval_vec
per_list_optimal<-rofanova(X_fdata,label_1,label_2,B = B,eff=eff,family="optimal",maxit=maxit,cores=cores)
pvalue_optimal_vec<-per_list_optimal$pval_vec



