#  Utils -------------------------------------------------------------------
sum_fdata<-function(x){
  out<-x
  out$data<-matrix(colSums(x$data),1)
  return(out)
}
diff_fdata<-function(x1,x2){
  out<-x1
  out$data<-x1$data-x2$data
  return(out)
}


diff_fdata_sur<-function(x1,x2){
  n<-dim(x1$data)[1]
  out<-x1
  for(ii in 1:n)out$data[ii,,]=x1$data[ii,,]-x2$data[ii,,]
  return(out)
}
sum_fdata_sur<-function(x1,x2){
  n<-dim(x1$data)[1]
  out<-x1
  for(ii in 1:n)out$data[ii,,]=x1$data[ii,,]+x2$data[ii,,]
  return(out)
}

func.mean_sur<-function (fdataobj)
{
  n<-dim(fdataobj$data)[1]
  pro<-array(NA,dim = c(1,dim(fdataobj$data)[2],dim(fdataobj$data)[3]))
  pro[1,,] <-apply(fdataobj$data , c(2,3), base::mean)
  fdataobj$data<-pro
  fdataobj$names$main <- "mean"
  fdataobj
}
func.var_sur<-function (fdataobj)
{

  n<-dim(fdataobj$data)[1]
  pro<-array(NA,dim = c(1,dim(fdataobj$data)[2],dim(fdataobj$data)[3]))
  pro[1,,] <-apply(fdataobj$data , c(2,3), stats::var)
  fdataobj$data<-pro
  fdataobj$names$main <- "var"
  fdataobj
}
frac_fdata_sur<-function(x1,x2){
  n<-dim(x1$data)[1]
  out<-x1
  for(ii in 1:n)out$data[ii,,]=x1$data[ii,,]/(x2$data[ii,,]+1e-15)
  return(out)
}
standardize_sur<-function(x,mu0,sig0=NA){
  if(is.na(sig0))sig0<-fdata(array(1,dim = c(1,dim(x$data)[2],dim(x$data)[3])),argvals = x$argvals)
  n<-dim(x$data)[1]
  mu0_rep<-x
  sig0_rep<-x

  for(ii in 1:n)mu0_rep$data[ii,,]=mu0$data[1,,]
  for(ii in 1:n)sig0_rep$data[ii,,]=sig0$data[1,,]
  frac_fdata_sur(diff_fdata_sur(x,mu0_rep),sig0_rep)

}
ex_fdata<-function(x,index){
  if(is.logical(index))index=which(index==TRUE)
  out<-x
  out$data<-x$data[index,,,drop=FALSE]
  return(out)

}
# acc_fdata_sur<-function(x){
#   out<-x
#   pro<-array(NA,dim = c(1,dim(x$data)[2],dim(x$data)[3]))
#   pro[1,,] <-apply(x$data , c(2,3), sum)
#   out$data<-pro
#   return(out)
# }

# int_sur<-function(x){
#
#   grid_s<-x$argvals[[1]]
#   grid_t<-x$argvals[[2]]
#   delta_s=1/(length(grid_s)-1)
#   delta_t=1/(length(grid_t)-1)
#   delta_s*delta_t*sum(x$data[1,,])
#
# }
# tra_data<-function(X_fdata){
#
#
#   out<-X_fdata
#
#   sss<-powerTransform(norm_fdata_c_sur(out))
#   out$data= bcPower( out$data+abs(min(out$data))+0.00001,0.5)
#   return(out)
#
# }
# alignment_fdata_sur<-function(X_fdata,ref){
#
#   if(is.numeric(ref)){
#     X_1<-ex_fdata(X_fdata,ref)
#     n=dim(X_fdata$data)[1]
#     n_minus<-(1:n)[-ref]
#   }
#   else{
#     X_1<-ref
#     n=dim(X_fdata$data)[1]
#     n_minus<-(1:n)
#   }
#   X_fdata_new<-X_fdata
#
#   for (kkk in 1:length(n_minus)) {
#
#     X_2<-ex_fdata(X_fdata,n_minus[kkk])
#     grid<-expand.grid(-4:4,-4:4)
#     norm_vec<-numeric()
#     X_2_new<-list()
#     for (ii in 1:length(grid[,1])) {
#       data_new<-X_2$data
#       dim_vec<-dim(data_new)
#       delta_ind<-as.numeric(grid[ii,])
#       if(delta_ind[1]>0){
#         data_new[1,(1+delta_ind[1]):dim_vec[2],]<-X_2$data[1,1:(dim_vec[2]-delta_ind[1]),]
#       }
#       else if (delta_ind[1]<0){
#         delta_ind[1]<-abs(delta_ind[1])
#         data_new[1,1:(dim_vec[2]-delta_ind[1]),]<-X_2$data[1,(1+delta_ind[1]):dim_vec[2],]
#       }
#
#       if(delta_ind[2]>0){
#         data_new[1,,(1+delta_ind[2]):dim_vec[3]]<-data_new[1,,1:(dim_vec[3]-delta_ind[2])]
#       }
#       else if (delta_ind[2]<0){
#         delta_ind[2]<-abs(delta_ind[2])
#         data_new[1,,1:(dim_vec[3]-delta_ind[2])]<-data_new[1,,(1+delta_ind[2]):dim_vec[3]]
#
#       }
#       X_2_new[[ii]]=X_2
#       X_2_new[[ii]]$data=data_new
#       norm_vec[ii]<-norm_fdata_c_sur(diff_fdata_sur(X_2_new[[ii]],X_1))
#     }
#
#     min_ind<-which(norm_vec==min(norm_vec))
#     X_fdata_new$data[n_minus[kkk],,]<-X_2_new[[min_ind]]$data
#   }
#
#   return(X_fdata_new)
#
# }
# alignment_fdata_sur_cv<-function(X_fdata){
#   n=dim(X_fdata$data)[1]
#   median_diff<-numeric()
#   par_fun<-function(ii) {
#     X_fdata_new<-alignment_fdata_sur(X_fdata,ii)
#     mean(parwise_diff_sur(X_fdata_new ))
#   }
#   mean_diff<-mclapply(1:n,par_fun,mc.cores=detectCores()-1)
#
#   ind_min<-max(which(unlist(mean_diff)==min(unlist(mean_diff))))
#   X_fdata_new<-alignment_fdata_sur(X_fdata,ind_min)
#   return(X_fdata_new)
# }
#
# der_sur<-function(X_fdata,nderiv=1,dim="s"){
#   n=dim(X_fdata$data)[1]
#   ss = X_fdata$argvals[[1]]
#   tt = X_fdata$argvals[[2]]
#
#   der_fdata<-X_fdata
#
#
#   for (ii in 1:n) {
#     cat(ii)
#     data<-X_fdata$data[ii,,]
#     ns <- nrow(data)
#     nt <- ncol(data)
#     if (dim == "s") {
#       res = matrix(NA, nrow = ns, ncol = nt)
#       for (i in 1:nt) {
#         a = diff(data[,i], differences = nderiv)/(ss[2:ns] -
#                                                     ss[1:(ns - 1)])^nderiv
#         ab = matrix(NA, ncol = ns, nrow = 2)
#         ab[1, 2:ns] = a
#         ab[2, 1:(ns - 1)] = a
#         res[,i ] = colMeans(ab, na.rm = TRUE)
#       }
#       der_fdata$data[ii,,]<-res
#     }
#     if (dim == "t") {
#       res = matrix(NA, nrow = ns, ncol = nt)
#       for (i in 1:ns) {
#         a = diff(data[i,], differences = nderiv)/(tt[2:nt] -
#                                                     tt[1:(nt - 1)])^nderiv
#         ab = matrix(NA, ncol = nt, nrow = 2)
#         ab[1, 2:nt] = a
#         ab[2, 1:(nt - 1)] = a
#         res[i, ] = colMeans(ab, na.rm = TRUE)
#       }
#       der_fdata$data[ii,,]<-res
#     }
#   }
#   return(der_fdata)
# }
#
# twosample_test_sur<-function(X_1,X_2,B=100,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,...){
#
#
#   n_1<-dim(X_1$data)[1]
#   n_2<-dim(X_2$data)[1]
#   n=n_1+n_2
#   X=fdata(abind(X_1$data,X_2$data,along = 1),argvals = X_1$argvals)
#   scale_1<-scale_fun_pw_sur(X_1,eff = eff,tol=tol, maxit = maxit)#4)
#   scale_2<-scale_fun_pw_sur(X_2,eff = eff,tol=tol, maxit = maxit )#4)
#   # scale_1<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
#   # scale_2<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
#   if(n_1==2|n_2==2){
#     group_mean_1<-FlocScaleM_sur(X_1,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
#     group_mean_2<-FlocScaleM_sur(X_2,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
#   }
#   else{
#     group_mean_1<-FlocScaleM_sur(X_1,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_1)$mu
#     group_mean_2<-FlocScaleM_sur(X_2,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_2)$mu
#   }
#   T_obs<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))
#   # # if(is.na(scale_glob))scale_glob<-scale_fun_pw_sur(X,eff = eff,tol=tol, maxit = maxit,... )#5
#   if(n_1==2|n_2==2){
#     per_fun<-function(kkk){
#       perm_comb<-sample(1:n,replace = F)
#       X_1per<-ex_fdata(X,perm_comb[1:n_1])
#       X_2per<-ex_fdata(X,perm_comb[(n_1+1):n])
#       group_mean_1<-FlocScaleM_sur(X_1per,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
#       group_mean_2<-FlocScaleM_sur(X_2per,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
#       T_i<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))
#
#       return(T_i)
#     }
#   }
#   else{
#     per_fun<-function(kkk){
#       perm_comb<-sample(1:n,replace = F)
#       X_1per<-ex_fdata(X,perm_comb[1:n_1])
#       X_2per<-ex_fdata(X,perm_comb[(n_1+1):n])
#       scale_1<-scale_fun_pw_sur(X_1per,eff = eff,tol=tol, maxit = maxit)#4)
#       scale_2<-scale_fun_pw_sur(X_2per,eff = eff,tol=tol, maxit = maxit )#4)
#       group_mean_1<-FlocScaleM_sur(X_1per,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_1)$mu
#       group_mean_2<-FlocScaleM_sur(X_2per,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_2)$mu
#       T_i<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))
#
#       return(T_i)
#     }
#   }
#   T_per<-unlist(mclapply(1:B,per_fun,mc.cores = detectCores()-1))
#
#   pvalue<-sum(T_per>=(T_obs+1e-10))/B
#   out<-list(pvalue=pvalue,
#             T_obs=T_obs,
#             T_per=T_per)
#   return(out)
# }
# pairwise_comparisons<-function(X_fdata,label_1,label_2,...) {
#
#   k_1=length(unique(label_1))
#   k_2=length(unique(label_2))
#
#   pairwise<-cbind(expand.grid(1:k_2,1:k_1)[2],expand.grid(1:k_2,1:k_1)[1])
#   mat<-matrix(0,dim(pairwise)[1],dim(pairwise)[1])
#   for (ii in 1:dim(pairwise)[1]) {
#     for (jj in 1:ii) {
#       cat(paste(as.numeric(pairwise[ii,])));cat(" - ");cat(as.numeric(pairwise[jj,]));cat("\n")
#       iii=pairwise[ii,1]; jjj=pairwise[ii,2]
#       kkk=pairwise[jj,1]; lll=pairwise[jj,2]
#       mat[ii,jj]<-twosample_test_sur(X_1=ex_fdata(X_fdata,label_1==iii&label_2==jjj),X_2=ex_fdata(X_fdata,label_1==kkk&label_2==lll),...)$pvalue
#     }
#   }
#   colnames(mat)<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))
#   rownames(mat)<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))
#
#
#   return(mat)
#
#
# }
# get_matrix_pairwise<-function(mod,label_1,label_2,fac=1){
#
#   k_1=length(unique(label_1))
#   k_2=length(unique(label_2))
#   int<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))
#   matrix_list<-list()
#   ind<-ifelse(fac==1,k_1,k_2)
#   for (ii in 1:ind) {
#     if(fac==1)     ind_ii<-which(substr(int, 1, 1)==ii)
#     else  ind_ii<-which(substr(int, 3, 3)==ii)
#     matrix_list[[ii]]<-mod[ind_ii,ind_ii]
#
#   }
#   return(matrix_list)
# }
