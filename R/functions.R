
# Simulate data -----------------------------------------------------------


#' @export
simulate_data<-function(scenario="one-way",mean="M1",con="C1",p=0.1,M=1,n_i=10,k_1=3,k_2=3,alpha=0.05,beta=0.05,sd=0.01,grid=seq(0,1,length.out = 30),c=1,err="s"){

  if(scenario=="one-way")
    return(simulate_data_oneway(mean=mean,con=con,p=p,M=M,n_i=n_i,k=k_1,sd=sd,grid=grid,c=c))
  if(scenario=="two-way")
    return(simulate_data_twoway(con=con,p=p,M=M,n_i=n_i,k_1=k_1,k_2=k_2,alpha=alpha,beta=beta,sd=sd,grid=grid,c=c))
  if(scenario=="one-way surface")
    return(simulate_data_oneway_sur(mean=mean,con=con,p=p,M=M,n_i=n_i,k=k_1,sd=sd,grid=grid,err=err))
  if(scenario=="two-way surface")
    return(simulate_data_twoway_sur(con=con,p=p,M=M,n_i=n_i,k_1=k_1,k_2=k_2,alpha=alpha,beta=beta,sd=sd,grid=grid,err=err))


}
simulate_data_oneway<-function(mean="M1",con="C1",p=0.1,M=1,n_i=10,k=3,sd=0.01,grid=seq(0,1,length.out = 30),c=1){

print("Simulated data one-way")
  if(mean=="M1"){
    mean_function<-function(t,i)t*(1-t)
  }
  else if(mean=="M2"){
    mean_function<-function(t,i)(t^(i))*(1-t)^(6-i)
  }
  else if(mean=="M3"){
    mean_function<-function(t,i)t^(i/5)*(1-t)^(6-(i/5))
  }


  if(con=="C0"){
    cont_function<-function(n_g,p_g,M_g,ii)0
  }
  else if(con=="C1"){
    cont_function<-function(n_g,p_g,M_g,ii)M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
  }
  else if(con=="C2"){
    cont_function<-function(n_g,p_g,M_g,ii){
      T_i<-runif(n_g,0,0.75)
      matrix_old<-M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
      matrix_new<-matrix(0,length(matrix_old[,1]),n_g)
      for (ii in 1:n_g) {
        ind<-(round(T_i[ii]*(length(matrix_old[,1])-1))+1)
        matrix_new[ind:length(matrix_old[,1]),ii]<-matrix_old[ind:length(matrix_old[,1]),ii]

      }
      return(matrix_new)
    }
  }
  else if(con=="C3"){
    cont_function<-function(n_g,p_g,M_g,ii){
      g<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      h<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=(sd+2)^2,theta=1))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*h
      return(out)
    }
  }
  else if(con=="C4"){
    cont_function<-function(n_g,p_g,M_g,ii)(-1)^(ii)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
  }
  else if(con=="C5"){
    cont_function<-function(n_g,p_g,M_g,ii){
      T_i<-runif(n_g,0,0.75)
      matrix_old<-(-1)^(ii)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      matrix_new<-matrix(0,length(matrix_old[,1]),n_g)
      for (ll in 1:n_g) {
        ind<-(round(T_i[ll]*(length(matrix_old[,1])-1))+1)
        matrix_new[ind:length(matrix_old[,1]),ll]<-matrix_old[ind:length(matrix_old[,1]),ll]

      }
      return(matrix_new)
    }
  }
  else if(con=="C6"){
    cont_function<-function(n_g,p_g,M_g,ii){
      g<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      h<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=(sd+2+(-1)^(ii))^2,theta=1))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*h
      return(out)
    }
  }
  if(con=="C3"|con=="C6"){
    data_list<-list()
    for (ll in 1:k) {

      data_list[[ll]]<-cont_function(n_i,p,M,ll)
    }

  }
  else{
    data_list<-list()
    for (ll in 1:k) {
      data_list[[ll]]<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)+cont_function(n_i,p,M,ll)
    }
  }

  data<-do.call("cbind", data_list)
  X_fdata<-fdata(t(data),argvals = grid)
  label<-rep(1:k,each=n_i)
  out=list(X_fdata=X_fdata,
           data=data,
           label=label,
           grid=grid)
  return(out)
}
simulate_data_twoway<-function(con="C1",n_i=10,k_1=3,k_2=3,p=0.1,M=1,alpha=0.05,beta=0.05,sd=0.01,grid=seq(0,1,length.out = 30),c=1){


  print("Simulated data Two-way")
  mean_function<-function(t,i,j)t*(1-t)
  f1_function<-function(t,i,j,alpha,beta)alpha*(-1)^(i)*abs(sin(4*pi*t))
  f2_function<-function(t,i,j,alpha,beta)beta*(-1)^(j)*ifelse(t>0.5,1,0)
  int_function<-function(t,i,j,alpha,beta)-f1_function(t,i,j,alpha,beta)*f2_function(t,i,j,alpha,beta)*ifelse(alpha>=0.25,1,0)


  Y<-function(t,i,j,alpha,beta)mean_function(t,i,j)+f1_function(t,i,j,alpha,beta)+f2_function(t,i,j,alpha,beta)+int_function(t,i,j,alpha,beta)



  if(con=="C0"){
    cont_function<-function(n_g,p_g,M_g,ii,jj)0
  }
  else if(con=="C1"){
    cont_function<-function(n_g,p_g,M_g,ii,jj)M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
  }
  else if(con=="C2"){
    cont_function<-function(n_g,p_g,M_g,ii,jj){
      T_i<-runif(n_g,0,0.75)
      matrix_old<-M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
      matrix_new<-matrix(0,length(matrix_old[,1]),n_g)
      for (ii in 1:n_g) {
        ind<-(round(T_i[ii]*(length(matrix_old[,1])-1))+1)
        matrix_new[ind:length(matrix_old[,1]),ii]<-matrix_old[ind:length(matrix_old[,1]),ii]

      }
      return(matrix_new)
    }
  }
  else if(con=="C3"){
    cont_function<-function(n_g,p_g,M_g,ii,jj){
      g<-Y(matrix(grid,length(grid),n_i),ii,jj,alpha,beta)+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      h<-Y(matrix(grid,length(grid),n_i),ii,jj,alpha,beta)+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=(sd+2)^2,theta=1))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*h
      return(out)
    }
  }
  else if(con=="C4"){
    cont_function<-function(n_g,p_g,M_g,ii,jj)(-1)^(ii)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
  }
  else if(con=="C5"){
    cont_function<-function(n_g,p_g,M_g,ii,jj){
      T_i<-runif(n_g,0,0.75)
      matrix_old<-(-1)^(ii)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      matrix_new<-matrix(0,length(matrix_old[,1]),n_g)
      for (ll in 1:n_g) {
        ind<-(round(T_i[ll]*(length(matrix_old[,1])-1))+1)
        matrix_new[ind:length(matrix_old[,1]),ll]<-matrix_old[ind:length(matrix_old[,1]),ll]

      }
      return(matrix_new)
    }
  }
  else if(con=="C6"){
    cont_function<-function(n_g,p_g,M_g,ii,jj){
      g<-Y(matrix(grid,length(grid),n_i),ii,jj,alpha,beta)+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      h<-Y(matrix(grid,length(grid),n_i),ii,jj,alpha,beta)+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=(sd+2+(-1)^(ii))^2,theta=1))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*h
      return(out)
    }
  }


  if(con=="C3"|con=="C6"){
    data_list<-list()
    kk=1
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        data_list[[kk]]<-cont_function(n_i,p,M,ii,jj)
        kk=kk+1
      }
    }

  }
  else{
    data_list<-list()
    kk=1
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        data_list[[kk]]<-Y(matrix(grid,length(grid),n_i),ii,jj,alpha,beta)+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)+cont_function(n_i,p,M,ii,jj)
        kk=kk+1
      }


    }
  }


  data<-do.call("cbind", data_list)
  X_fdata<-fdata(t(data),argvals = grid)
  label_1<-rep(1:k_1,each=n_i*k_2)
  label_2<-rep(rep(1:k_2,each=n_i),k_1)
  out=list(X_fdata=X_fdata,
           data=data,
           label_1=label_1,
           label_2=label_2,
           grid=grid)
  return(out)
}
simulate_data_oneway_sur<-function(mean="M1",con="C1",p=0.1,M=1,n_i=10,k=3,sd=0.01,grid=seq(0,1,length.out = 30),err="s"){

  print("Simulated data one-way surface")
  if(mean=="M1"){
    mean_function<-function(s,t,i)t*(1-t)+s*(1-s)
  }
  else if(mean=="M2"){
    mean_function<-function(s,t,i)(t^(i))*(1-t)^(6-i)+(s^(i))*(1-s)^(6-i)
  }


  if(con=="C0"){
    cont_function<-function(n_g,p_g,M_g,ii)0
   }
  else if(con=="C1"){
    cont_function<-function(n_g,p_g,M_g,ii)M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
  }

  expand_grid<-expand.grid(grid,grid)
  data<-array(0,c(k*n_i,length(grid),length(grid)))
  kk=1
  for (ll in 1:k) {

    for (mm in 1:n_i) {
      error<-sapply(1:length(grid),function(ii)rproc2fdata(1,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=1))$data)
      mat_err<-if(err=="s")error else t(error)
      data[kk,,] <-matrix(mean_function(expand_grid[,1],expand_grid[,2],ll),length(grid),length(grid))+mat_err+matrix(cont_function(1,p,M,ll),length(grid),length(grid))
      kk=kk+1
    }
  }


  X_fdata<-fdata(data,argvals = list(grid,grid))

  label<-rep(1:k,each=n_i)
  out=list(X_fdata=X_fdata,
           data=data,
           label=label,
           grid=grid)
  return(out)
}
simulate_data_twoway_sur<-function(con="C1",n_i=10,k_1=3,k_2=3,p=0.1,M=1,alpha=0.05,beta=0.05,sd=0.01,grid=seq(0,1,length.out = 30),err="s"){


  print("Simulated data Two-way surface")
  mean_function<-function(s,t,i,j)t*(1-t)+s*(1-s)
  f1_function<-function(s,t,i,j,alpha,beta)alpha*(-1)^(i)*(abs(sin(4*pi*t))+abs(sin(4*pi*s)))
  f2_function<-function(s,t,i,j,alpha,beta)beta*(-1)^(j)*ifelse(t>0.5,1,0)*ifelse(s>0.5,1,0)
  int_function<-function(s,t,i,j,alpha,beta)-f1_function(s,t,i,j,alpha,beta)*f2_function(s,t,i,j,alpha,beta)*ifelse(alpha>=0.25,1,0)
  Y<-function(s,t,i,j,alpha,beta)mean_function(s,t,i,j)+f1_function(s,t,i,j,alpha,beta)+f2_function(s,t,i,j,alpha,beta)+int_function(s,t,i,j,alpha,beta)
  if(con=="C0"){
    cont_function<-function(n_g,p_g,M_g)0
  }
  else if(con=="C1"){
    cont_function<-function(n_g,p_g,M_g)M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
  }
  expand_grid<-expand.grid(grid,grid)
  data<-array(0,c(k_1*k_2*n_i,length(grid),length(grid)))
  kk=1
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
    for (mm in 1:n_i) {
      error<-sapply(1:length(grid),function(ii)rproc2fdata(1,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      mat_err<-if(err=="s")error else t(error)
      data[kk,,] <-matrix(Y(expand_grid[,1],expand_grid[,2],ii,jj,alpha,beta),length(grid),length(grid))+mat_err+matrix(cont_function(1,p,M),length(grid),length(grid))
      kk=kk+1
    }
    }
  }
  X_fdata<-fdata(data,argvals = list(grid,grid))
  label_1<-rep(1:k_1,each=n_i*k_2)
  label_2<-rep(rep(1:k_2,each=n_i),k_1)
  out=list(X_fdata=X_fdata,
           data=data,
           label_1=label_1,
           label_2=label_2,
           grid=grid)
  return(out)
}


# Robust location/scale estimation -----------------------------------------

#si
FlocScaleM<-function (x, psi = "bisquare", eff = 0.95, maxit = 50, tol = 1e-04,mu0_g=NA,sig0_g=NA){
  kpsi <- switch(psi, bisquare = 1, huber = 2, optimal = 3, median=5,mean=6,hampel=7,8)
  if (kpsi == 8) {
    stop(paste0(psi, " - No such rho function"))
  }
  else if (kpsi==6){
    mu0 = func.med.FM(x,trim=0.1)
    sig0 = sqrt(func.trimvar.FM(x,trim=0.2))
    mu=func.mean(x)
    resu <- list(mu = mu,mu0=func.med.FM(x,trim=0.1),sig0=sig0)
  }
  else if (kpsi==5){

    ktun=keff=1

    if(is.na(mu0_g)) mu0=func.mean(x)
    else mu0=mu0_g
    if(is.na(sig0_g)) sig0=sqrt(func.var(x))
    else sig0=sig0_g

    if (max(sig0) < 1e-10) {
      mu = 0
      sigma = 0
    }
    else {
      resi_start =norm_fdata_c((x-mu0)/sig0 )
      ww_start= wfun_c(resi_start, kpsi,ktun)
      if(all(ww_start==0))mu0=func.trim.FM(x,trim=0.5)

      mu=iteration(x,mu0,sig0,kpsi,ktun,tol, maxit)

    }

    resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

  }
  else {
    family=psi
    kBis = c(3.44, 3.88, 4.685)
    kHub = c(0.732, 0.981, 1.34)
    kHampel<-rbind(c(0.9539156, 1.9078312, 3.8156624),c(NA,NA,NA),c(1.3521, 3.1549, 7.2112))
    kOpt<-c(0.87,0.94,1.06)
    efis = c(0.85, 0.9, 0.95)
    keff = match(eff, efis)
    if(kpsi==7)ktun<-kHampel[keff,]
    else if(kpsi==1)ktun<-kBis[keff]
    else if(kpsi==2)ktun<-kHub[keff]
    else if(kpsi==3)ktun<-kOpt[keff]
    if(is.na(ktun[1])){
      print(c(eff, " No such eff"))
      keff = 0
    }
    else {

      if(is.na(mu0_g)) mu0=func.mean(x)
      else mu0=mu0_g
      if(is.na(sig0_g)) sig0=sqrt(func.var(x))
      else sig0=sig0_g

      if (max(sig0) < 1e-10) {
        mu = 0
        sigma = 0
      }
      else {
        resi_start =norm_fdata_c((x-mu0)/sig0 )
        ww_start= robustbase::Mwgt(resi_start, ktun, family)
        if(all(ww_start==0))mu0=func.trim.FM(x,trim=0.5)

        mu=iteration_ho(x,mu0,sig0,matrix(ktun),family,tol, maxit)

      }
      resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

    }
  }

  return(resu)
}
FlocScaleM_sur<-function (x, psi = "bisquare", eff = 0.95, maxit = 50, tol = 1e-04,mu0_g=NULL,sig0_g=NULL){
  kpsi <- switch(psi, bisquare = 1, huber = 2, optimal = 3, median=5,mean=6,hampel=7,8)
  if (kpsi == 8) {
    stop(paste0(psi, " - No such rho function"))
  }
  else if (kpsi==6){
    mu0 = func.med.FM(x,trim=0.1)
    sig0 = sqrt(func.trimvar.FM(x,trim=0.2))
    mu=func.mean(x)
    resu <- list(mu = mu,mu0=func.med.FM(x,trim=0.1),sig0=sig0)
  }
  else if (kpsi==5){

    ktun=keff=1

    if(is.null(mu0_g)) mu0=func.mean_sur(x)
    else mu0=mu0_g
    if(is.null(sig0_g)) sig0=sqrt(func.var_sur(x))
    else sig0=sig0_g

    if (max(sig0) < 1e-10) {
      mu = 0
      sigma = 0
    }
    else {
      n<-dim(x$data)[1]
      st<-standardize_sur(x,mu0,sig0)
      resi_start =norm_fdata_c_sur(st )
      ww_start= wfun_c(resi_start, kpsi,ktun)
      if(all(ww_start==0))mu0=func.trim.FM(x,trim=0.5)
      mu=iteration_sur(x,mu0,sig0,kpsi,ktun,tol, maxit)
    }

    resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

  }
  else {
    family=psi
    kBis = c(3.44, 3.88, 4.685)
    kHub = c(0.732, 0.981, 1.34)
    kHampel<-rbind(c(0.9539156, 1.9078312, 3.8156624),c(NA,NA,NA),c(1.3521, 3.1549, 7.2112))
    kOpt<-c(0.87,0.94,1.06)
    efis = c(0.85, 0.9, 0.95)
    keff = match(eff, efis)
    if(kpsi==7)ktun<-kHampel[keff,]
    else if(kpsi==1)ktun<-kBis[keff]
    else if(kpsi==2)ktun<-kHub[keff]
    else if(kpsi==3)ktun<-kOpt[keff]
    if(is.na(ktun[1])){
      print(c(eff, " No such eff"))
      keff = 0
    }
    else {

      if(is.null(mu0_g)) mu0=func.mean_sur(x)
      else mu0=mu0_g
      if(is.null(sig0_g)) sig0=sqrt(func.var_sur(x))
      else sig0=sig0_g

      if (max(sig0) < 1e-10) {
        mu = 0
        sigma = 0
      }
      else {
        n<-dim(x$data)[1]
        st<-standardize_sur(x,mu0,sig0)
        resi_start =norm_fdata_c_sur(st)
        ww_start= robustbase::Mwgt(resi_start, ktun, family)
        if(all(ww_start==0))sig0=sqrt(func.var_sur(x))
        mu=iteration_ho_sur(x,mu0,sig0,matrix(ktun),family,tol, maxit )

      }
      resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

    }
  }

  return(resu)
}
rfun<-function (x,rho="bisquare",eff=0.95){
  kpsi <- switch(rho, bisquare = 1, huber = 2, optimal = 3,
                 median=5,mean=6,hampel=7,8)
  if (kpsi == 8) {
    stop(paste0(rho, " - No such rho function"))
  }
  else if(kpsi == 5)r=abs(x)
  else if(kpsi == 6)r=x^2
  else{
    kBis = c(3.44, 3.88, 4.685)
    kHub = c(0.732, 0.981, 1.34)
    kHampel<-rbind(c(0.9539156, 1.9078312, 3.8156624),c(NA,NA,NA),c(1.3521, 3.1549, 7.2112))
    kOpt<-c(0.87,0.94,1.06)
    efis = c(0.85, 0.9, 0.95)
    keff = match(eff, efis)
    if(kpsi==7)ktun<-kHampel[keff,]
    else if(kpsi==1)ktun<-kBis[keff]
    else if(kpsi==2)ktun<-kHub[keff]
    else if(kpsi==3)ktun<-kOpt[keff]
    if(is.na(ktun[1])){
      print(c(eff, " No such eff"))
      keff = 0
    }
    else {
      r=robustbase::Mpsi(x, cc = ktun, psi = rho, deriv = -1)
    }
  }
  return(r)
}
scale_fun_pw<-function(x,met="FMAD",...){

  if(met=="FMAD"){
    fun<-function(x,...){

      med<-FlocScaleM(x, psi = "median",...)$mu
      diff<-abs(x-med)
      MAD<-(1/0.675)*pw_median(diff)
      return(MAD)
    }
  }
  else if(met=="VAR"){
    fun<-function(x,...){
      var<-sqrt(func.var(x))
      return(var)
    }
  }
  scale<-fun(x,...)
  return(scale)
}
scale_fun_pw_sur<-function(x,...){

  med<-FlocScaleM_sur(x, psi = "median", ...)$mu
  diff<-abs(standardize_sur(x,med))
  MAD<-(1/0.675)*pw_median_sur(diff)
  return(MAD)

}
pw_median<-function(x){
  data<-x$data
  nvar<-dim(data)[2]
  grid<-x$argvals
  data_new<-sqrt((sapply(1:nvar, function(ii)median(data[,ii])^2)))
  fdata(data_new,argvals = grid)

}
pw_median_sur<-function(x){
  grid<-x$argvals
  data_new<-apply(x$data , c(2,3), median)
  fdata(data_new,argvals = grid)

}
scale_res_oneway_pw<-function(x,label,...){
  k=length(unique(label))
  grid<-x$argval
  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
  res_list<-lapply(1:k,function(ii)x[label==ii]-FlocScaleM(x[label==ii], psi = "median", sig0=sig0,  ...)$mu )
  res_data<-Reduce("rbind",lapply(1:k, function(ii)res_list[[ii]]$data))
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median(abs(res))
  return(med)
}
scale_res_oneway_pw_sur<-function(x,label,...){
  k=length(unique(label))

  grid<-x$argvals
  sig0<-fdata(array(1,dim = c(1,dim(x$data)[2],dim(x$data)[3])),argvals = x$argvals)
  res_list<-lapply(1:k,function(ii)center_sur(ex_fdata(x,label==ii),FlocScaleM_sur(ex_fdata(x,label==ii), psi = "median", sig0=sig0,...)$mu ))
  res_data<-abind::abind(lapply(1:k, function(ii)res_list[[ii]]$data),along = 1)
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median_sur(abs(res))
  return(med)
}
scale_res_twoway_pw<-function(x,label_1,label_2,...){
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
  grid<-x$argvals

  group_mean_ij<-list()
  for (ii in 1:k_1) {
    group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM(x[label_1==ii&label_2==jj],psi = "median", sig0=sig0,...)$mu  )
  }


  kkk=1
  res_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      res_list[[kkk]]<-x[label_1==ii&label_2==jj]-group_mean_ij[[ii]][[jj]]
      kkk=kkk+1
    }
  }
  res_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)res_list[[ii]]$data))
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median(abs(res))
  return(med)

}
scale_res_twoway_pw_sur<-function(x,label_1,label_2,...){
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  grid<-x$argvals
  sig0<-fdata(array(1,dim = c(1,dim(x$data)[2],dim(x$data)[3])),argvals = x$argvals)
  group_mean_ij<-list()
  for (ii in 1:k_1) {
    group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM_sur(ex_fdata(x,label_1==ii&label_2==jj),psi = "median", sig0_g = sig0)$mu  )
  }


  kkk=1
  res_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      res_list[[kkk]]<-center_sur(ex_fdata(x,label_1==ii&label_2==jj),group_mean_ij[[ii]][[jj]])
      kkk=kkk+1
    }
  }
  res_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)res_list[[ii]]$data),along=1)
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median_sur(abs(res))
  return(med)

}

# RoFanova ----------------------------------------------------------------

#' @title Robust Functional Analysis of Variance
#' @description Robust Functional Analysis of Variance (RoFANOVA) allows identifying the presence of significant differences, in terms of
#' functional mean, among groups of a functional data by being robust against the presence of outliers  (Centofanti et al., 2021).
#' @param X Either an object of class  \code{fdata} for monodimensional functional data  or an object of class \code{fdata2d} for bi-dimensional functional data.
#' @param label_1 A vector of containing group label corresponding to the first main effect.
#' @param label_2 A vector of containing group label corresponding to the first main effect. If it is NULL, the one-way RoFANOVA is performed.
#'  Otherwise, the two-way RoFANOVA with interaction is employed. Default is NULL.
#' @param B  The number of permutations used to approximate the p-value in the permutation test. Default is 1000.
#' @param cores If \code{cores}>1, then parallel computing is used, with \code{cores} cores. Default is 1.
#'  @param family The family of loss function for the calcualtion of the equivariant functional M-estimator. The values allowed are
#'  "bisquare" for the bisquare or Tukey's biweight family of loss functions;  "huber" for the the Huber's family of loss functions;
#'    "optimal" for the  optimal family of loss functions; "hampel" for the the Hampel's family of loss functions; "median" for the median loss function.
#'    A non-robust functional estimator of the mean based on the standard least squares loss function is used with the value "mean". Default is "bisquare".
#' @param eff Asymptotic efficiency of the equivariant functional M-estimator. When \code{family} is either "mean" or "median", \code{eff} is ignored.
#' @param mu0_g Initial estimate  used in re-weighted least-squares algorithm to compute the equivariant functional M-estimator.
#' If NULL the standard non-robust functional mean is used. Default is NULL.
#' @param scale Estimate of the standard error of \code{X}. If NUll, the functional normalized median absolute deviation estimator is used. Default is NULL.
#' @param maxit The maximum number of iterations allowed in the re-weighted least-squares algorithm to compute the equivariant functional M-estimator.
#' @param tol The tolerance for the stopping condition of the re-weighted least-squares algorithm to compute the equivariant functional M-estimator.
#' The algorithm stops when the relative variation of the weighted norm sum between two consecutive iterations is less than \code{tol}.
#' @return   A list containing the following arguments:
#' \code{mod} that is a list composed by
#' \itemize{
#' \item \code{data}: A list containing the vectorized form of \code{X}, \code{timeindex}, and \code{curve}. For functional data observed over a regular grid \code{timeindex} and \code{curve} are trivially obtained.
#'
#' \item \code{parameters}: A list containing all the estimated parameters.
#'
#' \item \code{vars}: A list containing results from the Expectation step of the ECM algorithm.
#'
#' \item \code{FullS}: The matrix of B-spline computed over \code{grid}.
#'
#' \item \code{grid}: The vector of time points where the curves are sampled.
#'
#' \item \code{W}: The basis roughness penalty matrix containing the inner products of pairs of basis function second derivatives.
#'
#' \item \code{AW_vec}: Vectorized version of the diagonal matrix used in the approximation of FAPFP.
#'
#' \item \code{P_tot}: Sparse Matrix used to compute all the pairwise comparisons in the FAPFP.
#'
#' \item \code{lambda_s}: Tuning parameter of the smoothness penalty.
#'
#' \item \code{lambda_l}: Tuning parameter of the FAPFP.
#'}
#'
#'A list, named \code{clus}, containing the following arguments:
#'\itemize{
#' \item \code{classes}: The vector of cluster membership.
#'
#' \item \code{po_pr}: Posterior probabilities of cluster membership.
#'}
#'
#'\code{mean_fd} The estimated cluster mean functions.
#'
#'\code{class} A label for the output type.
#'@seealso \code{\link{sasfclust_cv}}
#'
#' @export
#' @references
#' Centofanti, F., Lepore, A., & Palumbo, B. (2021).
#' Sparse and Smooth Functional Data Clustering.
#' \emph{arXiv preprint arXiv:2103.15224}.
#'
#' Ramsay, J., Ramsay, J., & Silverman, B. W. (2005). Functional Data Analysis. Springer Science & Business Media.
#' @examples
#' library(sasfunclust)
#' train<-simulate_data("Scenario I",n_i=20,var_e = 1,var_b = 0.5^2)
#' mod<-sasfclust(X=train$X,grid=train$grid,lambda_s = 10^-6,lambda_l =10,G = 2,maxit = 5,q=10)
#' plot(mod)
#' @importFrom matrixcalc vec
#' @importFrom fda create.bspline.basis fd plot.fd
rofanova<-function(X,label_1,label_2=NULL,B=100,cores=1,family="bisquare",eff=0.95,mu0_g=NULL,scale=NULL,maxit = 50, tol = 1e-04){

  if(length(dim(X_fdata$data))==2){
  if(is.null(label_2)[1]){
    return(rofanova_oneway_perm(X,label_1,B=B,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol,cores=cores))
  }
  else{
    return(rofanova_twoway_perm(X,label_1,label_2,B=B,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol,cores=cores))
  }
  }
  else if(length(dim(X_fdata$data))==3){
    if(is.null(label_2)[1]){
      return(rofanova_oneway_perm_sur(X,label_1,B=B,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol,cores=cores))
    }
    else{
      return(rofanova_twoway_perm_sur(X,label_1,label_2,B=B,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol,cores=cores))
    }

  }


}
rofanova_oneway<-function(X,label,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NULL,mu0_g=NULL){

  k=length(unique(label))
  n=length(label)
  grid<-X$argvals
  if(is.null(scale))scale<-scale_res_oneway_pw(X,label,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g )##median absolute value residual full model#146
  global_mean<-FlocScaleM(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale,mu0_g=mu0_g)$mu
  group_mean<-lapply(1:k,function(ii)FlocScaleM(X[label==ii],psi=family,eff = eff, maxit = maxit, tol = tol,sig0 = scale,mu0_g=mu0_g)$mu  )

  sum_1_i<-stdandar(X,global_mean,scale)
  norm_tot<-norm_fdata_c(sum_1_i)
  sum_1<-sum(sapply(1:n, function(ii)rfun(norm_tot[ii],rho = family,eff = eff)))


  res_list<-lapply(1:k,function(ii)X[label==ii]-group_mean[[ii]])
  red_list<-list()
  for (ll in 1:k) {
    norm_g<-norm_fdata_c(res_list[[ll]]/scale)
    red_list[[ll]]<-sapply(1:length(which(label==ll)), function(ii)rfun(norm_g[ii],rho = family,eff = eff))
  }
  sum_2<-sum(unlist(red_list))
  Tr<-(1/(k-1))*(sum_1-sum_2)

  out<-list(Tr=Tr,
            global_mean=global_mean,
            group_mean=group_mean,
            scale=scale,
            X_fdata=X,
            label=label,
            family=family)
  return(out)

}
rofanova_oneway_perm<-function(X,label,B=100,cores=1,eff=0.95,family="bisquare",mu0_g=NULL,scale=NULL,maxit = 50, tol = 1e-04){
  print("One-way RoFANOVA")
  n<-length(label)

  Tr_obs<-rofanova_oneway(X,label,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr
  perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
  per_fun<-function(kkk){
    perm_comb<-perm_mat[kkk,]
    X_per<-X[perm_comb]
    return(rofanova_oneway(X_per,label,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr)
  }
  if(cores==1){
    per_list<-lapply(1:B,per_fun)
  }
  else{
  if(.Platform$OS.type=="unix")
    per_list<-parallel::mclapply(1:B,per_fun,mc.cores = cores)
  else{
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("X","perm_mat","label","eff","family","mu0_g","scale","maxit", "tol"),envir = environment())
    parallel::clusterEvalQ(cl, library(rofanova))
    per_list<-parallel::parLapply(cl,1:B,per_fun)
    parallel::stopCluster(cl)
  }
}


  Tr_perm<-unlist(per_list)
  pval<-sum(Tr_perm>=Tr_obs)/B
  out<-list(pval=pval,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}
rofanova_twoway<-function(X,label_1,label_2=NULL,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NULL,mu0_g=NULL){



    k_1=length(unique(label_1))
    k_2=length(unique(label_2))
    n=length(label_1)
    grid<-X$argvals
    scale_res<-scale_res_twoway_pw(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g )
    scale_1<-scale_res_oneway_pw(X,label_1,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  )
    scale_2<-scale_res_oneway_pw(X,label_2,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  )
    if(is.null(scale))scale<-scale_fun_pw(X,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g)
    global_mean<-FlocScaleM(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale,mu0_g=mu0_g)$mu
    group_mean_1<-lapply(1:k_1,function(ii)FlocScaleM(X[label_1==ii],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_1,mu0_g=mu0_g)$mu )
    group_mean_2<-lapply(1:k_2,function(ii)FlocScaleM(X[label_2==ii],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_2,mu0_g=mu0_g)$mu  )
    group_mean_ij<-list()
    for (ii in 1:k_1) {
      group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM(X[label_1==ii&label_2==jj],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_res,mu0_g=mu0_g)$mu  )
    }


    kkk=1
    sum_full_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_full_list[[kkk]]<-(X[label_1==ii&label_2==jj]-group_mean_ij[[ii]][[jj]])/scale_res
        kkk=kkk+1
      }
    }
    sum_full_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data))
    sum_full_i<-fdata(sum_full_data,argvals =grid )
    norm_full<-norm_fdata_c(sum_full_i)
    sum_full<-sum(sapply(1:n, function(ii)rfun(norm_full[ii],rho = family,eff = eff)))

    ####FULL
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_red_list[[kkk]]<-(X[label_1==ii&label_2==jj]-global_mean)/scale_res
        kkk=kkk+1
      }
    }
    sum_red_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data))
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_full<-(1/((k_1-1)*(k_2-1)+(k_1-1)+(k_2-1)))*(sum_red-sum_full)


    ####Interaction
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_red_list[[kkk]]<-(X[label_1==ii&label_2==jj]-group_mean_1[[ii]]-group_mean_2[[jj]]+global_mean)/scale_res
        kkk=kkk+1
      }
    }
    sum_red_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data))
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_int<-(1/((k_1-1)*(k_2-1)))*(sum_red-sum_full)
    ####Factor_1
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_red_list[[kkk]]<-(X[label_1==ii&label_2==jj]-global_mean-group_mean_ij[[ii]][[jj]]+group_mean_1[[ii]])/scale_res
        kkk=kkk+1
      }
    }
    sum_red_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data))
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_f1<-(1/((k_1-1)))*(sum_red-sum_full)
    ####Factor_2
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_red_list[[kkk]]<-(X[label_1==ii&label_2==jj]-global_mean-group_mean_ij[[ii]][[jj]]+group_mean_2[[jj]])/scale_res
        kkk=kkk+1
      }
    }
    sum_red_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data))
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_f2<-(1/((k_2-1)))*(sum_red-sum_full)


    Tr_tot<-rbind(Tr_full,Tr_int,Tr_f1,Tr_f2)
    rownames(Tr_tot)<-c("MOD","INT","F1","F2")
    out<-list(Tr_tot=Tr_tot,
              global_mean=global_mean,
              group_mean_1=group_mean_1,
              group_mean_2=group_mean_2,
              group_mean_ij=group_mean_ij,
              scale_res=scale_res,
              X_fdata=X,
              label_1=label_1,
              label_2=label_2,
              family=family)
    return(out)

}
rofanova_twoway_perm<-function(X,label_1,label_2=NULL,B=100,cores=1,eff=0.95,family="bisquare",mu0_g=NULL,scale=NULL,maxit = 50, tol = 1e-04){


    print("Two-way RoFANOVA")
    n<-length(label_1)
    Tr_obs<-rofanova_twoway(X,label_1,label_2,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr_tot


    per_fun<-function(kkk){
      perm_comb<-sample(1:n,replace = F)
      X_per<-X[perm_comb]
      return(rofanova_twoway(X_per,label_1,label_2,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr_tot)
    }
    if(cores==1){
      per_list<-lapply(1:B,per_fun)
    }
    else{
    if(.Platform$OS.type=="unix")
      per_list<-parallel::mclapply(1:B,per_fun,mc.cores = detectCores())
    else{
      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, c("n","X","label_1","label_2","eff","family","mu0_g","scale","maxit", "tol"),envir = environment())
      parallel::clusterEvalQ(cl, library(rofanova))
      per_list<-parallel::parLapply(cl,1:B,per_fun)
      parallel::stopCluster(cl)
    }
    }


    pval_vec<-Tr_obs
    Tr_perm<-list()
    for (ii in 1:nrow(pval_vec)) {
      Tr_perm[[ii]]<-sapply(1:B,function(ll)per_list[[ll]][ii])
      pval_vec[ii]<-sum(Tr_perm[[ii]]>=Tr_obs[ii])/B
    }



  out<-list(pval_vec=pval_vec,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}

# RoFanova Surface----------------------------------------------------------------
rofanova_oneway_sur<-function(X,label,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NULL,mu0_g=NULL){

  k=length(unique(label))
  n=length(label)
  grid<-X$argvals
  if(is.null(scale))scale<-scale_res_oneway_pw_sur(X,label,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g )##median absolute value residual full model#146
  global_mean<-FlocScaleM_sur(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale,mu0_g=mu0_g)$mu
  group_mean<-lapply(1:k,function(ii)FlocScaleM_sur(ex_fdata(X,label==ii),psi=family,eff = eff, maxit = maxit, tol = tol,sig0 = scale,mu0_g=mu0_g)$mu  )

  sum_1_i<-stdandar_sur(X,global_mean,scale)
  norm_tot<-norm_fdata_c_sur(sum_1_i)
  sum_1<-sum(sapply(1:n, function(ii)rfun(norm_tot[ii],rho = family,eff = eff)))


  red_list<-list()
  for (ll in 1:k) {
    red_list[[ll]]<-stdandar_sur(ex_fdata(X,label==ll),group_mean[[ll]],scale)
  }
  sum_red_data<-abind::abind(lapply(1:(k), function(ii)red_list[[ii]]$data),along = 1)
  sum_red_i<-fdata(sum_red_data,argvals =grid )
  norm_red<-norm_fdata_c_sur(sum_red_i)
  sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
  Tr<-(1/((k-1)))*(sum_1-sum_red)

  out<-list(Tr=Tr,
            global_mean=global_mean,
            group_mean=group_mean,
            scale=scale,
            X_fdata=X,
            label=label,
            family=family)
  return(out)

}
rofanova_oneway_perm_sur<-function(X,label,B=100,cores=1,eff=0.95,family="bisquare",mu0_g=NULL,scale=NULL,maxit = 50, tol = 1e-04){
  print("One-way bivariate RoFANOVA")
  n<-length(label)

  Tr_obs<-rofanova_oneway_sur(X,label,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr
  perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
  per_fun<-function(kkk){
    perm_comb<-sample(1:n,replace = F)
    X_per<-ex_fdata(X,perm_comb)
    return(rofanova_oneway_sur(X_per,label,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr)
  }
  if(cores==1){
    per_list<-lapply(1:B,per_fun)
  }
  else{
    if(.Platform$OS.type=="unix")
      per_list<-parallel::mclapply(1:B,per_fun,mc.cores = cores)
    else{
      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, c("n", "X","label","eff","family","mu0_g","scale","maxit", "tol"),envir = environment())
      parallel::clusterEvalQ(cl, library(rofanova))
      per_list<-parallel::parLapply(cl,1:B,per_fun)
      parallel::stopCluster(cl)
    }
  }

  Tr_perm<-unlist(per_list)
  pval<-sum(Tr_perm>=Tr_obs)/B
  out<-list(pval=pval,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}
rofanova_twoway_perm_sur<-function(X,label_1,label_2=NA,B=100,cores=1,eff=0.95,family="bisquare",mu0_g=NULL,scale=NULL,maxit = 50, tol = 1e-04){

  print("Two-way bivariate RoFANOVA")
  n<-length(label_1)
  mod_obs<-rofanova_twoway_sur(X,label_1,label_2,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)
  Tr_obs<-mod_obs$Tr_tot

  per_fun<-function(kkk){
    perm_comb<-sample(1:n,replace = F)
    X_per<-ex_fdata(X,perm_comb)
    return(rofanova_twoway_sur(X_per,label_1,label_2,eff=eff,family=family,mu0_g=mu0_g,scale=scale,maxit = maxit, tol = tol)$Tr_tot)
  }
  if(cores==1){
    per_list<-lapply(1:B,per_fun)
  }
  else{
    if(.Platform$OS.type=="unix")
      per_list<-parallel::mclapply(1:B,per_fun,mc.cores = detectCores())
    else{
      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, c("n","X","label_1","label_2","eff","family","mu0_g","scale","maxit", "tol"),envir = environment())
      parallel::clusterEvalQ(cl, library(rofanova))
      per_list<-parallel::parLapply(cl,1:B,per_fun)
      parallel::stopCluster(cl)
    }
  }

  pval_vec<-Tr_obs
  Tr_perm<-list()
  for (ii in 1:nrow(pval_vec)) {
    Tr_perm[[ii]]<-sapply(1:B,function(ll)per_list[[ll]][ii])
    pval_vec[ii]<-sum(Tr_perm[[ii]]>=Tr_obs[ii])/B
  }


out<-list(pval_vec=pval_vec,
          Tr_obs=Tr_obs,
          Tr_perm=Tr_perm)
return(out)

}
rofanova_twoway_sur<-function(X,label_1,label_2=NA,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NULL,mu0_g=NULL){

    k_1=length(unique(label_1))
    k_2=length(unique(label_2))
    n=length(label_1)
    grid =X$argvals
    scale_res<-scale_res_twoway_pw_sur(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  )#5
    if(is.null(scale))scale<-scale_fun_pw_sur(X,eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  )#5
    scale_1<-lapply(1:k_1,function(ii)scale_fun_pw_sur(ex_fdata(X,label_1==ii),eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  ))#4)
    scale_2<-lapply(1:k_2,function(ii)scale_fun_pw_sur(ex_fdata(X,label_2==ii),eff = eff,tol=tol, maxit = maxit,mu0_g=mu0_g  ))#4
    scale_res_ij<-list()
    for (ii in 1:k_1) {
      scale_res_ij[[ii]]<-lapply(1:k_2,function(jj)scale_fun_pw_sur(ex_fdata(X,label_1==ii&label_2==jj), eff = eff, maxit = maxit, tol =tol)  )
    }
    global_mean<-FlocScaleM_sur(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0=scale,mu0_g=mu0_g)$mu
    group_mean_1<-lapply(1:k_1,function(ii)FlocScaleM_sur(ex_fdata(X,label_1==ii),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_1[[ii]],mu0_g=mu0_g )$mu )
    group_mean_2<-lapply(1:k_2,function(ii)FlocScaleM_sur(ex_fdata(X,label_2==ii),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_2[[ii]],mu0_g=mu0_g )$mu  )
    group_mean_ij<-list()
    for (ii in 1:k_1) {
      group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM_sur(ex_fdata(X,label_1==ii&label_2==jj),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_res_ij[[ii]][[jj]],mu0_g=mu0_g )$mu  )
    }


    kkk=1
    sum_full_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_full_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),group_mean_ij[[ii]][[jj]],scale_res)
        kkk=kkk+1
      }
    }
    sum_full_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
    sum_full_i<-fdata(sum_full_data,argvals =grid )
    norm_full<-norm_fdata_c_sur(sum_full_i)
    sum_full<-sum(sapply(1:n, function(ii)rfun(norm_full[ii],rho = family,eff = eff)))

    ####FULL
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_red_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),global_mean,scale_res)
        kkk=kkk+1
      }
    }
    sum_red_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c_sur(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_full<-(1/((k_1-1)*(k_2-1)+(k_1-1)+(k_2-1)))*(sum_red-sum_full)


    ####Interaction
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        diff_ele<-diff_fdata_sur(sum_fdata_sur(group_mean_1[[ii]],group_mean_2[[jj]]),global_mean)
        sum_red_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),diff_ele,scale_res)
        kkk=kkk+1
      }
    }
    sum_red_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c_sur(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_int<-(1/((k_1-1)*(k_2-1)))*(sum_red-sum_full)
    ####Factor_1
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        diff_ele<-diff_fdata_sur(sum_fdata_sur(global_mean,group_mean_ij[[ii]][[jj]]),group_mean_1[[ii]])
        sum_red_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),diff_ele,scale_res)
        kkk=kkk+1
      }
    }
    sum_red_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c_sur(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_f1<-(1/((k_1-1)))*(sum_red-sum_full)
    ####Factor_2
    kkk=1
    sum_red_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        diff_ele<-diff_fdata_sur(sum_fdata_sur(global_mean,group_mean_ij[[ii]][[jj]]),group_mean_2[[jj]])
        sum_red_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),diff_ele,scale_res)
        kkk=kkk+1
      }
    }
    sum_red_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
    sum_red_i<-fdata(sum_red_data,argvals =grid )
    norm_red<-norm_fdata_c_sur(sum_red_i)
    sum_red<-sum(sapply(1:n, function(ii)rfun(norm_red[ii],rho = family,eff = eff)))
    Tr_f2<-(1/((k_2-1)))*(sum_red-sum_full)


    Tr_tot<-rbind(Tr_full,Tr_int,Tr_f1,Tr_f2)
    rownames(Tr_tot)<-c("MOD","INT","F1","F2")
    out<-list(Tr_tot=Tr_tot,
              global_mean=global_mean,
              group_mean_1=group_mean_1,
              group_mean_2=group_mean_2,
              group_mean_ij=group_mean_ij,
              scale_res=scale_res,
              X_fdata=X,
              label_1=label_1,
              label_2=label_2,
              family=family)
    return(out)
}


# TWOWAY Fanova -----------------------------------------------------------
fanova_twoway<-function(X,label_1,label_2){

  n_ij=length(which(label_1==1&label_2==1))
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  n=length(label_1)

  global_mean<-func.mean(X)
  group_mean_1<-lapply(1:k_1,function(ii)func.mean(X[label_1==ii]))
  group_mean_2<-lapply(1:k_2,function(ii)func.mean(X[label_2==ii]))
  group_mean_ij<-list()
  for (ii in 1:k_1) {
    group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)func.mean(X[label_1==ii&label_2==jj])  )
  }


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      sum_full_list[[kkk]]<-(X[label_1==ii&label_2==jj]-group_mean_ij[[ii]][[jj]])^2
      kkk=kkk+1
    }
  }
  sum_full_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data))
  SSE<-sum_fdata(fdata(sum_full_data,argvals =grid ))
  # Mod ---------------------------------------------------------------------


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      sum_full_list[[kkk]]<-n_ij*(group_mean_ij[[ii]][[jj]]-global_mean)^2
      kkk=kkk+1
    }
  }
  sum_full_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data))
  SSmod<-sum_fdata(fdata(sum_full_data,argvals =grid ))
  T_numden_mod<-int.simpson(SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)),method = "TRAPZ")/int.simpson(SSE/(k_1*k_2*(n_ij-1)),method = "TRAPZ")
  T_full_mod<-int.simpson((SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)))/(SSE/(k_1*k_2*(n_ij-1))),method = "TRAPZ")
  T_mod<-c(numden=T_numden_mod,full=T_full_mod)

  # Int ---------------------------------------------------------------------

  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      sum_full_list[[kkk]]<-n_ij*(group_mean_ij[[ii]][[jj]]-group_mean_1[[ii]]-group_mean_2[[jj]]+global_mean)^2
      kkk=kkk+1
    }
  }
  sum_full_data<-Reduce("rbind",lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data))
  SSint<-sum_fdata(fdata(sum_full_data,argvals =grid ))
  T_numden_int<-int.simpson(SSint/((k_1-1)*(k_2-1)),method = "TRAPZ")/int.simpson(SSE/(k_1*k_2*(n_ij-1)),method = "TRAPZ")
  T_full_int<-int.simpson((SSint/((k_1-1)*(k_2-1)))/(SSE/(k_1*k_2*(n_ij-1))),method = "TRAPZ")
  T_int<-c(numden=T_numden_int,full=T_full_int)
  # F1 ----------------------------------------------------------------------
  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    sum_full_list[[kkk]]<-n_ij*k_2*(group_mean_1[[ii]]-global_mean)^2
    kkk=kkk+1

  }
  sum_full_data<-Reduce("rbind",lapply(1:(k_1), function(ii)sum_full_list[[ii]]$data))
  SS1<-sum_fdata(fdata(sum_full_data,argvals =grid ))
  T_numden_1<-int.simpson(SS1/((k_1-1)),method = "TRAPZ")/int.simpson(SSE/(k_1*k_2*(n_ij-1)),method = "TRAPZ")
  T_full_1<-int.simpson((SS1/((k_1-1)))/(SSE/(k_1*k_2*(n_ij-1))),method = "TRAPZ")
  T_1<-c(numden=T_numden_1,full=T_full_1)
  # F2 ----------------------------------------------------------------------

  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_2) {
    sum_full_list[[kkk]]<-n_ij*k_1*(group_mean_2[[ii]]-global_mean)^2
    kkk=kkk+1

  }
  sum_full_data<-Reduce("rbind",lapply(1:(k_2), function(ii)sum_full_list[[ii]]$data))
  SS2<-sum_fdata(fdata(sum_full_data,argvals =grid ))
  T_numden_2<-int.simpson(SS2/((k_2-1)),method = "TRAPZ")/int.simpson(SSE/(k_1*k_2*(n_ij-1)),method = "TRAPZ")
  T_full_2<-int.simpson((SS2/((k_2-1)))/(SSE/(k_1*k_2*(n_ij-1))),method = "TRAPZ")
  T_2<-c(numden=T_numden_2,full=T_full_2)


  Tr_tot<-rbind(T_mod,T_int,T_1,T_2)
  rownames(Tr_tot)<-c("MOD","INT","F1","F2")
  out<-list(Tr_tot=Tr_tot,
            global_mean=global_mean,
            group_mean_1=group_mean_1,
            group_mean_2=group_mean_2,
            group_mean_ij=group_mean_ij,
            X_fdata=X,
            label_1=label_1,
            label_2=label_2)
  return(out)
}
fanova_twoway_perm<-function(X,label_1,label_2,B=100,...){

  n<-length(label_1)

  Tr_obs<-fanova_twoway(X,label_1,label_2)$Tr_tot
  perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
  per_fun<-function(kkk){
    perm_comb<-perm_mat[kkk,]
    X_per<-X[perm_comb]
    return(fanova_twoway(X_per,label_1,label_2)$Tr_tot)
  }
  per_list<-mclapply(1:B,per_fun,mc.cores = detectCores())


  pval_mat<-Tr_obs
  Tr_perm<-list()
  for (ii in 1:nrow(pval_mat)) {
    Tr_perm[[ii]]<-list()
    for (jj in 1:ncol(pval_mat)) {
      Tr_perm[[ii]][[jj]]<-sapply(1:B,function(ll)per_list[[ll]][ii,jj])
      pval_mat[ii,jj]<-sum(Tr_perm[[ii]][[jj]]>=Tr_obs[ii,jj])/B
    }

  }

  out<-list(pval_mat=pval_mat,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}
fanova_twoway_sur<-function(X,label_1,label_2){

  n_ij=length(which(label_1==1&label_2==1))
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  n=length(label_1)
  grid=X$argvals
  global_mean<-func.mean_sur(X)
  group_mean_1<-lapply(1:k_1,function(ii)func.mean_sur(ex_fdata(X,label_1==ii)))
  group_mean_2<-lapply(1:k_2,function(ii)func.mean_sur(ex_fdata(X,label_2==ii)))
  group_mean_ij<-list()
  for (ii in 1:k_1) {
    group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)func.mean_sur(ex_fdata(X,label_1==ii&label_2==jj)))
  }


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      sum_full_list[[kkk]]<-center_sur(ex_fdata(X,label_1==ii&label_2==jj),group_mean_ij[[ii]][[jj]])^2
      kkk=kkk+1
    }
  }
  sum_full_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
  SSE<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))

  # Mod ---------------------------------------------------------------------


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      sum_full_list[[kkk]]<-n_ij*diff_fdata_sur(group_mean_ij[[ii]][[jj]],global_mean)^2
      kkk=kkk+1
    }
  }
  sum_full_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
  SSmod<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))
  T_numden_mod<-int_sur(SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)))/int_sur(SSE/(k_1*k_2*(n_ij-1)))
  T_full_mod<-int_sur(frac_fdata_sur(SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)),SSE/(k_1*k_2*(n_ij-1))))
  T_mod<-c(numden=T_numden_mod,full=T_full_mod)

  # Int ---------------------------------------------------------------------


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      ele<-sum_fdata_sur(diff_fdata_sur(group_mean_ij[[ii]][[jj]],sum_fdata_sur(group_mean_1[[ii]],group_mean_2[[jj]])),global_mean)
      sum_full_list[[kkk]]<-n_ij*(ele)^2
      kkk=kkk+1
    }
  }
  sum_full_data<-abind::abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
  SSint<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))
  T_numden_int<-int_sur(SSint/((k_1-1)*(k_2-1)))/int_sur(SSE/(k_1*k_2*(n_ij-1)))
  T_full_int<-int_sur(frac_fdata_sur(SSint/((k_1-1)*(k_2-1)),(SSE/(k_1*k_2*(n_ij-1)))))
  T_int<-c(numden=T_numden_int,full=T_full_int)
  # F1 ----------------------------------------------------------------------

  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    sum_full_list[[kkk]]<-n_ij*k_2*diff_fdata_sur(group_mean_1[[ii]],global_mean)^2
    kkk=kkk+1

  }
  sum_full_data<-abind::abind(lapply(1:(k_1), function(ii)sum_full_list[[ii]]$data),along = 1)
  SS1<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))
  T_numden_1<-int_sur(SS1/((k_1-1)))/int_sur(SSE/(k_1*k_2*(n_ij-1)))
  T_full_1<-int_sur(frac_fdata_sur(SS1/((k_1-1)),(SSE/(k_1*k_2*(n_ij-1)))))
  T_1<-c(numden=T_numden_1,full=T_full_1)
  # F2 ----------------------------------------------------------------------

  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_2) {
    sum_full_list[[kkk]]<-n_ij*k_1*diff_fdata_sur(group_mean_2[[ii]],global_mean)^2
    kkk=kkk+1

  }
  sum_full_data<-abind::abind(lapply(1:(k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
  SS2<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))
  T_numden_2<-int_sur(SS2/((k_2-1)))/int_sur(SSE/(k_1*k_2*(n_ij-1)))
  T_full_2<-int_sur(frac_fdata_sur(SS2/(k_2-1),(SSE/(k_1*k_2*(n_ij-1)))))
  T_2<-c(numden=T_numden_2,full=T_full_2)


  Tr_tot<-rbind(T_mod,T_int,T_1,T_2)
  rownames(Tr_tot)<-c("MOD","INT","F1","F2")
  out<-list(Tr_tot=Tr_tot,
            global_mean=global_mean,
            group_mean_1=group_mean_1,
            group_mean_2=group_mean_2,
            group_mean_ij=group_mean_ij,
            SSE=SSE,
            X_fdata=X,
            label_1=label_1,
            label_2=label_2)
  return(out)
}
fanova_twoway_perm_sur<-function(X,label_1,label_2,B=100,...){

  n<-length(label_1)

  Tr_obs<-fanova_twoway_sur(X,label_1,label_2)$Tr_tot
  perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
  per_fun<-function(kkk){
    perm_comb<-perm_mat[kkk,]
    X_per<-ex_fdata(X,perm_comb)
    return(fanova_twoway_sur(X_per,label_1,label_2)$Tr_tot)
  }
  per_list<-mclapply(1:B,per_fun,mc.cores = detectCores())


  pval_mat<-Tr_obs
  Tr_perm<-list()
  for (ii in 1:nrow(pval_mat)) {
    Tr_perm[[ii]]<-list()
    for (jj in 1:ncol(pval_mat)) {
      Tr_perm[[ii]][[jj]]<-sapply(1:B,function(ll)per_list[[ll]][ii,jj])
      pval_mat[ii,jj]<-sum(Tr_perm[[ii]][[jj]]>=Tr_obs[ii,jj])/B
    }

  }

  out<-list(pval_mat=pval_mat,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}
