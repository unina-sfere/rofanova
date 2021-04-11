
# Generate data -----------------------------------------------------------

simulate_data<-function(mean="M1",con="C1",p=0.1,M=1,n_i=10,k=3,sd=0.01,grid=seq(0,1,length.out = 30),c=1){


  if(mean=="M1"){
    mean_function<-function(t,i)t*(1-t)
  }
  else if(mean=="M2"){
    mean_function<-function(t,i)(t^(i))*(1-t)^(6-i)
  }
  else if(mean=="M3"){
    mean_function<-function(t,i)t^(i/5)*(1-t)^((6-i)/5)
  }
  else if(mean=="M4"){
    mean_function<-function(t,i)0.1*abs(sin(4*pi*t))
  }

  else if(mean=="M5"){
    mean_function<-function(t,i)i*0.05*abs(sin(4*pi*t))
  }
  else if(mean=="M6"){
    mean_function<-function(t,i)i*0.025*abs(sin(4*pi*t))
  }
  else if(mean=="M7"){
    mean_function<-function(t,i)1+i/50
  }

  if(con=="C0"){
    cont_function<-function(n_g,p_g,M_g,ll)0
  }
  else if(con=="C1"){
    cont_function<-function(n_g,p_g,M_g,ll)M_g*matrix(rbinom(n_g,1,p_g)*sample(c(-1,1),n_g,replace = T),length(grid),n_g,byrow = T)
  }
  else if(con=="C2"){

    cont_function<-function(n_g,p_g,M_g,ll){
      T_i<-runif(n_g,0,1)
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
    cont_function<-function(n_g,p_g,M_g,ll){
      g<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      Y<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.001))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*Y
      return(out)
    }
  }
  else if(con=="C4"){
    cont_function<-function(n_g,p_g,M_g,ll)(-1)^(ll)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
  }
  else if(con=="C5"){
    cont_function<-function(n_g,p_g,M_g,ll){
      T_i<-runif(n_g,0,1)
      matrix_old<-(-1)^(ll)*M_g*matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      matrix_new<-matrix(0,length(matrix_old[,1]),n_g)
      for (ii in 1:n_g) {
        ind<-(round(T_i[ii]*(length(matrix_old[,1])-1))+1)
        matrix_new[ind:length(matrix_old[,1]),ii]<-matrix_old[ind:length(matrix_old[,1]),ii]

      }
      return(matrix_new)
    }
  }
  else if(con=="C6"){
    cont_function<-function(n_g,p_g,M_g,ll){
      g<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.00001))$data)
      Y<-mean_function(matrix(grid,length(grid),n_i),matrix(ll,length(grid),n_i))+t(rproc2fdata(n_i,t=grid,sigma="vexponential",par.list = list(scale=sd^2,theta=0.0001*10^(ll)))$data)
      ber<-matrix(rbinom(n_g,1,p_g),length(grid),n_g,byrow = T)
      out<-(1-ber)*g+ber*Y
      return(out)
    }
  }

  else if(con=="C7"){
    cont_function<-function(t,i,ll)0.1*abs(sin(4*pi*t))
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
  label<-rep(1:k,each=n_i)
  out=list(data=data,
           label=label,
           grid=grid)
  return(out)
}

simulate_data_twoway<-function(mean="M1",con="C1",n_i=10,k_1=3,k_2=3,p=0.1,M=1,alpha=0.05,beta=0.05,sd=0.01,grid=seq(0,1,length.out = 30),c=1){



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
  label_1<-rep(1:k_1,each=n_i*k_2)
  label_2<-rep(rep(1:k_2,each=n_i),k_1)
  out=list(data=data,
           label_1=label_1,
           label_2=label_2,
           grid=grid)
  return(out)
}

# Robust locatio/scale estimation -----------------------------------------


FlocScaleM<-function (x, psi = "bisquare", eff = 0.95, maxit = 50, tol = 1e-04,mu0_g=NA,sig0_g=NA,cpp=TRUE)
{
  kpsi <- switch(psi, bisquare = 1, huber = 2, optimal = 3,
                 modopt = 4, median=5,mean=6,hampel=7,8)
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
      if(cpp==TRUE){
        mu=iteration(x,mu0,sig0,kpsi,ktun,tol, maxit)
      }
      else{
        dife = 1e+10
        iter = 0
        while (dife > tol & iter < maxit) {
          iter = iter + 1
          resi =norm_fdata_c((x-mu0)/sig0 )
          ww = wfun(resi, kpsi,ktun)
          prod<-x
          prod$data<-diag(as.numeric(ww))%*%x$data
          mu = sum_fdata(prod)/sum(ww)
          resi_new =norm_fdata((x-mu)/sig0)
          dife = abs(mean(resi_new) - mean(resi))/mean(resi)
          mu0 = mu
        }



      }
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
        ww_start= Mwgt(resi_start, ktun, family)
        if(all(ww_start==0))mu0=func.trim.FM(x,trim=0.5)
        if(cpp==TRUE){
          mu=iteration_ho(x,mu0,sig0,matrix(ktun),family,tol, maxit)
        }
        else{
          dife <- 1e+10
          iter <- 0
          while (dife > tol && iter < maxit) {
            cat(iter)
            iter <- iter + 1
            st<-stdandar(x, mu0,sig0)
            resi <- norm_fdata_c(st)
            ww <- Mwgt(resi, ktun, family)
            prod<-x
            prod$data<-diag(as.numeric(ww))%*%x$data
            mu = sum_fdata_c(prod)/sum(ww)
            st2<-stdandar(x, mu,sig0)
            resi_new =norm_fdata_c(st2)
            dife = abs(mean(resi_new) - mean(resi))/mean(resi)
            mu0 <- mu
          }
        }
      }
      resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

    }
  }

  return(resu)
}
FlocScaleM_sur<-function (x, psi = "bisquare", eff = 0.95, maxit = 50, tol = 1e-04,mu0_g=NA,sig0_g=NA,cpp=TRUE)
{
  kpsi <- switch(psi, bisquare = 1, huber = 2, optimal = 3,
                 modopt = 4, median=5,mean=6,hampel=7,8)
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

    if(is.na(mu0_g)) mu0=func.mean_sur(x)
    else mu0=mu0_g
    if(is.na(sig0_g)) sig0=sqrt(func.var_sur(x))
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
      if(cpp==TRUE){
        mu=iteration_sur(x,mu0,sig0,kpsi,ktun,tol, maxit)
      }
      else{
        dife = 1e+10
        iter = 0
        while (dife > tol & iter < maxit) {
          cat(iter)
          iter = iter + 1
          st<-stdandar_sur(x,mu0,sig0)
          resi =norm_fdata_c_sur(st )
          ww = wfun(resi, kpsi,ktun)
          prod<-x
          prod$data<-as.numeric(ww)*x$data
          mu = acc_fdata_sur(prod)/sum(ww)
          st_new<-stdandar_sur(x,mu,sig0)
          resi_new =norm_fdata_c_sur(st_new)
          dife = abs(mean(resi_new) - mean(resi))/mean(resi)
          mu0 = mu
        }



      }
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

      if(is.na(mu0_g)) mu0=func.mean_sur(x)
      else mu0=mu0_g
      if(is.na(sig0_g)) sig0=sqrt(func.var_sur(x))
      else sig0=sig0_g

      if (max(sig0) < 1e-10) {
        mu = 0
        sigma = 0
      }
      else {
        n<-dim(x$data)[1]
        st<-standardize_sur(x,mu0,sig0)
        resi_start =norm_fdata_c_sur(st)
        ww_start= Mwgt(resi_start, ktun, family)
        if(all(ww_start==0))sig0=sqrt(func.var_sur(x))
        if(cpp==TRUE){
          mu=iteration_ho_sur(x,mu0,sig0,matrix(ktun),family,tol, maxit )
        }
        else{
          dife <- 1e+10
          iter <- 0
          while (dife > tol && iter < maxit) {
            cat(iter)
            iter <- iter + 1
            st<-standardize_sur(x,mu0,sig0)
            resi =norm_fdata_c_sur(st )
            ww <- Mwgt(resi, ktun, family)
            prod<-x
            for(ii in 1:n) prod$data[ii,,]<-ww[ii]*x$data[ii,,]
            mu = acc_fdata_sur(prod)/sum(ww)
            st_new<-standardize_sur(x,mu,sig0)
            resi_new =norm_fdata_c_sur(st_new)
            dife = abs(mean(resi_new) - mean(resi))/mean(resi)
            mu0 <- mu
          }
        }
      }
      resu <- list(mu = mu,mu0=mu0_g,sig0=sig0_g)

    }
  }

  return(resu)
}
wfun<-function (x, k,ktun)
{
  if (k == 1)  ww = (1 - (x/ktun)^2)^2 * (abs(x) <= ktun)
  else if(k==2) ww = (abs(x) <= ktun) + (abs(x) > ktun)*ktun/(abs(x) + 1e-20)
  else if(k==5) ww = 1/(abs(x) + 1e-10)

  return(ww)
}
wfun_der<-function (x, k,ktun)
{
  if (k == 1)
    ww = (1 - (x/ktun)^2)^2 * (abs(x) <= ktun)
  else if(k==2) ww = (abs(x) <= ktun) + (abs(x) > ktun)*ktun/(abs(x) + 1e-20)
  else if(k==5) ww = 1/(abs(x) + 1e-10)
  return(ww)
}
rfun<-function (x,rho="bisquare",eff=0.95)
{
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
      r=Mpsi(x, cc = ktun, psi = rho, deriv = -1)
    }
  }
  return(r)
}
rfun_pw<-function (x, rho="bisquare",eff=0.95)
{
  out<-x
  data<-x$data
  for (ii in 1:dim(data)[1]) {
    for (jj in 1:dim(data)[2]) {
      data[ii,jj]<-rfun(x$data[ii,jj],rho = rho,eff)
    }

  }
  out$data<-data
  return(out)
}
psif_der_pw<-function (x, rho="bisquare",eff=0.9)
{
  out<-x
  data<-x$data
  for (ii in 1:dim(data)[1]) {
    for (jj in 1:dim(data)[2]) {
      data[ii,jj]<-psif_der(x$data[ii,jj],rho,eff)
    }

  }
  out$data<-data
  return(out)
}
psif_pw<-function (x, rho="bisquare",eff=0.9)
{
  out<-x
  data<-x$data
  for (ii in 1:dim(data)[1]) {
    for (jj in 1:dim(data)[2]) {
      data[ii,jj]<-psif_2(x$data[ii,jj],rho,eff)
    }

  }
  out$data<-data
  return(out)


  scale<-fun(x)
  return}
psif<-function (x, k)
{
  return(x * wfun(x, k))
}
psif_2<-function (x, rho,eff=0.9)
{

  kpsi <- switch(rho, bisquare = 1, huber = 2, optimal = 3,
                 modopt = 4, median=5,mean=6, 7)

  if(kpsi==5)kpsi=3
  kBis = c(3.44, 3.88, 4.685)
  kHub = c(0.732, 0.981, 1.34)
  kmedian=c(1,1,1)
  kk = rbind(kBis, kHub,kmedian)
  efis = c(0.85, 0.9, 0.95)
  keff = match(eff, efis)
  ktun = kk[kpsi, keff]
  return(x * wfun(x, kpsi,ktun))
}
psif_der<-function (x, rho,eff=0.9)
{

  kpsi <- switch(rho, bisquare = 1, huber = 2, optimal = 3,
                 modopt = 4, median=5,mean=6, 7)
  if(kpsi==5)kpsi=3
  kBis = c(3.44, 3.88, 4.685)
  kHub = c(0.732, 0.981, 1.34)
  kmedian=c(1,1,1)
  kk = rbind(kBis, kHub,kmedian)
  efis = c(0.85, 0.9, 0.95)
  keff = match(eff, efis)
  ktun = kk[kpsi, keff]
  return(x * wfun_der(x, kpsi,ktun))
}
scale_fun<-function(x,met="FMAD",...){

  if(met=="FMAD"){
    fun<-function(x,...){

      med<-FlocScaleM(x, psi = "median", eff = 0.95, maxit = 100, tol = 1e-20,...)$mu
      diff<-abs(x-med)
      MAD<-FlocScaleM(diff, psi = "median", eff = 0.95, maxit = 100, tol = 1e-20,...)$mu
      return(MAD)
    }
  }
  else if(met=="VAR"){
    fun<-function(x,...){
      var<-sqrt(func.var(x))
      return(var)
    }
  }
  scale<-fun(x)
  return(scale)
}
scale_fun_pw<-function(x,met="FMAD",...){

  if(met=="FMAD"){
    fun<-function(x,...){

      med<-FlocScaleM(x, psi = "median", eff = 0.95, maxit = 100, tol = 1e-20,...)$mu
      diff<-abs(x-med)
      MAD<-pw_median(diff)
      return(MAD)
    }
  }
  else if(met=="VAR"){
    fun<-function(x,...){
      var<-sqrt(func.var(x))
      return(var)
    }
  }
  scale<-fun(x)
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
scale_res_oneway<-function(x,label,...){
  k=length(unique(label))

  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
  res_list<-lapply(1:k,function(ii)x[label==ii]-FlocScaleM(x[label==ii], psi = "median", sig0=sig0,  ...)$mu )
  res_data<-Reduce("rbind",lapply(1:k, function(ii)res_list[[ii]]$data))
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.8)*FlocScaleM(abs(res), psi = "median",sig0=sig0,  ...)$mu
  return(med)
}
scale_res_twoway<-function(x,label_1,label_2,...){
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
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
  med<-(1/0.8)*FlocScaleM(abs(res), psi = "median",  sig0=sig0,...)$mu

  return(med)

}
scale_res_oneway_pw<-function(x,label,...){
  k=length(unique(label))

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
  res_data<-abind(lapply(1:k, function(ii)res_list[[ii]]$data),along = 1)
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median_sur(abs(res))
  return(med)
}
scale_res_twoway_pw<-function(x,label_1,label_2,...){
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
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
  res_data<-abind(lapply(1:(k_1*k_2), function(ii)res_list[[ii]]$data),along=1)
  res<-fdata(res_data,argvals =grid )
  med<-(1/0.675)*pw_median_sur(abs(res))
  return(med)

}
scale_res_oneway_trim<-function(x,label,trim=0.2){
  k=length(unique(label))


  data_list<-lapply(1:k,function(ii)x[label==ii]-func.trim.FM(x[label==ii],trim=trim) )
  data<-Reduce("rbind",lapply(1:k, function(ii)data_list[[ii]]$data))
  fdata_i<-fdata(data,argvals =grid )

  med<-sqrt(func.trimvar.FM(fdata_i,trim=trim))

  return(med)
}
scale_res_twoway_trim<-function(x,label_1,label_2,trim=0.2){
  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  sig0<-fdata(rep(1,length(x$data[1,])),argvals = x$argvals)
  group_mean_ij<-list()
  for (ii in 1:k_1) {
    group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)func.trim.FM(x[label_1==ii&label_2==jj],trim=trim))
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
  med<-sqrt(func.trimvar.FM(res,trim=trim))

  return(med)

}
norm_fdata<-function(x){
  sqrt(int.simpson(x^2))
}

# RoFanova ----------------------------------------------------------------

RoFanova<-function(X,label,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NA,...){

  k=length(unique(label))
  n=length(label)
  if(is.na(scale))scale<-scale_res_oneway_pw(X,label,eff = eff,tol=tol, maxit = maxit,... )##median absolute value residual full model#146
  global_mean<-FlocScaleM(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale,...)$mu
  group_mean<-lapply(1:k,function(ii)FlocScaleM(X[label==ii],psi=family,eff = eff, maxit = maxit, tol = tol,sig0 = scale,...)$mu  )

  sum_1_i<-stdandar(X,global_mean,scale)
  norm_tot<-norm_fdata_c(sum_1_i)
  sum_1<-sum(sapply(1:n, function(ii)rfun(norm_tot[ii],rho = family,eff = eff)))


  res_list<-lapply(1:k,function(ii)X[label==ii]-group_mean[[ii]])
  sss<-list()
  for (ll in 1:k) {
    norm_g<-norm_fdata_c(res_list[[ll]]/scale)
    sss[[ll]]<-sapply(1:n_i, function(ii)rfun(norm_g[ii],rho = family,eff = eff))
  }
  sum_2<-sum(unlist(sss))
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

#' @export
RoFanova_perm<-function(X,label,B=100,...){

  n<-length(label)

  Tr_obs<-RoFanova(X,label,...)$Tr
  perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
  per_fun<-function(kkk){
    perm_comb<-perm_mat[kkk,]
    X_per<-X[perm_comb]
    return(RoFanova(X_per,label,...)$Tr)
  }
  per_list<-lapply(1:B,per_fun)#,mc.cores = 1)
  Tr_perm<-unlist(per_list)
  pval<-sum(Tr_perm>=Tr_obs)/B
  out<-list(pval=pval,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}
RoFanova_twoway<-function(X,label_1,label_2=NA,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale_glob=NA,...){

  if(is.na(label_2[1])){
    return(RoFanova(X,label_1,family=family,eff = eff, maxit = maxit, tol =tol,...))
  }
  else{

    k_1=length(unique(label_1))
    k_2=length(unique(label_2))
    n=length(label_1)
    # scale_res<-scale_res_twoway(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )
    # scale_1<-scale_res_oneway(X,label_1,eff = eff,tol=tol, maxit = maxit,... )
    # scale_2<-scale_res_oneway(X,label_2,eff = eff,tol=tol, maxit = maxit,... )
    # if(is.na(scale_glob))scale_glob<-scale_fun(X,eff = eff,tol=tol, maxit = maxit,... )
    scale_res<-scale_res_twoway_pw(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )
    scale_1<-scale_res_oneway_pw(X,label_1,eff = eff,tol=tol, maxit = maxit,... )
    scale_2<-scale_res_oneway_pw(X,label_2,eff = eff,tol=tol, maxit = maxit,... )
    if(is.na(scale_glob))scale_glob<-scale_fun_pw(X,eff = eff,tol=tol, maxit = maxit,... )
    # scale_res<-scale_res_twoway(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit )
    # scale_1<-fdata(rep(1,length(X$data[1,])),argvals =X$argvals)
    # scale_2<-fdata(rep(1,length(X$data[1,])),argvals =X$argvals)
    # scale_glob<-fdata(rep(1,length(X$data[1,])),argvals =X$argvals)
    # scale_res<-scale_res_twoway_trim(X,label_1,label_2,trim=0.2)
    # scale_1<-scale_res_oneway_trim(X,label_1,trim=0.2 )
    # scale_2<-scale_res_oneway_trim(X,label_2,trim=0.2 )
    # scale_glob<-sqrt(func.trimvar.FM(X,trim=0.2))

    global_mean<-FlocScaleM(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_glob,...)$mu
    group_mean_1<-lapply(1:k_1,function(ii)FlocScaleM(X[label_1==ii],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_1,...)$mu )
    group_mean_2<-lapply(1:k_2,function(ii)FlocScaleM(X[label_2==ii],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_2,...)$mu  )
    group_mean_ij<-list()
    for (ii in 1:k_1) {
      group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM(X[label_1==ii&label_2==jj],psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_res,...)$mu  )
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
}
RoFanova_twoway_perm<-function(X,label_1,label_2=NA,B=100,...){

  if(is.na(label_2[1])){
    cat("one-way")
    return(RoFanova_perm(X,label_1,B=B,...))
  }
  else{
    n<-length(label_1)
    Tr_obs<-RoFanova_twoway(X,label_1,label_2,...)$Tr_tot


    per_fun<-function(kkk){
      perm_comb<-sample(1:n,replace = F)
      X_per<-X[perm_comb]
      return(RoFanova_twoway(X_per,label_1,label_2,...)$Tr_tot)
    }
    per_list<-mclapply(1:B,per_fun,mc.cores = detectCores())


    pval_vec<-Tr_obs
    Tr_perm<-list()
    for (ii in 1:nrow(pval_vec)) {
      Tr_perm[[ii]]<-sapply(1:B,function(ll)per_list[[ll]][ii])
      pval_vec[ii]<-sum(Tr_perm[[ii]]>=Tr_obs[ii])/B
    }

  }

  out<-list(pval_vec=pval_vec,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
  return(out)

}

# RoFanova ----------------------------------------------------------------

# RoFanova<-function(X,label,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale=NA,...){
#
#   k=length(unique(label))
#   n=length(label)
#   if(is.na(scale))scale<-scale_res_oneway_pw(X,label,eff = eff,tol=tol, maxit = maxit,... )##median absolute value residual full model#146
#   global_mean<-FlocScaleM(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale,...)$mu
#   group_mean<-lapply(1:k,function(ii)FlocScaleM(X[label==ii],psi=family,eff = eff, maxit = maxit, tol = tol,sig0 = scale,...)$mu  )
#
#   sum_1_i<-stdandar(X,global_mean,scale)
#   norm_tot<-norm_fdata_c(sum_1_i)
#   sum_1<-sum(sapply(1:n, function(ii)rfun(norm_tot[ii],rho = family,eff = eff)))
#
#
#   res_list<-lapply(1:k,function(ii)X[label==ii]-group_mean[[ii]])
#   sss<-list()
#   for (ll in 1:k) {
#     norm_g<-norm_fdata_c(res_list[[ll]]/scale)
#     sss[[ll]]<-sapply(1:n_i, function(ii)rfun(norm_g[ii],rho = family,eff = eff))
#   }
#   sum_2<-sum(unlist(sss))
#   Tr<-(1/(k-1))*(sum_1-sum_2)
#
#   out<-list(Tr=Tr,
#             global_mean=global_mean,
#             group_mean=group_mean,
#             scale=scale,
#             X_fdata=X,
#             label=label,
#             family=family)
#   return(out)
#
# }
# RoFanova_perm<-function(X,label,B=100,...){
#
#   n<-length(label)
#
#   Tr_obs<-RoFanova(X,label,...)$Tr
#   perm_mat<-t(sapply(1:B,function(ii)sample(1:n,replace = F)))
#   per_fun<-function(kkk){
#     perm_comb<-perm_mat[kkk,]
#     X_per<-X[perm_comb]
#     return(RoFanova(X_per,label,...)$Tr)
#   }
#   per_list<-mclapply(1:B,per_fun,mc.cores = detectCores())
#   Tr_perm<-unlist(per_list)
#   pval<-sum(Tr_perm>=Tr_obs)/B
#   out<-list(pval=pval,
#             Tr_obs=Tr_obs,
#             Tr_perm=Tr_perm)
#   return(out)
#
# }
RoFanova_twoway_sur<-function(X,label_1,label_2=NA,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,scale_glob=NA,...){

  if(is.na(label_2[1])){
    return(RoFanova(X,label_1,family=family,eff = eff, maxit = maxit, tol =tol,...))
  }
  else{

    k_1=length(unique(label_1))
    k_2=length(unique(label_2))
    n=length(label_1)
    grid =X$argvals
    # scale_res<-scale_res_twoway(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )
    # scale_1<-scale_res_oneway(X,label_1,eff = eff,tol=tol, maxit = maxit,... )
    # scale_2<-scale_res_oneway(X,label_2,eff = eff,tol=tol, maxit = maxit,... )
    # if(is.na(scale_glob))scale_glob<-scale_fun(X,eff = eff,tol=tol, maxit = maxit,... )
    scale_res<-scale_res_twoway_pw_sur(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )#5
    # scale_1<-scale_res_oneway_pw_sur(X,label_1,eff = eff,tol=tol, maxit = maxit,... )#4
    # scale_2<-scale_res_oneway_pw_sur(X,label_2,eff = eff,tol=tol, maxit = maxit,...)#4
    if(is.na(scale_glob))scale_glob<-scale_fun_pw_sur(X,eff = eff,tol=tol, maxit = maxit,... )#5
    # scale_res<-scale_res_twoway_pw_sur(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )
    # scale_1<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
    # scale_2<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
    # scale_glob<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
    # scale_res<-scale_res_twoway_trim(X,label_1,label_2,trim=0.2)
    # scale_1<-scale_res_oneway_trim(X,label_1,trim=0.2 )
    # scale_2<-scale_res_oneway_trim(X,label_2,trim=0.2 )
    # scale_glob<-sqrt(func.trimvar.FM(X,trim=0.2))
    # scale_res<-scale_res_twoway_pw_sur(X,label_1,label_2,eff = eff,tol=tol, maxit = maxit,... )#5
    scale_1<-lapply(1:k_1,function(ii)scale_fun_pw_sur(ex_fdata(X,label_1==ii),eff = eff,tol=tol, maxit = maxit,... ))#4)
    scale_2<-lapply(1:k_2,function(ii)scale_fun_pw_sur(ex_fdata(X,label_2==ii),eff = eff,tol=tol, maxit = maxit,... ))#4
    # # if(is.na(scale_glob))scale_glob<-scale_fun_pw_sur(X,eff = eff,tol=tol, maxit = maxit,... )#5
    scale_res_ij<-list()
    for (ii in 1:k_1) {
      scale_res_ij[[ii]]<-lapply(1:k_2,function(jj)scale_fun_pw_sur(ex_fdata(X,label_1==ii&label_2==jj), eff = eff, maxit = maxit, tol =tol)  )
    }
    global_mean<-FlocScaleM_sur(X,psi = family, eff = eff, maxit = maxit, tol =tol,sig0=scale_glob,...)$mu
    group_mean_1<-lapply(1:k_1,function(ii)FlocScaleM_sur(ex_fdata(X,label_1==ii),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_1[[ii]],...)$mu )
    group_mean_2<-lapply(1:k_2,function(ii)FlocScaleM_sur(ex_fdata(X,label_2==ii),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_2[[ii]],...)$mu  )
    group_mean_ij<-list()
    for (ii in 1:k_1) {
      group_mean_ij[[ii]]<-lapply(1:k_2,function(jj)FlocScaleM_sur(ex_fdata(X,label_1==ii&label_2==jj),psi = family, eff = eff, maxit = maxit, tol =tol,sig0 = scale_res_ij[[ii]][[jj]],...)$mu  )
    }


    kkk=1
    sum_full_list<-list()
    for (ii in 1:k_1) {
      for (jj in 1:k_2) {
        sum_full_list[[kkk]]<-stdandar_sur(ex_fdata(X,label_1==ii&label_2==jj),group_mean_ij[[ii]][[jj]],scale_res)
        kkk=kkk+1
      }
    }
    sum_full_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
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
    sum_red_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
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
    sum_red_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
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
    sum_red_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
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
    sum_red_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_red_list[[ii]]$data),along = 1)
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
}
RoFanova_twoway_perm<-function(X,label_1,label_2=NA,B=100,...){

  if(is.na(label_2[1])){
    cat("one-way")
    return(RoFanova_perm(X,label_1,B=B,...))
  }
  else{
    n<-length(label_1)
    mod_obs<-RoFanova_twoway_sur(X,label_1,label_2,...)
    Tr_obs<-mod_obs$Tr_tot


    per_fun<-function(kkk){
      perm_comb<-sample(1:n,replace = F)
      X_per<-ex_fdata(X,perm_comb)
      mod<-RoFanova_twoway_sur(X_per,label_1,label_2,...)

      out<-list(mod$Tr_tot)
      return(out)
    }
    per_list_out<-pbmclapply(1:B,per_fun,mc.cores = detectCores()-1)
    per_list<-lapply(per_list_out,"[[",1)

    pval_vec<-Tr_obs
    Tr_perm<-list()
    for (ii in 1:nrow(pval_vec)) {
      Tr_perm[[ii]]<-sapply(1:B,function(ll)per_list[[ll]][ii])
      pval_vec[ii]<-sum(Tr_perm[[ii]]>=Tr_obs[ii])/B
    }






  }

  out<-list(pval_vec=pval_vec,
            Tr_obs=Tr_obs,
            Tr_perm=Tr_perm)
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
      sum_full_list[[kkk]]<-(group_mean_ij[[ii]][[jj]]-global_mean)^2
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
  sum_full_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
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
  sum_full_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
  SSmod<-acc_fdata_sur(fdata(sum_full_data,argvals =grid ))
  T_numden_mod<-int_sur(SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)))/int_sur(SSE/(k_1*k_2*(n_ij-1)))
  T_full_mod<-int_sur(frac_fdata_sur(SSmod/((k_1-1)+(k_2-1)+(k_1-1)*(k_2-1)),SSE/(k_1*k_2*(n_ij-1))))
  T_mod<-c(numden=T_numden_mod,full=T_full_mod)

  # Int ---------------------------------------------------------------------


  kkk=1
  sum_full_list<-list()
  for (ii in 1:k_1) {
    for (jj in 1:k_2) {
      ele<-sum_fdata_sur(diff_fdata(group_mean_ij[[ii]][[jj]],sum_fdata_sur(group_mean_1[[ii]],group_mean_2[[jj]])),global_mean)
      sum_full_list[[kkk]]<-n_ij*(ele)^2
      kkk=kkk+1
    }
  }
  sum_full_data<-abind(lapply(1:(k_1*k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
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
  sum_full_data<-abind(lapply(1:(k_1), function(ii)sum_full_list[[ii]]$data),along = 1)
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
  sum_full_data<-abind(lapply(1:(k_2), function(ii)sum_full_list[[ii]]$data),along = 1)
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
# Varie -------------------------------------------------------------------


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



func.mean_sur<-function (fdataobj)
{
  n<-dim(fdataobj$data)[1]
  pro<-array(NA,dim = c(1,dim(fdataobj$data)[2],dim(fdataobj$data)[3]))
  pro[1,,] <-apply(fdataobj$data , c(2,3), mean)
  fdataobj$data<-pro
  fdataobj$names$main <- "mean"
  fdataobj
}
func.var_sur<-function (fdataobj)
{

  n<-dim(fdataobj$data)[1]
  pro<-array(NA,dim = c(1,dim(fdataobj$data)[2],dim(fdataobj$data)[3]))
  pro[1,,] <-apply(fdataobj$data , c(2,3), var)
  fdataobj$data<-pro
  fdataobj$names$main <- "var"
  fdataobj
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
frac_fdata_sur<-function(x1,x2){
  n<-dim(x1$data)[1]
  out<-x1
  for(ii in 1:n)out$data[ii,,]=x1$data[ii,,]/(x2$data[ii,,]+1e-15)
  return(out)
}
acc_fdata_sur<-function(x){
  out<-x
  pro<-array(NA,dim = c(1,dim(x$data)[2],dim(x$data)[3]))
  pro[1,,] <-apply(x$data , c(2,3), sum)
  out$data<-pro
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
int_sur<-function(x){

  grid_s<-x$argvals[[1]]
  grid_t<-x$argvals[[2]]
  delta_s=1/(length(grid_s)-1)
  delta_t=1/(length(grid_t)-1)
  delta_s*delta_t*sum(x$data[1,,])

}
tra_data<-function(X_fdata){


  out<-X_fdata

  sss<-powerTransform(norm_fdata_c_sur(out))
  out$data= bcPower( out$data+abs(min(out$data))+0.00001,0.5)
  return(out)

}
alignment_fdata_sur<-function(X_fdata,ref){

  if(is.numeric(ref)){
    X_1<-ex_fdata(X_fdata,ref)
    n=dim(X_fdata$data)[1]
    n_minus<-(1:n)[-ref]
  }
  else{
    X_1<-ref
    n=dim(X_fdata$data)[1]
    n_minus<-(1:n)
  }
  X_fdata_new<-X_fdata

  for (kkk in 1:length(n_minus)) {

    X_2<-ex_fdata(X_fdata,n_minus[kkk])
    grid<-expand.grid(-4:4,-4:4)
    norm_vec<-numeric()
    X_2_new<-list()
    for (ii in 1:length(grid[,1])) {
      data_new<-X_2$data
      dim_vec<-dim(data_new)
      delta_ind<-as.numeric(grid[ii,])
      if(delta_ind[1]>0){
        data_new[1,(1+delta_ind[1]):dim_vec[2],]<-X_2$data[1,1:(dim_vec[2]-delta_ind[1]),]
      }
      else if (delta_ind[1]<0){
        delta_ind[1]<-abs(delta_ind[1])
        data_new[1,1:(dim_vec[2]-delta_ind[1]),]<-X_2$data[1,(1+delta_ind[1]):dim_vec[2],]
      }

      if(delta_ind[2]>0){
        data_new[1,,(1+delta_ind[2]):dim_vec[3]]<-data_new[1,,1:(dim_vec[3]-delta_ind[2])]
      }
      else if (delta_ind[2]<0){
        delta_ind[2]<-abs(delta_ind[2])
        data_new[1,,1:(dim_vec[3]-delta_ind[2])]<-data_new[1,,(1+delta_ind[2]):dim_vec[3]]

      }
      X_2_new[[ii]]=X_2
      X_2_new[[ii]]$data=data_new
      norm_vec[ii]<-norm_fdata_c_sur(diff_fdata_sur(X_2_new[[ii]],X_1))
    }

    min_ind<-which(norm_vec==min(norm_vec))
    X_fdata_new$data[n_minus[kkk],,]<-X_2_new[[min_ind]]$data
  }

  return(X_fdata_new)

}
alignment_fdata_sur_cv<-function(X_fdata){
  n=dim(X_fdata$data)[1]
  median_diff<-numeric()
  par_fun<-function(ii) {
    X_fdata_new<-alignment_fdata_sur(X_fdata,ii)
    mean(parwise_diff_sur(X_fdata_new ))
  }
  mean_diff<-mclapply(1:n,par_fun,mc.cores=detectCores()-1)

  ind_min<-max(which(unlist(mean_diff)==min(unlist(mean_diff))))
  X_fdata_new<-alignment_fdata_sur(X_fdata,ind_min)
  return(X_fdata_new)
}

der_sur<-function(X_fdata,nderiv=1,dim="s"){
  n=dim(X_fdata$data)[1]
  ss = X_fdata$argvals[[1]]
  tt = X_fdata$argvals[[2]]

  der_fdata<-X_fdata


  for (ii in 1:n) {
    cat(ii)
    data<-X_fdata$data[ii,,]
    ns <- nrow(data)
    nt <- ncol(data)
    if (dim == "s") {
      res = matrix(NA, nrow = ns, ncol = nt)
      for (i in 1:nt) {
        a = diff(data[,i], differences = nderiv)/(ss[2:ns] -
                                                    ss[1:(ns - 1)])^nderiv
        ab = matrix(NA, ncol = ns, nrow = 2)
        ab[1, 2:ns] = a
        ab[2, 1:(ns - 1)] = a
        res[,i ] = colMeans(ab, na.rm = TRUE)
      }
      der_fdata$data[ii,,]<-res
    }
    if (dim == "t") {
      res = matrix(NA, nrow = ns, ncol = nt)
      for (i in 1:ns) {
        a = diff(data[i,], differences = nderiv)/(tt[2:nt] -
                                                    tt[1:(nt - 1)])^nderiv
        ab = matrix(NA, ncol = nt, nrow = 2)
        ab[1, 2:nt] = a
        ab[2, 1:(nt - 1)] = a
        res[i, ] = colMeans(ab, na.rm = TRUE)
      }
      der_fdata$data[ii,,]<-res
    }
  }
  return(der_fdata)
}

twosample_test_sur<-function(X_1,X_2,B=100,family="bisquare",eff = 0.95, maxit = 100, tol = 1e-30,...){


  n_1<-dim(X_1$data)[1]
  n_2<-dim(X_2$data)[1]
  n=n_1+n_2
  X=fdata(abind(X_1$data,X_2$data,along = 1),argvals = X_1$argvals)
  scale_1<-scale_fun_pw_sur(X_1,eff = eff,tol=tol, maxit = maxit)#4)
  scale_2<-scale_fun_pw_sur(X_2,eff = eff,tol=tol, maxit = maxit )#4)
  # scale_1<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
  # scale_2<-fdata(array(1,dim = c(1,dim(X$data)[2],dim(X$data)[3])),argvals =X$argvals)
  if(n_1==2|n_2==2){
    group_mean_1<-FlocScaleM_sur(X_1,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
    group_mean_2<-FlocScaleM_sur(X_2,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
  }
  else{
    group_mean_1<-FlocScaleM_sur(X_1,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_1)$mu
    group_mean_2<-FlocScaleM_sur(X_2,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_2)$mu
  }
  T_obs<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))
  # # if(is.na(scale_glob))scale_glob<-scale_fun_pw_sur(X,eff = eff,tol=tol, maxit = maxit,... )#5
  if(n_1==2|n_2==2){
    per_fun<-function(kkk){
      perm_comb<-sample(1:n,replace = F)
      X_1per<-ex_fdata(X,perm_comb[1:n_1])
      X_2per<-ex_fdata(X,perm_comb[(n_1+1):n])
      group_mean_1<-FlocScaleM_sur(X_1per,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
      group_mean_2<-FlocScaleM_sur(X_2per,psi = family, eff = eff, maxit = maxit, tol =tol)$mu
      T_i<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))

      return(T_i)
    }
  }
  else{
    per_fun<-function(kkk){
      perm_comb<-sample(1:n,replace = F)
      X_1per<-ex_fdata(X,perm_comb[1:n_1])
      X_2per<-ex_fdata(X,perm_comb[(n_1+1):n])
      scale_1<-scale_fun_pw_sur(X_1per,eff = eff,tol=tol, maxit = maxit)#4)
      scale_2<-scale_fun_pw_sur(X_2per,eff = eff,tol=tol, maxit = maxit )#4)
      group_mean_1<-FlocScaleM_sur(X_1per,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_1)$mu
      group_mean_2<-FlocScaleM_sur(X_2per,psi = family, eff = eff, maxit = maxit, tol =tol,sig0_g = scale_2)$mu
      T_i<-as.numeric( int_sur(abs(diff_fdata_sur(group_mean_1,group_mean_2))))

      return(T_i)
    }
  }
  T_per<-unlist(mclapply(1:B,per_fun,mc.cores = detectCores()-1))

  pvalue<-sum(T_per>=(T_obs+1e-10))/B
  out<-list(pvalue=pvalue,
            T_obs=T_obs,
            T_per=T_per)
  return(out)
}
pairwise_comparisons<-function(X_fdata,label_1,label_2,...) {

  k_1=length(unique(label_1))
  k_2=length(unique(label_2))

  pairwise<-cbind(expand.grid(1:k_2,1:k_1)[2],expand.grid(1:k_2,1:k_1)[1])
  mat<-matrix(0,dim(pairwise)[1],dim(pairwise)[1])
  for (ii in 1:dim(pairwise)[1]) {
    for (jj in 1:ii) {
      cat(paste(as.numeric(pairwise[ii,])));cat(" - ");cat(as.numeric(pairwise[jj,]));cat("\n")
      iii=pairwise[ii,1]; jjj=pairwise[ii,2]
      kkk=pairwise[jj,1]; lll=pairwise[jj,2]
      mat[ii,jj]<-twosample_test_sur(X_1=ex_fdata(X_fdata,label_1==iii&label_2==jjj),X_2=ex_fdata(X_fdata,label_1==kkk&label_2==lll),...)$pvalue
    }
  }
  colnames(mat)<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))
  rownames(mat)<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))


  return(mat)


}
get_matrix_pairwise<-function(mod,label_1,label_2,fac=1){

  k_1=length(unique(label_1))
  k_2=length(unique(label_2))
  int<-levels(interaction(1:k_1,1:k_2,lex.order = TRUE))
  matrix_list<-list()
  ind<-ifelse(fac==1,k_1,k_2)
  for (ii in 1:ind) {
    if(fac==1)     ind_ii<-which(substr(int, 1, 1)==ii)
    else  ind_ii<-which(substr(int, 3, 3)==ii)
    matrix_list[[ii]]<-mod[ind_ii,ind_ii]

  }
  return(matrix_list)
}
