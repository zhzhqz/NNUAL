NNUAL<-function(ZZnchu,TTnchu,YYnchu,ntrain,K_n,d)
{
n_zong<-length(ZZnchu)#####sample size+test set sample size
n=ntrain########sample size
q_n<-23################the number of spline basis functions for eigenfunctions
gq_n=5#################the number of knots for link
pesai1q_n=2################the number of knots for pedai1
pesai2q_n=2################the number of knots for pedai1

tuning_lambda<-0.6##########penalty parameter for SCAD 
tuning_a<-3.7####tunning parameter for SCAD
v<-0.7########tunning parameter for w
################function F
yhat<-matrix(NA,200,3)

################spline basis for eigenfunctions
kontx=(seq(0,1,length=(q_n-2)))[2:(q_n-3)]
splinex<-function(t) bs(t,df=NULL,kontx,degree=3,intercept=T,Boundary.knots=c(0,1))
work<-matrix(0,q_n,q_n)
for (i in 1:10000)
{work=work+splinex(i/10000)[1,]%o%splinex(i/10000)[1,]
}
working=work/10000
ortho_matrix<-solve(sqrtm(working)$B)
splinex_pca<- function(t) splinex(t)%*%ortho_matrix

ZZ<-list()
TT<-list()
YY<-list()
ni_zong<-rep(NA,n_zong)
for (i in 1:n_zong) {
  ZZ[[i]]<-ZZnchu[[i]] ####################functional covariate observation
  TT[[i]]<-TTnchu[[i]]###################obvserved time
  YY[[i]]<-YYnchu[[i]]###################response
  ni_zong[i]<-length(TT[[i]])
}
ni<-ni_zong[1:n]
yy<-matrix(NA,n,1)
for (i in 1:n)
{
  yy[i,]<-YY[[i]]
}
y_test<-matrix(NA,(n_zong-n),1)
for (i in (n+1):n_zong)
{
  y_test[i-n,]<-YY[[i]]
}

alpha<-matrix(0,q_n,1)
Gamma<-matrix(0,K_n,q_n)
UU<-matrix(0,n,K_n)
BB<-matrix(NA,d,K_n)
BBxia<-matrix(NA,d,(K_n-2))
Delta<-rep(NA,(gq_n+4))
theta1<-rep(NA,(pesai1q_n+4))
theta2<-rep(NA,(pesai2q_n+4))

 #######################initialization
samp=seq(0,1,length=100)
samp_grid<-matrix(samp,n_zong,100,byrow=T)
desig<-splinex_pca(samp)
project<-solve(t(desig)%*%desig)%*%t(desig)
splines <- create.bspline.basis(rangeval=c(0, 1), nbasis = 20,norder=4)
zzgrid<-matrix(NA,n_zong,100)

for(i in 1:n_zong)
{
  zzgrid[i,]<-eval.fd(samp,Data2fd(unlist(ZZ[[i]]), argvals=unlist(TT[[i]]), basisobj=splines))
}
fdata<-Data2fd(t(zzgrid), argvals=t(samp_grid), basisobj=splines)
fpca<-pca.fd(fdata,nharm=K_n,harmfdPar=fdPar(fdata),centerfns = TRUE)
theta<-project%*%(eval.fd(samp,fpca$meanfd))
Gamma<-t(project%*%eval.fd(samp,fpca$harmonics)%*%diag((fpca$values[1:K_n])^(0.5)))
UUtt<-matrix(NA,n_zong,K_n)
for (i in 1:n_zong)
{
  BX<-splinex_pca(TT[[i]])
  a<-t(Gamma%*%t(BX))
  b<-c(ZZ[[i]]-BX%*%alpha)
  UUtt[i,]<-solve((t(a)%*%a),t(a)%*%b)+rnorm(K_n,0,0.05)
}

UU<-UUtt[1:n,]
UU_test<-UUtt[(n+1):n_zong,]
XX<-UU

dat<-data.frame(y=yy,x1=XX[,1],x2=XX[,2],x3=XX[,3],x4=XX[,4],x5=XX[,5],x6=XX[,6],x7=XX[,7],x8=XX[,8],x9=XX[,9],x10=XX[,10],x11=XX[,11],x12=XX[,12],x13=XX[,13],x14=XX[,14],x15=XX[,15],x16=XX[,16],x17=XX[,17],x18=XX[,18],x19=XX[,19],x20=XX[,20])
bc <- glm(y ~ x1+ x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15+ x16+ x17+ x18+ x19+ x20,data=dat)
ymm<-bc$fitted.values
beta<-matrix(NA,d,K_n)
beta[1,]<-(bc$coefficients)[2:(K_n+1)]
beta[2,]<-(bc$coefficients)[2:(K_n+1)]
beta<-beta/sqrt(sum(beta^2))
BB<-beta

zhibiao_init<-UU%*%t(BB)
#####################spline basis for pesai1
bg<-c(min(zhibiao_init[,1]),max(zhibiao_init[,1]))
pesai1_kontx=seq(bg[1],bg[2],length=(2+pesai1q_n))[2:(1+pesai1q_n)]
Boundary=bg+c(-1,1)
splinex_pesai1<-function(t) bs(t,df=NULL,pesai1_kontx,degree=3,intercept=T,Boundary.knots=Boundary)

samp=seq(bg[1],bg[2],length=300)
desig<-splinex_pesai1(zhibiao_init[,1])
project<-ginv(t(desig)%*%desig)%*%t(desig)
theta1<-project%*%ymm
pesai1_hat<-function(x) splinex_pesai1(x)%*%theta1

#####################spline basis for pesai2
bg<-c(min(zhibiao_init[,2]),max(zhibiao_init[,2]))
Boundary=bg+c(-1,1)
pesai2_kontx=seq(bg[1],bg[2],length=(2+pesai2q_n))[2:(1+pesai2q_n)]
Boundaryg<-bg+c(-1,1)
splinex_pesai2<-function(t)   bs(t,df=NULL,pesai1_kontx,degree=3,intercept=T,Boundary.knots=Boundary) #bs#(t,df=NULL,pesai2_kontx,degree=3,intercept=T,Boundary.knots=Boundary)

samp=seq(bg[1],bg[2],length=300)
desig<-splinex_pesai2(zhibiao_init[,2])
project<-ginv(t(desig)%*%desig)%*%t(desig)
theta2<-project%*%(ymm-pesai1_hat(zhibiao_init[,1]))
pesai2_hat<-function(x) splinex_pesai2(x)%*%theta2

#####################spline basis for link function
g_splineX=pesai1_hat(zhibiao_init[,1])+pesai2_hat(zhibiao_init[,2])
bg<-c(min(g_splineX),max(g_splineX))
g_kontx=seq(bg[1],bg[2],length=(2+gq_n))[2:(1+gq_n)]
Boundaryg<-bg+c(-1,1)
splinex_Delta<-function(t) bs(t,df=NULL,g_kontx,degree=3,intercept=T,Boundary.knots=Boundaryg)

samp=seq(bg[1],bg[2],length=300)
desig<-splinex_Delta(g_splineX)
project<-ginv(t(desig)%*%desig)%*%t(desig)
Delta<-project%*%ymm
ghat<-function(x) splinex_Delta(x)%*%Delta

SStep<-0
alpha_1<-alpha
Gamma_1<-Gamma
UU_1<-UU
BB_1<-BB
Delta_1<-Delta
theta1_1<-theta1
theta2_1<-theta2
dif=10
dd<-seq(0.000001,1,length=100)
mubiao=10

#######################################begin iteration
while((SStep<20)&(dif>10^(-3))){
  SStep<-SStep+1
  mubiao1=mubiao
  ############################update for meanfunction, eigenfucntions
  alpha_matrix<-matrix(0,q_n,q_n)
  alpha_vector<-matrix(0,q_n,1)
  
  Gamma_vector<-array(0,dim=c(q_n,1,K_n))
  Gamma_matrix<-array(0,dim=c(q_n,q_n,K_n))
  for (i in 1:n)
  {
    BX<-splinex_pca(TT[[i]])
    VV1<-t(BX)%*%ZZ[[i]]
    VV2<-t(BX)%*%BX
    alpha_vector<-alpha_vector+VV1-VV2%*%t(UU_1[i,]%*%Gamma_1)
    alpha_matrix<-alpha_matrix+VV2
    for (k in 1:K_n)
    {
      Gamma_vector[,,k]<-Gamma_vector[,,k]+(VV1-VV2%*%alpha_1-VV2%*%t(UU_1[i,]%*%Gamma_1-UU_1[i,k]%*%Gamma_1[k,]))*UU_1[i,k]
      Gamma_matrix[,,k]<-Gamma_matrix[,,k]+VV2*UU_1[i,k]*UU_1[i,k]
    }
  }
  alpha<-solve(alpha_matrix)%*%alpha_vector
  for (k in 1:K_n)
  {
    Gamma[k,]<-t(solve(Gamma_matrix[,,k])%*%Gamma_vector[,,k])
  }
  alpha_1<-alpha
  Gamma_1<-Gamma
  
  ############################define some derivatives
  pesai1dao<-rep(NA,n)
  pesai2dao<-rep(NA,n)
  daoshu<-matrix(NA,n,2)
  youhuaBB<-matrix(NA,d,K_n)
  zhibiao<-matrix(NA,n,d)
  varY_inverse<-rep(NA,n)
  
  zhibiao_var<-UU_1%*%t(BB_1)
  g_var=pesai1_hat(zhibiao_var[,1])+pesai2_hat(zhibiao_var[,2])
  varY_inverse<-1/(ghat(g_var)*(1-ghat(g_var)))
  gdao<-function(x) fderiv(ghat,x)
  pesai1dao<-function(x) fderiv(pesai1_hat,x)
  pesai2dao<-function(x) fderiv(pesai2_hat,x)
  
  daoshu<-function(x)#
  {
    youhuaBB<-t(matrix(x,K_n,d))
    zhibiao<-UU_1%*%t(youhuaBB)
    t1=zhibiao[,1]
    t2=zhibiao[,2]
    t=pesai1_hat(zhibiao[,1])+pesai2_hat(zhibiao[,2])
    varY_inverse<-1/(ghat(t)*(1-ghat(t)))
    daoshu<-matrix(varY_inverse*(yy-ghat(t))*gdao(t),n,2,byrow=F)*cbind(pesai1dao(t1),pesai2dao(t2))
    return(daoshu)
  }
  
  
  ############################update for score  
  daoshuU<-matrix(NA,n,2)
  BBchu<-matrix(t(BB_1),(K_n*d),1)
  daoshuU<-daoshu(BBchu)
  w<-(min(ni))^(-v)
  for (i in 1:n)
  {
    BX<-splinex_pca(TT[[i]])
    VV1<-t(BX)%*%ZZ[[i]]
    VV2<-t(BX)%*%BX
    UU[i,]<-solve(Gamma_1%*%VV2%*%t(Gamma_1))%*%(t(BB_1)%*%daoshuU[i,]/(2*w)+Gamma_1%*%(VV1-VV2%*%alpha_1))
  }

  ############################update for B
  penalty_dao<-function(t)
  {
    ifelse (t<=tuning_lambda,tuning_lambda,max(0,tuning_a*tuning_lambda-t)/(tuning_a-1))
  }
  BB_colnorm<-sqrt(colSums(BB^2))
  GB<-diag(penalty_dao(BB_colnorm)/BB_colnorm)
  miniBB<-function(x)
  {
    youhuaBB<-t(matrix(x,K_n,d))
    zhibiao<-UU_1%*%t(youhuaBB)
    t1=zhibiao[,1]
    t2=zhibiao[,2]
    t=pesai1_hat(zhibiao[,1])+pesai2_hat(zhibiao[,2])
    varY_inverse<-1/(ghat(t)*(1-ghat(t)))
    return(sum((t(daoshu(youhuaBB))%*%UU_1/n-youhuaBB%*%diag(penalty_dao(BB_colnorm)/sqrt(colSums(youhuaBB^2))))^2))
  }
  nonzero_index<-which(abs(BB[1:d,1:K_n])>0.005) 
  if(SStep>1)
  {   zero_index<-which(abs(BB[1:d,1:K_n])<0.02)
  BB[zero_index]<-10^(-5)
  nonzero_index<-which(abs(BB[1:d,1:K_n])>0.02)
  }
  BBjchu<-matrix(BB,(K_n*d),1)
  yijie<-function(x) jacobian(miniBB,x)
  newton<-ginv(100*(jacobian(yijie,BBjchu))[nonzero_index,nonzero_index])%*%(matrix(yijie(BBjchu),(K_n*d),1))[nonzero_index]
  BB[nonzero_index]<-BB[nonzero_index]-newton
  if(SStep>1)
    BB[zero_index]<-10^(-5)
  zero_index<-which(abs(BB[1:d,1:K_n])<0.02)
  BB[zero_index]<-10^(-5)
  
  for (j in 1:d)
  {
    BB[j,]<-BB[j,]*sum(abs(beta[j,]))/sum(abs(BB[j,]))
  } 
  
  ############################update for link function
  zhibiao_var<-UU_1%*%t(BB)
  g_var=pesai1_hat(zhibiao_var[,1])+pesai2_hat(zhibiao_var[,2])
  varY_inverse<-1/(ghat(g_var)*(1-ghat(g_var)))
    Delta<-ginv(t(splinex_Delta(g_var))%*%diag(c(varY_inverse))%*%splinex_Delta(g_var))%*%t(splinex_Delta(g_var))%*%diag(c(varY_inverse))%*%ymm
  ghat<-function(x) splinex_Delta(x)%*%Delta
  
  zhibiao_var<-UU_1%*%t(BB)
  g_var=pesai1_hat(zhibiao_var[,1])+pesai2_hat(zhibiao_var[,2])
  varY_inverse<-1/(ghat(g_var)*(1-ghat(g_var)))
  
  ############################update for pesai1
  miniPesai1<-function(x) 
  {
    zhibiao<-UU_1%*%t(BB)
    t1=zhibiao[,1]
    pesai1_hatyouhua<- splinex_pesai1(zhibiao[,1])%*%x
    t2=zhibiao[,2]
    t= pesai1_hatyouhua+pesai2_hat(zhibiao[,2])
    daoshu<-t(splinex_pesai1(t1))%*%(varY_inverse*(yy-ghat(t))*gdao(t))
    return(sum((daoshu/n)^2))
  }
 
  ei<-diag(1,length(theta1))
  theta1_chu<-matrix(theta1,length(theta1),1)
  theta1_yijie<-function(x) jacobian(miniPesai1,x)
  erjie1<-jacobian(theta1_yijie,theta1_chu)
  newton<-ginv(erjie1)%*%(matrix(theta1_yijie(theta1_chu),length(theta1),1))
  theta1<-theta1-newton
  
  ############################update for pesai2
  miniPesai2<-function(x) 
  {
    zhibiao<-UU_1%*%t(BB)
    t1=zhibiao[,1]
    t2=zhibiao[,2]
    pesai2_hatyouhua<- splinex_pesai2(zhibiao[,2])%*%x
    t=pesai1_hat(zhibiao[,1])+pesai2_hatyouhua
    daoshu<-t(splinex_pesai2(t2))%*%(varY_inverse*(yy-ghat(t))*gdao(t))
    return(sum((daoshu/n)^2))
  }
  ei<-diag(1,length(theta2))
  theta2_chu<-matrix(theta2,length(theta2),1)
  
  theta2_yijie<-function(x) jacobian(miniPesai2,x)
  erjie2<-jacobian(theta2_yijie,theta2_chu)
  newton<-ginv(erjie2)%*%(matrix(theta2_yijie(theta2_chu),length(theta2),1))
  theta2<-theta2-newton
  
  difff<-c(mean((alpha_1-alpha)^2)/mean((alpha_1)^2) ,mean((Gamma_1-Gamma)^2)/mean((Gamma_1)^2),mean((UU_1-UU)^2)/mean((UU_1)^2) ,mean((BB_1-BB)^2)/mean((BB_1)^2)
           ,mean((theta1_1-theta1)^2)/mean((theta1_1)^2),mean((theta2_1-theta2)^2)/mean((theta2_1)^2))
  dif<-max(difff[4])
  alpha_1<-alpha
  Gamma_1<-Gamma
  BB_1<-BB    
  UU_1<-UU
  Delta_1<-Delta
  theta1_1<-theta1
  theta2_1<-theta2
  ghat<-function(x) splinex_Delta(x)%*%Delta_1
  pesai1_hat<-function(x) splinex_pesai1(x)%*%theta1_1
  pesai2_hat<-function(x) splinex_pesai2(x)%*%theta2_1
  samp=seq(-3,5,length=300)
}
#########################the end of iteration


########################################prediction for test set
zhibiao_cv<-UU_test%*%t(BB)
aa<-ghat(pesai1_hat(zhibiao_cv[,1])+pesai2_hat(zhibiao_cv[,2]))
aa[which(aa>1)]<-1
aa[which(aa<0)]<-0
########################prediction error
#Updated_Est<-list(scores=UU,B=BB,ypre=mean((aa-y_test)^2))
splinesg <- create.bspline.basis(rangeval=c(min(g_splineX), max(g_splineX)), nbasis = 10,norder=4)
gfini<-Data2fd((splinex_Delta(g_splineX)%*%Delta_1), argvals=g_splineX, basisobj=splinesg)
gf<-fd(gfini$coef,splinesg)

splinespsai1 <- create.bspline.basis(rangeval=c(min(zhibiao_init[,1]),max(zhibiao_init[,1])), nbasis = 10,norder=4)
psaifini1<-Data2fd(cbind(splinex_pesai1(zhibiao_init[,1])%*%theta1_1,splinex_pesai2(zhibiao_init[,2])%*%theta2_1),argvals=zhibiao_init, basisobj=splinespsai1)
psaif1<-fd(psaifini1$coef,splinespsai1)

Updated_Est<-list(scores=UU,B=BB,ypre=mean((aa-y_test)^2),gf=gf,psaif1=psaif1)
return(Updated_Est)

}
