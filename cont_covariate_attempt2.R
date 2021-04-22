rm(list=ls())
library(Matrix)
library(magic)
library(psych)
library(ggplot2)
library(reshape2)
library(MASS)
library(pracma)
library(latex2exp)
library(ks)
###################################################################

logit <- function(x) {if ((x==0)|(x==1)) return(0) else return(log(x/(1-x))) }
Afunc <- function(x) {if (x==0) return(-.125) else return(-tanh(x/2)/(4*x))}
Cfunc <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
#####FUNCTION DEFINITION#####
#############################
#############################
ELBO_calculator=function(y,X_mat,S_sq,mu,alpha, sigmasq, sigmabeta_sq, true_pi, W,n,p){
  mu=matrix(mu,p,1)
  alpha=matrix(alpha,p,1)
  s=matrix(S_sq,p,1)
  mu_alpha= matrix(mu*alpha,p,1)
  W=matrix(W,n,1)
  t1=-sum(W*(y- X_mat%*%mu_alpha)^2)/(2*sigmasq)
  t2=- sum(W*((X_mat)^2%*%(alpha*(mu^2 + s) - alpha^2*mu^2)))/(2*sigmasq)
  t3= sum(alpha*((1+ log(s))))/2
  t4= -sum(alpha*log((alpha+0.000001)/true_pi) + (1-alpha)*log((1-alpha+0.000001)/(1-true_pi)))
  t5= - sum(alpha*((mu^2 + s)/(2*sigmasq*sigmabeta_sq) + log(sigmasq*sigmabeta_sq)/2))
  t6=sum(0.5*log(1/(2*pi*sigmasq)))
  t=t1+t2+t3+t4+t5+t6
  return(t)
}
#############################
#############################
#############################
######### Data generation ############
n <- 180
p <- 4
MAXITER=2
STR=1
in_pr_13=matrix(0,MAXITER,n)
in_pr_12=in_pr_13
Var_cont=function(z){
  pr=matrix(0,p+1,p+1)
  diag(pr)=2
  #  pr[1,2]=STR*((z>0) && (z< .33)) + (STR - STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (0)*((z>0.66) && (z<1))
  #  pr[1,3]=0*((z>0) && (z< .33)) + (STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (STR)*((z>0.66) && (z<1))
  pr[2,3]=STR
  pr[1,2]=STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56))*((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3]=0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56))*((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))
  
  pr[2,1]=pr[1,2]
  pr[3,1]=pr[1,3]
  pr[3,2]=pr[2,3]
  
  
  Var=solve(pr)
  return(Var)
}
sensitivity_20=matrix(0,MAXITER,1)
specificity_20=sensitivity_20
sensitivity_90=sensitivity_20
specificity_90=sensitivity_20
sensitivity_160=sensitivity_20
specificity_160=sensitivity_20
for(overiter in 1:MAXITER){
  # Z=c(runif(n/3,-1,-0.33),runif(n/3, -0.23,0.33),runif(n/3,0.43,1))
  Z=c(seq(-0.99,-0.331,(-.331+.99)/59),seq(-0.229,0.329,(.329+.229)/59),seq(0.431,.99,(.99-.431)/59))
  # Z=seq(0.01,.99,.98/(n-1))
  Z=matrix(Z,n,1)
  X=matrix(0,n,p+1)
  for(i in 1:n){
    X[i,]=  mvrnorm(1,rep(0,p+1),Var_cont(Z[i]))
  }
  
  ######### Generating the covariates ##########
  
  
  beta=matrix(0,n,p)
  resp_index=1;# The index we consider as response
  mylist <- rep(list(beta),p+1)
  data_mat=X
  
  
  
  
  
  ###############################################
  for(resp_index in 1:(p+1)){
    
    
    ##############
    y=data_mat[,resp_index];
    
    X_mat=data_mat[,-resp_index]
    X_vec <- matrix(0,n*p,1)
    #X<- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
    
    X<- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
    
    for(i in 1:n){
      for(j in 1:p){
        k=p*(i-1)+1
        X[i,k+j-1]=X_mat[i,j]
        X_vec[k+j-1]=X[i,k+j-1]
      }
    }
    ELBO_LBit=rep(0,10000)
    Big_diag_mat <- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
    for(i in 1:n){
      k=p*(i-1)
      for(j in 1:p){
        Big_diag_mat[i,k+j]=1
      }
    }
    A_xi=rep(1,n)
    #X <- matrix(rnorm(n*p,0,1), ncol=p)
    q=matrix(2,n,1)
    
    
    # for(i in 1:n){
    #   q[i]=floor()
    #   for(j in 1:q[i]){
    #     beta[(i-1)*p+j] <- 4 # runif(1,3,3.8) ###c(rep(3.5,q),rep(0,p-q))
    #   }
    # }
    #beta ### Done
    sigmasq=1
    E <- rnorm(n,0,sigmasq)
    #snr=sd(X%*%as.vector(beta0))/sd(sigmasq*E)
    
    #y
    
    XtX=t(X)%*%X
    
    DXtX=diag(XtX)
    DXtX_rep=rep(DXtX,p); dim(as.matrix(DXtX_rep))
    DXtX_mat=matrix(DXtX_rep,n*p,p,byrow=FALSE)
    Diff_mat=XtX-diag(DXtX)
    
    
    tt1=kde(Z[,1],eval.points=Z[,1])
    h1=tt1$h
    # tt2=kde(Z[,2],eval.points=Z[,2])
    # h2=tt2$h
    
    D=matrix(1,n,n)
    for(i in 1:n){
      for(j in 1:n){
        # D[j,i]= dnorm(sqrt(norm((Z[i,1]-Z[j,1])/h1,"2")^2 + norm((Z[i,2]-Z[j,2])/h2,"2")^2),0,1)
        D[j,i]= dnorm(norm(Z[i,]-Z[j,],"2"),0,0.56)
      }
    }
    for(i in 1:n){
      D[,i]=n*(D[,i]/sum(D[,i]))
      #    D[,i]=1
    }
    #true_lambda=rep(-log((p-q)/q)/p,p)
    true_lambda=0.5*rep(1,p)
    L0=0.5
    lambda_mean=true_lambda##rep(0,p) ###rnorm(p,0,4)
    lambda_var=.001*diag(p)
    mu0_lambda<-L0*rep(1,p)## rep(0,p)
    Sigma0_lambda=lambda_var###diag(p)
    alpha= rep(0.2,n*p)
    sigmabeta_sq=3
    mu=rep(0,p)
    true_pi=0.5
    
    
    ###########
    
    
    
    y_long_vec=as.vector(t(y%*%matrix(1,1,p)))
    Xty=t(X)%*%as.vector(y)
    beta_mat=matrix(beta,n,p,byrow=TRUE)
    mu_mat=beta_mat
    
    
    D_long=matrix(0,n*p,n)
    for( i in 1:n){
      D_long[,i]=matrix(t(D[,i]%*%matrix(1,1,p)),n*p,1)
    }
    
    
    S_sq=matrix(sigmasq*(DXtX + 1/sigmabeta_sq)^(-1),n,p)
    Sigma_xi=diag(p)
    S0=solve(lambda_var)
    iter=1
    ###############################################
    
    
    
    
    ind_vec=seq(0,(n-1)*p,by=p)
    Ind_mat=matrix(0,n,p)
    for(j in 1:p){
      Ind_mat[,j]=ind_vec+j
    }
    Big_ind=matrix(0,n*p,p)
    Big_ind_1=matrix(0,n*p,p)
    for(j in 1:p){
      Big_ind[Ind_mat[,j],j]=1
      Big_ind_1[Ind_mat[,j],j]=0
    }
    
    DXtX_Big_ind=DXtX_mat*Big_ind
    ###################################################
    
    cov_vsvb= function(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq,sigmabeta_sq,lambda_mean,lambda_var,L0,true_pi){
      
      thres=1e-7
      tol=1e-9
      mu0_lambda<-L0*rep(1,p)
      msg <- function(s, ...)
      {
        time <- format(Sys.time(), "%X")
        cat(sprintf("%s %s\n", time, s))
      }
      
      
      change_alpha <- rep(0.001,n*p) #alpha_new - alpha_int
      
      #change_ELBO_LB = ELBO_LB - ELBO_LB_int
      #  sigmabeta_sq=5*sigmasq
      max_iter <- 100
      iter=1
      Mu_vec=matrix(rep(mu,n),n*p,1)
      while(sqrt(sum(change_alpha^2))>tol & iter<max_iter){
        #    while(iter<max_iter){
        alpha_int=alpha
        # alpha_int
        # ELBO_LB0=ELBO_LB
        #  msg(sprintf("iteration %d; ",iter))
        
        alpha_mat=matrix(alpha,n,p,byrow=TRUE)
        
        alpha_vec=matrix(alpha,n*p,1,byrow=TRUE)
        for(i in 1:n){
          S_sq[i,]=sigmasq*(t(DXtX_Big_ind)%*%D_long[,i] + 1/sigmabeta_sq)^(-1)
        }
        
        S_sq_vec=matrix(t(S_sq),n*p,1)
        
        for(i in 1:n){
          
          y_XW=y_long_vec*X_vec*D_long[,i]
          y_XW_mat=matrix(y_XW,n,p,byrow=TRUE)
          
          X_mu_alpha=X_vec*Mu_vec*alpha_vec
          xmualpha_mat=t(matrix(X_mu_alpha,p,n))%*%(matrix(1,p,p)-diag(rep(1,p)))
          XW_mat=matrix(X_vec*D_long[,i],n,p,byrow=TRUE)*xmualpha_mat
          
          
          
          
          
          
          mu_mat[i,]=(t(y_XW_mat)%*%matrix(1,n,1)-(t(XW_mat)%*%matrix(1,n,1)))*(S_sq[i,]/sigmasq) ### ### CAVI updation of variational parameter mu
        }
        Mu_vec=matrix(t(mu_mat),n*p,1)
        #vec_1=as.matrix(t(Big_diag_mat))%*%(as.matrix(Z)%*%as.matrix(mu0_lambda))  ### dim()   dim(mu0_lambda)
        vec_1=log(true_pi/(1-true_pi))
        
        vec_2=as.matrix(0.5*log(S_sq_vec/(sigmasq*sigmabeta_sq)))
        vec_3=as.matrix(Mu_vec^2/(2*S_sq_vec))
        # (vec_1+vec_2+vec_3)[1:p]
        unlogitalpha=vec_1+vec_2+vec_3
        # thres=10^{-9}
        lthres=logit(thres)
        uthres=logit(1-thres)
        indlarge=which(unlogitalpha > uthres)
        indsmall=which(unlogitalpha < lthres)
        unlogitalpha[indlarge]<-uthres
        unlogitalpha[indsmall]<-lthres
        alpha[which(unlogitalpha>9)]=1
        alpha[which(unlogitalpha<=9)]=1/(1+ exp(-unlogitalpha[which(unlogitalpha<=9)])) ### ### CAVI updation of variational parameter alpha
        
        
        e=0
        for(i in 1:n){
          e=e+ELBO_calculator(y,X_mat,S_sq[i,],mu_mat[i,],alpha_mat[i,],sigmasq, sigmabeta_sq, true_pi, D[,i],n,p )
        }
        ELBO_LB= e
        # # print(ELBO_LB)
        # # change_ELBO_LB = ELBO_LB - ELBO_LB_int
        # # print(change_ELBO_LB)
        #
        ###################################################################################
        alpha_new <- alpha
        change_alpha <-alpha_new - alpha_int
        
       
        
        ELBO_LBit[iter]=ELBO_LB
        iter=iter+1
      }
      ELBO_LBit=ELBO_LBit[1:(iter-1)]
      list(var.alpha=alpha, var.mu=mu_mat, var.S_sq=S_sq, var.mu0_lambda=mu0_lambda, var.Sigma0_lambda=Sigma0_lambda, var.elbo=ELBO_LB,var.elboit=ELBO_LBit)
    }
    
    
    
    
    
    candL=seq(0.1,0.9,.2)# Different values of hyperparameter true_pi
    #candL=0.5
    like=rep(0,length(candL))
    elb=like
    
    est_pi=rep(0,n)
    est_q=est_pi
    beta_matr=matrix(0,n,p)
   
    ####################tuning hyperparameters using Carbonetto Stephens##################################
    idmod=varbvs(X_mat,y,Z=Z[,1],verbose=FALSE)
    sigmasq=mean(idmod$sigma)
    pi_est=mean(1/(1+exp(-idmod$logodds)))
    sigmavec=c(0.01,0.05,0.1,0.5,1,3,7,10)
    elb1=matrix(0,length(sigmavec),1)
    for(j in 1:length(sigmavec)){
      res=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq,sigmavec[j],lambda_mean,lambda_var,0,pi_est)
      elb1[j]=res$var.elbo
      
    }
    sigmabeta_sq=sigmavec[which.max(elb1)]
    
    result=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq, sigmabeta_sq,lambda_mean,lambda_var,0,pi_est)
    incl_prob=result$var.alpha
    mu0_val=result$var.mu0_lambda
    
    #
    
    
    
    heat_alpha=matrix(incl_prob,n,p,byrow=TRUE)
    mylist[[resp_index]]=heat_alpha
  }
  alph=matrix(0,p+1,p+1)
  #For first subject
  SUBJECT=20
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0
  
  heat_alpha=alph
  
  
  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
      a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }
  
  heat_alpha=a
  
  
  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))
  
  selection1=1*(heat_alpha>0.5)
  
  
  sensitivity_20[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_20[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)
  
  SUBJECT=90
  
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0
  
  heat_alpha=alph
  
  
  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
      a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }
  
  heat_alpha=a
  
  
  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))
  
  selection1=1*(heat_alpha>0.5)
  
  
  sensitivity_90[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_90[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)
  
  SUBJECT=160
  
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0
  
  heat_alpha=alph
  
  
  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
      a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }
  
  heat_alpha=a
  
  
  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))
  
  selection1=1*(heat_alpha>0.5)
  
  
  sensitivity_160[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_160[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)
  for(SUBJECT in 1:n){
    alph=matrix(0,p+1,p+1)
    
    
    for(i in 1:(p+1)){
      alph[i,-i]=mylist[[i]][SUBJECT,];
    }
    in_pr_13[overiter,SUBJECT]=alph[1,3]
    in_pr_12[overiter,SUBJECT]=alph[1,2]
  }
  
}
#################BELOW IS FOR VISUALIZATION#####################

SUBJECT=20
for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
diag(beta)=0

data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence Structure")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig


heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[j,i]=a[i,j]
  }
}

#heat_alpha=a

data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability for Subject 20")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 20")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig
###########
SUBJECT=90
for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
diag(beta)=0

data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence Structure")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig


heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[j,i]=a[i,j]
  }
}

#heat_alpha=a

data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability For Subject 90")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 90")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig


###########
SUBJECT=160

for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
diag(beta)=0
data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence Structure")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[j,i]=a[i,j]
  }
}

#heat_alpha=a
data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability for Subject 160")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()
fig

alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

#selection1=1*(heat_alpha>alphvec[kmeans(alphvec,c(min(alphvec),max(alphvec)))$size[1]])
selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") + 
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))  
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 160")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   
                   panel.background = element_blank())

fig=fig+ coord_equal()
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  
                  axis.title = element_text(size=30, face = "bold"))
fig=fig + theme(plot.title = element_text(hjust = 0.5))
fig = fig + theme( legend.title = element_text(face = "bold",size = 25), 
                   
                   legend.text = element_text(face="bold",size = 25),
                   
                   legend.key.size = unit(2,'lines'))

fig

data=melt(t(sensitivity_20))
fig=ggplot(data, aes(x=sensitivity_20)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 20")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_20))
fig=ggplot(data, aes(x=specificity_20)) + 
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 20")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig

data=melt(t(sensitivity_90))
fig=ggplot(data, aes(x=sensitivity_90)) + 
  geom_histogram(aes(y=..density..), colour="black",bins=10, fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 90")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_90))
fig=ggplot(data, aes(x=specificity_90)) + 
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 90")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig
data=melt(t(sensitivity_160))
fig=ggplot(data, aes(x=sensitivity_160)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 160")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_160))
fig=ggplot(data, aes(x=specificity_160)) + 
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 160")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"), 
                  axis.title = element_text(size=30, face = "bold"))
fig

al=matrix(0,1,n)

for(SUBJECT in 1:n){
  alph=matrix(0,p+1,p+1)
  
  
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  al[SUBJECT]=alph[1,3]
}

grd=seq(1,n,1)
grd=matrix(grd,n,1)
al=matrix(al,n,1)
data=as.data.frame(cbind(grd,al))
fig=ggplot(data, aes(x=grd, y=al)) +
  geom_line()+
  geom_point()
fig=fig+labs( x=expression(bold("Subject Index")), y=expression(bold("Inclusion Probability")) )
fig=fig +  theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA),
  plot.title = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
  axis.title.x = element_text(size=16, face="bold", colour = "black"),    
  axis.title.y = element_text(size=16, face="bold", colour = "black"),    
  axis.text.x = element_text(size=16, face="bold", colour = "black"), 
  axis.text.y = element_text(size=16, face="bold", colour = "black"),
  strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
  strip.text.y = element_text(size = 14, face="bold", colour = "black"),
  axis.line.x = element_line(color="black", size = 0.3),
  axis.line.y = element_line(color="black", size = 0.3),
  panel.border = element_rect(colour = "black", fill=NA, size=0.3) 
)
fig=fig+ ggtitle(TeX('Inclusion Probability of edge between $x_1$ and $x_3$ versus covariate'))
fig
#save(mylist,file="converged_alpha_per_regression_continuouscov_p4.RData")