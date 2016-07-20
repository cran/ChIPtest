TS_kernel <-
function(data,band,quantile)
{
  
  sigma=0;Xg=0;
  n=dim(data)[2];
  hwidth=band/n;
  x=c(1:n)/n
  Snw <- matrix(0, nrow = n, ncol = n)
  ## Calculate the hat matrix for the Nadaraya-Watson kernel estimator
  In <- diag(rep(1, n)) ## identity matrix
  for(j in 1:n){
    y <- In[,j]
    Snw[,j] <- ksmooth(x, y, kernel = "normal", bandwidth = hwidth, x.points = x)$y
  }
  #ad_df=2*sum(diag(Snw))-sum(diag(Snw%*%t(Snw)));
  ad_df=sum(diag(Snw))
  
  for (k in 1:dim(data)[1])
  {
    x=c(1:n)/n
    diff41=data[k,]; 
    # kernel smoothing
    hwidth=band/n;
    d41=ksmooth(x,diff41,"normal",bandwidth=hwidth)
    
    #sigma1_ut=mad(y1$y-L1)  #0
    sigma <- sqrt(sum((d41$y-diff41)^2) / (length(diff41)-ad_df))
    Xg[k]=sigma;
  }
  
  Xg[which(Xg==0)]=min(Xg[which(Xg>0)]); 
  
  #
  hist(Xg^2,xlab=expression(paste(hat(sigma)^2)), main="")
  #abline(v=quantile(Xg, quantile),col=2)
  CHQBC_1_adjB=Xg;
  
  if (quantile > 0)
  {
    CHQBC_1_adjB=Xg+quantile(Xg, quantile)
  }
  
  Ts_yvec=0;
  Sign=0;
  # Start Kernel Smoothing #
  
  
  for (k in 1:dim(data)[1])
  {
    x=c(1:n)/n
    # normalized variance
    diff=data[k,]/CHQBC_1_adjB[k];  
    
    hwidth=band/n;   
    #y1=ksmooth(x,L1,"normal",bandwidth=hwidth)
    #y4=ksmooth(x,L4,"normal",bandwidth=hwidth)
    d41=ksmooth(x,diff,"normal",bandwidth=hwidth)
    diffy4y1=d41$y
    ####  Build the test statistics #####
    # minimax test
    #library(kernlab)
    #rbfdot(sigma=1)  
    ytest.vec=diffy4y1
    Ts_yvec[k]=mean(ytest.vec^2)
    Sign[k]=sign(sum(ytest.vec))
  }
  
  Tsb = 1/(2*(band*0.37)*sqrt(pi));
  vu=sqrt(2*pi)/(2*pi*n*(band*0.37))
  #delta=1/(sqrt(2)*n)
  #d=1/(sqrt(2*pi)*band*0.37/n)
  Amax=Snw%*%Snw;
  eigenvalue=as.numeric(eigen(Amax)$values)
  d= sum(eigenvalue)^2/sum(eigenvalue^2)
  delta=  sum(eigenvalue)/(n*d)  
  
  Test.adj=Ts_yvec
  TS_kn=(((Test.adj/(delta*d))^(1/3)-(1-2/(9*d)))/sqrt(2/(9*d)))
  
  return(list("TS"=TS_kn,"TS_sign"=Sign, "Tmean"=Ts_yvec))
}
