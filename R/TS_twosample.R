TS_twosample <-
function(data1,data4,tao,band,quant)
{
  n=dim(data1)[2];
  hwidth=band/n;
  x=c(1:n)/n
  Snw <- matrix(0, nrow = n, ncol = n)
  ## Calculate the hat matrix for the Nadaraya-Watson kernel estimator
  In <- diag(rep(1, n)) ## identity matrix
  for(j in 1:n){
    y <- In[,j]
    Snw[,j] <- ksmooth(x, y, kernel = "normal", bandwidth = hwidth, x.points = x)$y
  }
 
  ad_df=sum(diag(Snw))
  
  
  Sev=0;Suv=0;Dsum=0;Deql=0;Dnun=0;
  sigma1=0;sigma4=0;M=0;sigma41=0;
  sigma_1_sq=0;sigma_4_sq=0;sigma1_fg=0;
  sigma4_fg=0;sig4sig1=0;sigma_unequal=0;
  Ts_yvec=0;pairT=0;pairT.sigma=0;  
  sigma=0;Xg=0;
  
  for (k in 1:dim(data1)[1])
  {
    n=dim(data1)[2] 
    L1=data1[k,];
    L4=data4[k,];  

    # kernel smoothing test
    Diff=L4-L1;  
    d41=ksmooth(x,Diff,"normal",bandwidth=hwidth)  
    Xg[k]=sqrt(sum((d41$y-Diff)^2)/(n-ad_df))
    Ts_yvec[k]=mean((d41$y)^2) 
    
    # Assumption Free nonparametric test
    D41_upp=c(Diff[2:n])
    D41_low=c(Diff[1:(n-1)])
    M[k]=sum(D41_upp*D41_low)/(n-1)
    
    # equal variance
    sigma1[k]=mean((L1[2:(length(L1))]-L1[1:(length(L1)-1)])^2)/2
    sigma4[k]=mean((L4[2:(length(L4))]-L4[1:(length(L4)-1)])^2)/2
    sigma41[k]=(sigma1[k]+sigma4[k])^2+4*M[k]*(sigma1[k]+sigma4[k])
    
    # unequal variance
    sigma_1_sq[k]=sum((L1[2:(length(L1)-2)]-L1[1:(length(L1)-3)])^2*(L1[4:(length(L1))]-L1[3:(length(L1)-1)])^2)/(4*(n-3))
    sigma_4_sq[k]=sum((L4[2:(length(L4)-2)]-L4[1:(length(L4)-3)])^2*(L4[4:(length(L4))]-L4[3:(length(L4)-1)])^2)/(4*(n-3))    
    sigma1_fg[k]=sum((L4[1:(n-2)]-L1[1:(n-2)])*(c(L4[1],L4[1:(n-3)])-c(L1[1],L1[1:(n-3)]))*(L1[3:n]-L1[2:(n-1)])^2)/(2*(n-3))
    sigma4_fg[k]=sum((L4[1:(n-2)]-L1[1:(n-2)])*(c(L4[1],L4[1:(n-3)])-c(L1[1],L1[1:(n-3)]))*(L4[3:n]-L4[2:(n-1)])^2)/(2*(n-3))
    sig4sig1[k]=sum((L1[2:n]-L1[1:(n-1)])^2*(L4[2:n]-L4[1:(n-1)])^2)/(4*(n-3))          
    sigma_unequal[k]=(sigma_1_sq[k]+4*sigma1_fg[k])+(sigma_4_sq[k]+4*sigma4_fg[k])+2*sig4sig1[k]
    
    Sev[k]=sigma41[k]  
    Suv[k]=sigma_unequal[k]
  }
  
  Suv[which(Suv <= 0)] = min(Suv[which(Suv > 0)])
  Sev[which(Sev <= 0)] = min(Sev[which(Sev > 0)])
  Xg[which(Xg <= 0)] = min(Xg[which(Xg > 0)])

  Sev_a0=sqrt(Sev);
  Suv_a0=sqrt(Suv);
  CHQBC_1_adjB=Xg;

  if (quant[1] > 0)
  {
    Sev_a0=sqrt(Sev)+quantile(sqrt(Sev), quant[1])
  }
  if (quant[2] > 0)
  {
    Suv_a0=sqrt(Suv)+quantile(sqrt(Suv), quant[2])
  }
  if (quant[3] > 0)
  {
    CHQBC_1_adjB=Xg+quantile(Xg, quant[3])
  }
  
Dsum=M
Deql=(Dsum-tao)/(Sev_a0/sqrt(n-1))
Dnun=(Dsum-tao)/(Suv_a0/sqrt(n-1))

Tsb = 1/(2*(band*0.37)*sqrt(pi));
vu=sqrt(2*pi)/(2*pi*n*(band*0.37))
Amax=Snw%*%Snw;
eigenvalue=as.numeric(eigen(Amax)$values)
d= sum(eigenvalue)^2/sum(eigenvalue^2)
delta=  sum(eigenvalue)/(n*d)  
  
Test.adj=(Ts_yvec)/(CHQBC_1_adjB^2)
TS_kn=(((Test.adj/(delta*d))^(1/3)-(1-2/(9*d)))/sqrt(2/(9*d)))

return(a=list("sigma1"=sigma1,"sigma4"=sigma4,"TS_kn"=TS_kn,"Ts_yvec"=Ts_yvec, "Dsum"=Dsum,"Deql"=Deql, "Dnun"=Dnun, "Sev"=Sev, "Suv"=Suv, "Xg"=Xg))
}
