est.c <-
function(data1, data4, max1=5, max4=5)
{
  M=0;Dmean=0;
  dat1=apply(data1,1,max)
  dat4=apply(data4,1,max)
  loc41=which(dat1 <= max1 & dat4 <= max4)

  for (k in 1:dim(data1)[1])
  {
    n=dim(data1)[2] 
    L1=data1[k,];
    L4=data4[k,];    
    
    Diff=L4-L1;
    # add a small constant
    D41_upp=c(Diff[2:n])
    D41_low=c(Diff[1:(n-1)])
    M[k]=sum(D41_upp*D41_low)/n
    Dmean[k]=mean(Diff)
  }   

tao=mean(M[loc41])
return(tao)
}
