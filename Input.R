#############################################################################
################## Input of Image ###########################################

Input <- function(y,k=1)
{
  n=dim(y)[1]
  m=dim(y)[2]
  y5=y;
  k2=k;
  i=0;
  while(k>0)
  {
    y1=y5;
    y2=y5;
    y3=y5;
    y4=y5;
    
    lm1=lm1_Image(y5[max(1,(n-k2)/2-2*k2+1+i):((n-k2)/2+i),((m-k2)/2+1):((m+k2)/2)])
    lm2=lm1_Image(t(y5[((n-k2)/2+1):((n+k2)/2),max(1,(m-k2)/2-2*k2+1+i):((m-k2)/2+i)]))
    lm3=lm1_Image(t(y5[((n-k2)/2+1):((n+k2)/2),min(m,(m+k2)/2+2*k2-i):((m+k2)/2+1-i)]))
    lm4=lm1_Image(y5[min(n,(n+k2)/2+2*k2+1-i):((n+k2)/2+1-i),((m-k2)/2):((m+k2)/2)])
    
    coef1=lm1$coef
    coef2=lm2$coef
    coef3=lm3$coef
    coef4=lm4$coef
    for(j in 1:k)
    {
      y1[(n-k)/2+1,(m-k)/2+j] <- coef1[1] + coef1[2]*y1[(n-k)/2,(m-k)/2+j] + coef1[3]*y1[(n-k)/2+1,(m-k)/2+j-1] + coef1[4]*y1[(n-k)/2,(m-k)/2+j-1]
      y2[(n-k)/2+j,(m-k)/2+1] <- coef2[1] + coef2[2]*y2[(n-k)/2+j-1,(m-k)/2+1] + coef2[3]*y2[(n-k)/2+j,(m-k)/2] + coef2[4]*y2[(n-k)/2+j-1,(m-k)/2]
      y3[(n-k)/2+j,(m+k)/2] <- coef3[1] + coef3[2]*y3[(n-k)/2+j-1,(m+k)/2] + coef3[3]*y3[(n-k)/2+j,(m+k)/2+1] + coef3[4]*y3[(n-k)/2+j-1,(m+k)/2+1]
      y4[(n+k)/2,(m-k)/2+j] <- coef4[1] + coef4[2]*y4[(n+k)/2+1,(m-k)/2+j] + coef4[3]*y4[(n+k)/2,(m-k)/2+j-1] + coef4[4]*y4[(n+k)/2+1,(m-k)/2+j-1]
    }
    y5[(n-k)/2+1,((m-k)/2+1):((m+k)/2)] <- y1[(n-k)/2+1,((m-k)/2+1):((m+k)/2)]
    y5[((n-k)/2+1):((n+k)/2),(m-k)/2+1] <- y2[((n-k)/2+1):((n+k)/2),(m-k)/2+1]
    y5[((n-k)/2+1):((n+k)/2),(m+k)/2] <- y3[((n-k)/2+1):((n+k)/2),(m+k)/2]
    y5[(n+k)/2,((m-k)/2+1):((m+k)/2)] <- y4[(n+k)/2,((m-k)/2+1):((m+k)/2)]
    
    y5[(n-k)/2+1,(m-k)/2+1] <- mean( c(y1[(n-k)/2+1,(m-k)/2+1], y2[(n-k)/2+1,(m-k)/2+1] ))
    y5[(n-k)/2+1,(m+k)/2] <- mean( c(y1[(n-k)/2+1,(m+k)/2], y3[(n-k)/2+1,(m+k)/2] ))
    y5[(n+k)/2,(m-k)/2+1] <- mean( c(y2[(n+k)/2,(m-k)/2+1], y4[(n+k)/2,(m-k)/2+1] ))
    y5[(n+k)/2,(m+k)/2] <- mean( c(y3[(n+k)/2,(m+k)/2], y4[(n+k)/2,(m+k)/2] ))
    
    k=k-2
    
#    print(k)
  }
  i=i+1
  y5
}



##########################################################################
################ lm de una imagen 3 covariables ##########################

lm1_Image=function(imagen)
{
  tam=dim(imagen)
  Im=as.matrix(imagen)
  Im1=Im[2:tam[1],2:tam[2]]
  Im11=Im[1:(tam[1]-1),2:tam[2]]
  Im12=Im[2:tam[1],1:(tam[2]-1)]
  Im13=Im[1:(tam[1]-1),1:(tam[2]-1)]
  lm_image=lm(c(Im1)~c(Im11)+c(Im12)+c(Im13))
  lm_image
}
