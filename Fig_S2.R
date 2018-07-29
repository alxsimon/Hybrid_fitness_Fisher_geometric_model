# Script to produce Fig. S4

oup <- function(t,n,nu,lambda,sigma,x0){
  dw  <- rnorm(n, 0, sqrt(t/n))
  dt  <- t/n
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + lambda*(nu-x[i-1])*dt + sigma*dw[i-1]
  }
  return(x);
}

d <- 100
s <- 0.1
N <- 1e4

mmar <- c(4,5,1,0.5)
# Mutational variance
v <- 2*s*s

getd <- function(x)
	return(x^2)


for (whichfig in 1:4){
    if(whichfig==1)
    {	
    	p <- rnorm(d+1,mean=0,sd=s)
    	mylims <- c(0,0.25)*1.1
    	fm0 <- 0
    }
    if(whichfig==2)
    {
    	p <- rcauchy(d+1)*5e-3
    	mylims <- c(0,0.25)*1.1
    	fm0 <- 0
    }
    if(whichfig==3)
    {
    	P1p <- oup(t=1,n=d/2,nu=0,lambda=10,sigma=sqrt(s),x0=0.5)
    	P2p <- oup(t=1,n=d/2,nu=0,lambda=10,sigma=sqrt(s),x0=0.5)
    	mylims <- c(0,0.8)*1.1
    	p <- c(P1p[(d/2):1],P2p)
    	fm0 <- 0
    }
    if(whichfig==4)
    {
    	P1p <- cumsum(rnorm(d/2,mean=0,sd=s))
    	P2p <- cumsum(rnorm(d/2,mean=0,sd=s))
    	mylims <- c(0,0.8)*1.1
    	p <- c(P1p[(d/2):1],0,P2p)
    	fm0 <- 1
    }
    
    m <- diff(p)
    z1 <- p[1]
    z2 <- p[d+1]
    # the d mutational changes
    xx <- (0:d)/d
    t <- seq(0,1,length.out=d+1)
    #graphics.off()
    #quartz(width=5,height=3.5*(4/3))
    layout(mat=matrix(c(1,1,2,3),2,2,byrow=T),heights=c(1,0.6))
    mmar <- c(4,5,1,0.5); par(mar=mmar)
    
    plot(range(xx),mylims,col='white',xlab=expression(paste('hybrid index, ', italic(h))),ylab=expression(paste('relative breakdown, ',italic(f[h]))))
    mtext(side=3,adj=0,padj=2,paste0(' (',letters[whichfig],')'))
    lines(xx,getd(z1+c(0,cumsum(m)))/var(m)/d)
    a <- matrix(NA,N,d+1)
    for(i in 1:N)
    {
    	rm <- sample(m)
    	a[i,] <- getd(z1+c(0,cumsum(rm)))
    }
    e <- colMeans(a)
    lines(xx,e/d/var(m))
    fmal <- (((1-t)*z1+t*z2)^2)/d/var(m)
    f0 <- t*(1-t)
    #fm0 <- fm0* ((1-t)*(1-t)*z1*z1 + t*t*z2*z2)/d/var(m)
    fm0 <- fm0* 0.5*((1-t)*(1-t) + t*t)
    lines(xx,f0+fm0,col='blue',lty='dotted')
    lines(xx,f0+fmal,col='red',lty='dotted')
    
    hist(m/sd(m),col='lightgrey',border=F,main='',xlab=expression(paste('mutational effect, ',italic(m[i]/root(v)))),ylab='Pr. density',freq=F,xlim=4*c(-1,1),ylim=c(0,dnorm(0,0,1)))
    lines(seq(-4,4,0.01),dnorm(seq(-4,4,0.01),mean=0,sd=1),col='blue',lty='dotted')
    box()
    
    V <- apply(a,MARGIN=2,FUN=var); V <- V/d/var(m)/d/var(m)
    plot(range(xx),range(V)*1.1,col='white',xlab=expression(italic(h)),ylab=expression(paste('Var(',italic(f[h]),')')))
    lines(xx,V)
    lines(xx,2*f0*(f0+2*fm0),col='blue',lty='dotted')
    lines(xx,2*f0*(f0+2*fmal),col='red',lty='dotted')
}

# > devtools::session_info()
# Session info --------------------------------------------------------------------------------------
# setting  value                       
# version  R version 3.5.0 (2018-04-23)
# system   x86_64, darwin15.6.0        
# ui       RStudio (1.1.447)           
# language (EN)                        
# collate  fr_FR.UTF-8                 
# tz       Europe/London               
# date     2018-06-28                  
# 
# Packages ------------------------------------------------------------------------------------------
# package   * version date       source        
# base      * 3.5.0   2018-04-24 local         
# compiler    3.5.0   2018-04-24 local         
# datasets  * 3.5.0   2018-04-24 local         
# devtools    1.13.5  2018-02-18 CRAN (R 3.5.0)
# digest      0.6.15  2018-01-28 CRAN (R 3.5.0)
# graphics  * 3.5.0   2018-04-24 local         
# grDevices * 3.5.0   2018-04-24 local         
# memoise     1.1.0   2017-04-21 CRAN (R 3.5.0)
# methods   * 3.5.0   2018-04-24 local         
# stats     * 3.5.0   2018-04-24 local         
# tools       3.5.0   2018-04-24 local         
# utils     * 3.5.0   2018-04-24 local         
# withr       2.1.2   2018-03-15 CRAN (R 3.5.0)
# yaml        2.1.19  2018-05-01 CRAN (R 3.5.0)