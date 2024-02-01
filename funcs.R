getcov<-function(SEV){
  options(warn=-1)
  Sigma<-matrix(rep(NA,9),ncol=3)
  while(is.na(Sigma[1,1]==TRUE)){
    if (((3*SEV)/(4*pi))^2/3 < 1){
      lambda1<-runif(1,1,((3*SEV)/(4*pi))^2)
    }else {
      lambda1<-runif(1,((3*SEV)/(4*pi))^2/3,((3*SEV)/(4*pi))^2)
    }
    lambda2<-runif(1,(((3*SEV)/(4*pi))^2)/lambda1,((3*SEV)/(4*pi)))
    lambda3<-(((3*SEV)/(4*pi))^2)/(lambda1*lambda2)
    u<-diag(c(lambda1,lambda2,lambda3))
    Q<-qr.Q(qr(runif(3,-0.8,.8)),complete=TRUE)
    Sigma <- Q %*% u %*% t(Q)
  }
  options(warn=0)
  return(Sigma)
}

SEV.function<-function(sigma) {
  if (eigen(sigma)$values[3] > 0)
  {
    axises<-sqrt(eigen(sigma)$values)
  }
  else {
    axises<-NA
  }
  out<-list()
  out$a<-axises[1]
  out$b<-axises[2]
  out$c<-axises[3]
  out$SEV<-prod(axises,pi,4/3)
  return(out$SEV)
}

overlap<-function(mu1, mu2, sigma1, sigma2){
  axises<-sqrt(eigen(sigma1)$values)
  
}

pairwise<-function(ngroups){
  groups<-seq(1,ngroups)
  max.compare<-((ngroups-1)*ngroups)/2
  i=0
  counter=ngroups
  pairs<-c(0,0)
  while (counter>0){
    i=1+i
    counter=counter-1
    if (counter!=0){
      temp1<-rep(groups[i],counter)
      temp2<-seq(i+1,ngroups)
      pairs<-rbind(pairs,cbind(temp1,temp2))
      if (i==1){
        pairs<-pairs[-1,]
      }
    }}
  colnames(pairs)<-c('group1','group2')
  return(pairs)
}

CubeCount3D <- function(mu1,mu2,S1,S2,HowFine=1,ShowProgress=FALSE){
  # set-up stuff:
  MVNdist <- function(mu,S){function(x){t(x-mu)%*%solve(S)%*%(x-mu)}}
  OneStep <- 10^(-HowFine);
  CubeVol <- OneStep^3;
  CubeCt <- 0;
  A1  <- eigen(S1)$vect
  dD1 <- sqrt((eigen(S1)$val)^2)
  MVNdist2 <- MVNdist(mu2,S2);
  # the loops:
  zlist <- seq(from=0,to=dD1[3]^.5,by=OneStep)
  zlist <- if(length(zlist)>1){c(rev(-zlist),zlist[2:length(zlist)])}
  for(Z in zlist){
    yBd <- (dD1[2]*(1-Z^2/dD1[3]))^.5
    ylist <- seq(from=0,to=yBd,by=OneStep)
    ylist <- if(length(ylist)>1){c(rev(-ylist),ylist[2:length(ylist)])}
    for(Y in ylist){
      xBd <- (dD1[1]*(1-Y^2/dD1[2]-Z^2/dD1[3]))^.5
      xlist <- seq(from=0,to=xBd,by=OneStep)
      xlist <- if(length(xlist)>1){c(rev(-xlist),xlist[2:length(xlist)])}
      for(X in xlist){
        if(MVNdist2(A1%*%c(X,Y,Z)+mu1)<1){CubeCt=CubeCt+1}
      }
    }
    # This will show a "progress report" in the console window as this script is running:
    # (0=start, 1=end)
    if(ShowProgress){print((.5*(Z+dD1[3]^.5))/dD1[3]^.5)}
  }
  return(CubeCt*CubeVol)
} 

#Calculate euclidian distance between two 3 demensional points
distance<-function(a,b){
  l<-sqrt((a[,1]-b[,1])^2+(a[,2]-b[,2])^2+(a[,3]-b[,3])^2)
  return(l)
}

#calculate the mode of a posterior distrobution
mode <- function(s) {
  d <- density(s,na.rm=TRUE)
  d$x[which.max(d$y)]
}
