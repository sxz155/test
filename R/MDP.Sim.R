# unlink("/Library/Frameworks/R.framework/Versions/3.2/Resources/library/PFDA", recursive = TRUE)
#' Releasing 4 Differential Private RKHS smoothing means of a simulated dataset
#'
#' This function create 4 DP RKHS smoothing means from an existing data set
#' with known eigenvalues and eigenvectors
#' @param alpha,beta Privacy parameters, real numbers
#' @param N real vector 4*1, number of observations
#' @param n real vector 4*1, number of grid points
#' @param e.val.x real vector n*1, eigenvalues
#' @param e.vec.x real valued matrix n*N, eigenvectors
#' @param tau range of the uniform distribution in KL expansion
#' @param phi real number, penalty parameter
#' @param mu real vector n*1, initial mean vector
#' @param e.val.z real valued matrix n*N, eigenvectors of noise
#' @param pow smoothing parameter, e.val.x_{i}=i^{-pow}
#' @param ro range parameter in kernel, real number
#' @return four \link{DP_Sim} output each of them includin "f.tilda","delta","f" and "X"
#' @export
DP4.Sim=function(alpha=rep(1,4),beta=rep(0.1,4),kernel=c("Exp","M3/2","M5/2","Gau"),
              n=rep(100,4),N=rep(25,4),
              phi=rep(0.01,4),ro=rep(0.2,4),tau=rep(0.4,4),
              case=rep(2,4),pow=rep(4,4),mu=matrix(rep(0,4*100),ncol = 4)){
if(length(Reduce(intersect,list(length(alpha),length(beta),length(kernel),
                                length(n),length(N),
                                length(phi),length(ro),length(tau),
                                length(case),length(pow),dim(mu)[2])))==1){
  M=length(alpha)
  A=list("DP1"=NA,"DP2"=NA,"DP3"=NA,"DP4"=NA,
         "alpha"=alpha,"beta"=beta,"kernel"=kernel,"phi"=phi,
         "ro"=ro,"pow"=pow,"mu"=mu,"n"=n,"N"=N,"tau"=tau,"case"=case,"pow"=pow)
  for(i in 1:M){
    A[[i]]=DP.Sim(alpha =alpha[i],beta = beta[i],kernel=kernel[i],phi=phi[i],ro=ro[i],
              n=n[i],N=N[i],tau=tau[i],case=case[i],pow=pow[i],mu=mu[,i])
  }
  return(A)
}
else{
  stop("length of the vectors are not equal")
}
}
