# unlink("", recursive = TRUE)
#' Releasing four Differential Private RKHS smoothing means of a dataset
#'
#' This function create four DP RKHS smoothing means from an existing data set
#' based on any arbitrary 4 by 1 vectors for \code{alpha}, \code{beta}, \code{kernel},
#' \code{phi} and \code{ro}.
#' @param Data a matrix which the of interest curves are located in columns
#' @param grid grid (x-axis) for each curve, default is equally espaced between 0 and 1.
#' @param alpha,beta vector of privacy parameters
#' @param phi vector of penalty parameters
#' @param ro vector of kernel range parameters
#' @param drop.col TRUE/FALSE for cleaning the Data, deleting Columns/Rows including missing values
#' @return four \link{DP.RKHS} output each of them includin "f.tilda","delta","f" and "X"
#' @export
MDP.RKHS=function(grid="FALSE",Data,alpha=rep(1,1),beta=rep(0.1,1),
                 kernel=c("Gau"),phi=rep(0.01,1),ro=rep(0.2,1),col.drop=TRUE){

 if(length(Reduce(intersect,list(length(alpha),length(beta),
                                 length(kernel),length(phi),length(ro))))==1){

  M=length(alpha)
  A=append(as.list(rep(NA,M)),list("alpha"=alpha,"beta"=beta,"kernel"=kernel,"phi"=phi,"ro"=ro,"M"=M))
for(i in 1:M){
     A[[i]]=DP.RKHS(grid=grid,Data=Data,alpha =alpha[i],beta = beta[i],kernel=kernel[i],
                    phi=phi[i],ro=ro[i],col.drop=col.drop)
}
  names(A)=c(paste("DP", 1:M, sep = ""),"alpha","beta","kernel","phi","ro","M")
  return(A)

}
  else{
    stop("length of the vectors are not equal")
  }
  }

# Y=MDP.RKHS(Data = Data,alpha=rep(1,4),beta = rep(0.1,4),kernel = c("Exp","M3/2","Gau","M5/2"),
#            phi = rep(0.01,4),ro = rep(0.2,4))
# plotMDP(MDP = Y,extra_range = matrix(rep(c(0,0.4),4),nrow=4))
