#unlink("/Library/Frameworks/R.framework/Versions/3.3/Resources/library/PFDA", recursive = TRUE)
#' Plotting 4 RKHS smoothing mean and its DP version simultaneously
#'
#' This function plots 4 different figures which the original data, the RKHS smoothing mean and its differential private
#' version are simply including in each of them
#' @param DP output of \link{DP.RKHS} or \link{DP.Sim} list
#' @export
plotMDP=function(MDP,text.main="Default",
              legend.cex="Default",legend.loc="Default",seg.lin=0.75,
              text.font=2,text.width=0.15,legend.size=0.57,xlab="", ylab="",extra_range="Default",
              lty.p=3,lwd.p=3,par="Default"){

  if(is.character(text.main)==TRUE){text.main=rep("",MDP$M)}
  if(is.character(legend.cex)==TRUE){legend.cex=rep(1.4,MDP$M)}
  if(is.character(legend.loc)==TRUE){legend.loc=rep("topleft",MDP$M)}
  if(is.character(extra_range)==TRUE){extra_range=matrix(0,MDP$M,2)}

  if(is.character(par)==TRUE){
    par.row=floor(sqrt(MDP$M))
    par.col=floor((MDP$M+par.row-1)/par.row)
  }
  else{
    par.row=par[1]
    par.col=par[2]
  }

  if(length(Reduce(intersect,list(length(text.main),length(legend.cex),length(legend.loc))))==1){
  par(mfrow=c(par.row,par.col),font=2,cex=legend.size)
  for(i in 1:(MDP$M)){
    plotDP(DP=MDP[[i]],
       text.main=text.main[i],legend.cex=legend.cex[i],legend.loc=legend.loc[i],seg.lin=seg.lin,
       text.font=text.font,text.width=text.width,xlab=xlab,ylab=ylab,extra_range=extra_range[i,],
       lty.p=lty.p,lwd.p=lwd.p)
  }
  }
    else{
      stop("length of the vectors are not equal")
    }
  }

# library(fda)
# library(refund)
# Data=growth$hgtf
# Data=DTI$cca
# Data=t(Data)
# dim(Data)
# AA=MDP.RKHS(grid = as.numeric(row.names(Data)),Data = Data)
# str(AA)
# plotMDP(MDP = AA,legend.size=0.52,seg.lin=1.1)
