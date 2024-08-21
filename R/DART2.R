#library(qlcMatrix)  
library(ape)
library(cubature)
#library(xtable)



find.p <- function(Tp,lTT,alpha=0.05,length.out=1000,prevrej=c(0,0)){
  #### function for finding the rejection threshold on layer \ell/
  m <- length(Tp)
  m0 <- sum(!is.na(Tp))
  mall <- sum(lTT[which(!is.na(Tp))])
  lowlim <- 1/(m0*max(sqrt(log(m0)),1/alpha))
  cl <- lowlim; cr <- 1
  c0=seq(cl,alpha,length.out=length.out)
  Tp <- matrix(Tp,m,1)
  TTc0 <- apply(Tp,1,'<',c0)
  TTc0[,which(is.na(TTc0[1,]))] <- 0
  
  if(is.null(dim(lTT))){
    TTc0 <- TTc0%*%lTT
  }else{
    TTc0 <- TTc0%*%diag(lTT)
  }
  
  estfdr1 <- prevrej[1]+mall*c0
  estfdr2 <- prevrej[2]+(pmax(apply(TTc0,1,sum,na.rm=TRUE),1))
  indsmall <- which(estfdr1/estfdr2<=alpha)
  if(length(indsmall)==0){
    out=cl;indx=1
  }else{
    indx=max(indsmall)
    out=c0[indx]
  }
  return(list(out=out,prevrej=c(estfdr1[indx],estfdr2[indx])))
}



#' DART2 function based on existed aggregation tree. This function does not construct the aggregation tree. Instead, the aggregation tree must be provided as an input parameter to the function in the form of a list of matrices.
#'
#' @importFrom stats uniroot pchisq splinefunH na.omit integrate pchisq pnorm qnorm dhyper glm model.matrix
#' @import tidyverse dplyr
#' @import Matrix
#' @importFrom reshape2 colsplit
#' @importFrom data.table ":=" uniqueN data.table setkey setDT rleid
#' @param alpha (numeric) Vector of desired FDR, default value is 0.05.
#' @param Llist (list of matrices) Aggregation tree presented in the form of matrix list.
#' @param T1 (vector) Input test statistics
#'
#' @return (list of vectors) Index of rejected nodes on each layers
#' @export
DART2_list <- function(alpha=0.05,Llist,T1){
  #### MAIN function of DART2 framwork when the Aggregation tree is existed in the form of a list of matrix.
  ## arguments:
  # alpha: (numeric) Vector of desired FDR
  # Llist: (list of matrixes) Aggregation tree presented in the form of a list of matrix.
  # T1: (vector) Input Z-values transformed from p-values
  ## output: (list of vectors)
  # Index of rejected nodes on each layers
  Rej <- NULL
  nnodes <- m <- length(c(T1))
  S=length(Llist)
  TT <- ifelse(T1==Inf,100,
               ifelse(T1==-Inf,-100,T1))
  
  pval1 <- 2*pnorm(abs(TT),lower.tail=FALSE)
  
  layer1.stats=find.p(pval1,rep(1,m),alpha=alpha,prevrej=c(0,0))
  
  cl1=layer1.stats$out
  prevrejs <- layer1.stats$prevrej
  
  rejindex <- which(pval1<cl1)
  Rej <- data.frame("Feature"=rejindex,"Layer"=1,"Size"=1,
                    "Node"=seq(rejindex),
                    "That"=qnorm(cl1/2,lower.tail=FALSE))
  
  rejs2 <- NULL
  s=2
  J22 <- NULL # store the dynamic number of child 
  LI32 <- diag(rep(1,m))
  while(s<=S){
    I3 <- rep(1,ncol(Llist[[s-1]]))
    I3[rejindex] <- 0
    I3[which(J22==0)] <- 0 # delete the node whose dencendent are all rejected
    I3 <- diag(I3) # diagonal indicates the existing nodes on layer s-1 (1 exist, 0 detected)
    
    LI3 <- Llist[[s-1]]%*%I3
    LI32 <- LI3%*%LI32
    
    J22 <- rowSums(LI3)
    
    if(sum(J22>=2)>0){
      testindex <- which(J22>=2) #index of elements that is in the random aggregation tree
      obsval3 <- LI32%*%t(TT)
      
      J32 <- rowSums(LI32)
      
      pval3 <- 2*pnorm(abs(obsval3[testindex]/sqrt(J32[testindex])),lower.tail=FALSE)
      
      layer3.stats=find.p(pval3,J32[testindex],
                          alpha=alpha,
                          prevrej=c(0,0))
      
      cl3 <- layer3.stats$out
      prevrejs <- layer3.stats$prevrej
      rejindex <- testindex[pval3<cl3]
      
      if(length(rejindex) %in% c(0,1)){
        Rej3 <- NULL
        if(length(rejindex)==1){
          rejm <- LI32[rejindex,]
          Rej3 <- data.frame("Feature"=which(rejm==1),"Layer"=s,"Size"=J32[rejindex],"Node"=1,"That"=qnorm(cl3/2,lower.tail=FALSE))
        }
      }else{
        rejm <- LI32[rejindex,]
        rejm2 <- t(t(rejm)*seq(ncol(rejm)))
        rejsize <- rejm*J32[rejindex]
        rejnode <- rejm*seq(rejindex)
        Rej3 <- data.frame("Feature"=c(rejm2),"Layer"=s,"Size"=c(rejsize),
                           "Node"=c(rejnode),
                           "That"=qnorm(cl3/2,lower.tail=FALSE))%>%
          filter(Feature!=0)
      }
      Rej <- rbind(Rej,Rej3)
      s=s+1
    }else{
      break
    }

  }
  Rej0 <- Rej%>%mutate(T1=TT[Feature])%>%
    mutate(Test=as.integer(abs(T1)>max(That/sqrt(Size),
                                       qnorm(alpha/2,lower.tail=FALSE))))%>%
    group_by(Node,Layer)%>%
    mutate(Test=as.integer(Test==1|abs(T1)==max(abs(T1))))%>%
    ungroup()
  return(Rej0)
}


#' DART2 function based on ordered hypotheses. The function will automatically construct the aggregation tree based on the order of hypotheses.
#'
#' @importFrom stats uniroot pchisq splinefunH na.omit integrate pchisq pnorm qnorm dhyper glm model.matrix
#' @import tidyverse dplyr
#' @import Matrix qlcMatrix
#' @importFrom reshape2 colsplit
#' @importFrom data.table ":=" uniqueN data.table setkey setDT rleid
#' @param alpha (numeric) Vector of desired FDR, default value is 0.05.
#' @param out (data.frame) Data frame with variables 1) "hypo": index of hypothesis; 2) "pval": p-values for each hypothesis; 3)"z": test statistics for each hypothesis; 4) "orderp": the order of hypothesis based on external information.
#' @param L (numeric) Maximum number of layers of the aggregation tree.
#'
#' @return (data.frame) Data frame informing rejected nodes on each layer.
#' @export
DART2_order <- function(out,alpha=0.05,L=5){
  struct_map <- data.frame(out)%>%filter(!is.na(pval))%>%
    rename(Gene=hypo)%>%
    arrange(orderp)%>%
    mutate(Z=as.numeric(z),pvals=as.numeric(pval),
           L1=rleid(Gene))%>%
    dplyr::select(-c("pval","z"))
  
  if(!("Group"%in%colnames(struct_map))){struct_map$Group <- 1}
  
  for(l in 2:L){
    Lm1 <- paste0("L",l-1)
    Ll <- paste0("L",l)
    Lrl <- paste0("Lr",l)
    struct_map <- struct_map%>%group_by(Group)%>%
      mutate(Lg=rleid(get(Lm1)))%>%
      mutate(!!Ll:=ceiling(Lg/2))%>%
      ungroup()%>%
      mutate(!!Ll:=rleid(paste(Group,get(Ll),sep="_")))%>%
      group_by_at(Ll,.add=TRUE) %>%
      mutate(!!Lrl:=ifelse(uniqueN(get(Lm1))>=2,get(Ll),NA))%>%
      ungroup()%>%dplyr::select(-c("Lg"))
  }
  
  smap <- struct_map%>%mutate(phat=find.p(Tp=pvals,lTT=rep(1,n()),alpha=alpha)$out)
  
  Rej_DART2 <- NULL
  if(min(smap$pvals,na.rm=TRUE)<=unique(smap$phat)){
    Rej_DART2 <- smap%>%filter(pvals<=phat)%>%
      mutate("Layer"=1,"Size"=1,"That"=qnorm(phat/2,lower.tail=FALSE))%>%
      dplyr::select(c("L1","Layer","Size","That"))
  }
  smap <- smap%>%filter(pvals>phat)%>%dplyr::select(-c("phat"))
  
  for(l in 2:L){
    Lm1 <- paste0("L",l-1)
    Ll <- paste0("L",l)
    Lrl <- paste0("Lr",l)
    
    smap <- smap%>%group_by_at(Ll)%>%
      mutate(Zl=sum(Z)/sqrt(n()),
             pl=2*pnorm(abs(Zl),lower.tail=FALSE),
             Size=uniqueN(L1))%>%
      ungroup()
    
    smap1 <- smap%>%dplyr::select(c(Ll,all_of(Lrl),"pl","Size","Zl"))%>%
      distinct_all()%>%
      filter(!is.na(get(Lrl)))%>%
      mutate(phat=find.p(Tp=pl,lTT=rep(1,n()),alpha=alpha/max(smap$Size))$out)%>%
      dplyr::select(-c(Lrl,"pl","Zl","Size"))%>%
      full_join(smap,by=Ll)
    
    if(min(smap1$pl,na.rm=TRUE)<=unique(na.omit(smap1$phat))){
      Rej_DART2 <- smap1%>%filter(pl<=phat&!is.na(phat))%>%
        mutate("Layer"=l,"That"=qnorm(phat/2,lower.tail=FALSE))%>%
        mutate(test=as.integer(abs(Z)>max(That/sqrt(Size),
                                          qnorm(alpha/(2*max(smap$Size)),lower.tail=FALSE))))%>%
        group_by_at(Lrl)%>%
        mutate(test=as.integer(test==1|(abs(Z)==max(abs(Z)))))%>%
        ungroup()%>%
        filter(test==1,!is.na(get(Lrl)))%>%
        dplyr::select(c("L1","Layer","Size","That"))%>%
        bind_rows(Rej_DART2)
      
    }
    smap <- smap1%>%filter(pl>phat|is.na(phat))%>%
      dplyr::select(-c("phat"))
  }  
  
  if(is.null(Rej_DART2)){
    output=NULL
  }else{
    output=Rej_DART2%>%left_join(struct_map)%>%
      mutate(Method="DART2")%>%
      dplyr::select(c("Gene","pvals","Z","Layer","Method"))
  }
  
  
  return(output)
}

