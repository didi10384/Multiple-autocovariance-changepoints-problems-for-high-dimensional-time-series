#summary of algorithm 
library(foreach)
library(doParallel)
library(rlist)
library(FarmTest)
#Y: multivariate time series: T*p 
#window: window size for comparison 
#step: step size for checkpoints 
#lag.max: maximum lag considered 
#nper: number of permutation for finding threshold
#nblock: number of blocks in block permutation for half window (total blocks are 2*nblock)
#norm.type: norm type for comparison (max, Frobenius, Spectral)

#function:change-point detection in three norms (peak detection of test statistics after smoothing)


our.algorithm=function(Y, window, step, lag.max, nper,nblock){
  
  candidates=statistics.all(Y,lag.max,window,step)
  
  sig.max=after.threshold(Y,candidates$S,window,nper,candidates$maxnorm,"Max",lag.max,nblock)
  sig.F=after.threshold(Y,candidates$S,window,nper,candidates$Fnorm,"Frobenius",lag.max,nblock)
  sig.spe=after.threshold(Y,candidates$S,window,nper,candidates$Snorm,"Spectral",lag.max,nblock)
  
  cp.max=binarySeg(nrow(Y),candidates$S,sig.max,window)
  cp.F=binarySeg(nrow(Y),candidates$S,sig.F,window)
  cp.spe=binarySeg(nrow(Y),candidates$S,sig.spe,window)
  
  return (list(cp.max=cp.max,cp.F=cp.F,cp.spe=cp.spe))
}

#function: use binary segmentation to detect all change points
####input
#Time: length of original time series 
#checkpoints: all checkpoints within original time series
#signals: signals corresponding to checkpoints
#window: window size

####output
#CPs: all change points

binarySeg=function(Time,checkpoints,signals,window){
  
  CPs=NULL                         #all detected change points 
  
  intervals=list(list(chk=checkpoints, 
                      sig=signals, 
                      start=1,
                      end=Time))   #list of all undetected intervals 
  
  ninterval=1                      #number of interval
  
  #while there exist undetected intervals
  while (ninterval > 0){
    
    interval=intervals[[1]]        #current interval
    CP=detection(interval, window) #change point detection for current interval
    intervals=intervals[-1]        #remove current interval from list of undetected intervals
    
    #if there is one change point in current interval, segmentation 
    if (!is.null(CP)){
      
      CPs=c(CPs,CP) 
      two=segment(interval,CP,window)
      
      #if there are valid segmented intervals with length greater than 2*window
      if (length(two)>0){
        for(l in 1:length(two)) {intervals=list.append(intervals,two[[l]])}
      }
    }
    
    ninterval=length(intervals)
  }    
  
  return (CPs)
}

######function: segment an interval into two
######input 
#interval 
#1. chk: location of checkpoints within the interval (at original time series)
#2. sig: signal strength of checkpoints within the interval 
#3. start: location of interval start point at original time series (boundary or previous change point)
#4. end: location of interval end point at original time series (boundary or previous change point)

#CP: location of change point at original time series 

segment=function(interval,CP,window){
  
  intervals=list()                    #list of segmented intervals if their lengths are greater than 2*window 
  
  mid=which.min(abs(interval$chk-CP)) #index of the change point in the checkpoints vector
  
  #left interval 
  if (CP-interval$start>2*window){
    interval1=list(chk=interval$chk[1:mid],
                   sig=interval$sig[1:mid],
                   start=interval$start,
                   end=CP)
    
    intervals=list(interval1)
  }
  #right interval   
  if (interval$end-CP>2*window){
    total=length(interval$chk)        #total number of checkpoints 
    interval2=list(chk=interval$chk[mid:total],
                   sig=interval$sig[mid:total],
                   start=CP,
                   end=interval$end)
    intervals=list.append(intervals,interval2)
  }
  return (intervals)
}

#function: return the location of a change point by finding the checkpoint with maximum signal 
#####input:
#interval 
#1. chk: location of checkpoints within the interval (at original time series)
#2. sig: signal strength of checkpoints within the interval 

#####output
#CP: location of change point at original time series 

detection=function(interval,window){
  #identify checkpoints that have a distance of at least a window size from both boundaries  
  #sequential match between location and signal of one checkpoint
  s=which.min(abs(interval$chk-(interval$start+window)))
  e=which.min(abs(interval$chk-(interval$end-window)))
  
  chk=interval$chk[s:e]
  sig=interval$sig[s:e]
  
  #if the maximum signal is greater than 0, a change point is found, otherwise no change point
  if (max(sig)>0){
    return (CP=chk[which.max(sig)])
  }
}

#outer product of time series Y
outer=function(Y,lag){
  Time =nrow(Y)
  p = ncol(Y)
  Y.1 = Y[1:(Time-lag),]
  Y.2 = Y[(1+lag):Time,]
  
  outer=array(0,c(Time-lag,p,p))
  for (i in 1:(Time-lag)){
    outer[i,,]=Y.1[i,]%*%t(Y.2[i,])
  }
  return (outer=outer)
}

#function: calculate test statistics for each check point before permutation
#results 1. S: locations in the range [window,Time-window]
#results 2. test statistics as the norm difference

statistics.all=function(Y,lag.max,window,step){
  
  Time=nrow(Y)
  #first and last candidates
  first=window
  last=Time-window
  
  #all candidate change point
  S=seq(first,last,by=step)
  maxnorm=matrix(0, nrow=length(S),ncol=lag.max+1)
  Fnorm=matrix(0, nrow=length(S),ncol=lag.max+1)
  Snorm=matrix(0, nrow=length(S),ncol=lag.max+1)
  
  for(i in 1:length(S)){
    s=S[i]
    
    #split one time series into two
    start=s-window+1
    end=s+window
    
    Y1=Y[start:s,]
    Y2=Y[(s+1):end,]
    
    #outer product at each check point for each lag
    for (j in 1:(lag.max+1)){
      
      outer1=outer(Y1,j-1)
      outer2=outer(Y2,j-1)
      
      #sample autocovariance given lag (sample mean of outer product)
      Sigma1=apply(outer1,c(2,3),huber.mean)
      Sigma2=apply(outer2,c(2,3),huber.mean)
      
      maxnorm[i,j]=norm(Sigma1-Sigma2,"m")
      Fnorm[i,j]=norm(Sigma1-Sigma2,"f")
      Snorm[i,j]=norm(Sigma1-Sigma2,"2")
    }
  }
  
  #output
  #S: vector of all checkpoints
  #Norm differences at all checkpoints
  return (candidates=list(S=S,maxnorm=maxnorm, Fnorm=Fnorm, Snorm=Snorm))
}

###function: return index for permuted outer product
perm.index=function(one,two,bsize,nblock){
  #start and end index for the ith block is (i-1)*bsize+1 and i*bsize
  one.start=(one-1)*bsize+1
  one.end=one*bsize
  
  two.start=(two-1)*bsize+1
  two.end=two*bsize
  
  
  id1=seq(one.start[1],one.end[1])  #half index for permuted Sigma1
  for (i in 2:nblock%/%2){
    id1=c(id1,seq(one.start[i],one.end[i]))
  }
  
  
  id2=seq(1,nblock*bsize)           #half index for permuted Sigma2
  id2=setdiff(id2,id1)
  
  id3=seq(two.start[1],two.end[1])  #another half index for permuted Sigma1
  for (i in 2:nblock%/%2){
    id3=c(id3,seq(two.start[i],two.end[i]))
  }
  
  id4=seq(1,nblock*bsize)           #another half index for permuted Sigma2
  id4=setdiff(id4,id3)
  
  return (id=list(l1=id1,r1=id2,l2=id3,r2=id4))
}

###function: return the final test statistics(signal) for each checkpoint 
after.threshold=function(Y,S,window,nper,norm,norm.type,lag.max,nblock){
  
  signal=foreach(i=1:length(S),.combine='c') %dopar% {  #signal strength for each checkpoint
    s=S[i]                                            #location of checkpoint 
    start=s-window+1
    end=s+window
    
    Y1=Y[start:s,]
    Y2=Y[(s+1):end,]
    
    Time=nrow(Y1)                                     #number of time points for half window
    p=ncol(Y1)
    
    permuted=matrix(0,nrow=lag.max+1,ncol=nper)
    
    for(r in 1:nper){
      #for each checkpoint, all lag share the same set of permutations
      set.seed(93*i+r)
      bsize=floor((Time-lag.max)/nblock)            #block size 
      #sample nblock/2 from the left and right window respectively 
      one=sample(seq(1,nblock),size=nblock%/%2)
      two=sample(seq(1,nblock),size=nblock%/%2)
      
      id=perm.index(one,two,bsize,nblock)
      
      for (j in 1:(lag.max+1)){
        #block permutation is conducted at the sequence of outer product for theoretical analysis purpose 
        outer1=outer(Y1,j-1) 
        outer2=outer(Y2,j-1)
        
        #Sigma1 is the mean of half left outer product and half right outer product 
        Sigma1=(apply(outer1[id$l1,,],c(2,3),mean)+apply(outer2[id$l2,,],c(2,3),mean))/2
        Sigma2=(apply(outer1[id$r1,,],c(2,3),mean)+apply(outer2[id$r2,,],c(2,3),mean))/2
        
        if (norm.type=="Max"){
          permuted[j,r]=norm(Sigma1-Sigma2,"m")
        }
        if (norm.type=="Frobenius"){
          permuted[j,r]=norm(Sigma1-Sigma2,"f")
        }
        if(norm.type=="Spectral"){
          permuted[j,r]=norm(Sigma1-Sigma2,"2")
        }
      }
    }
    diff=norm[i,]-apply(permuted,2,max)                 #signal for all lag at ith checkpoint 
    result=max(diff[diff>=0])                           #signal for the ith checkpoint 
    result
  }
  return(signal=signal)
}


