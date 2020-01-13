library(weights)
library(R.utils)

dischist= function (x) {
    hist(x,
         breaks=c((min(x)-0.5):(max(x)+0.5))
         )
}

wtd.dischist=function(x,w) {
    wtd.hist(x,weight=w,
             breaks=c((min(x)-0.5):(max(x)+0.5))
             )
}

## Calculates differences of subsequent columns in data frame df.
coldiffs=function(df) {
    data.frame(
        t(
            diff(
                t(
                    as.matrix(df)
                )
            )
        )
    )
}

## Frequency table of x generated using counts. Multiple occurences of elements are allowed. x and counts must have the same lengths.
aggregateData=function (x,counts) {
    lapply(
        x,
        function(col) aggregate(counts,by=list(col),FUN=function(i) sum(as.numeric(i)))
    )
}

## Collapses what is obtained from aggregateData in one list
collapseAggregateData=function(aggdf) {
    do.call(rbind,aggdf[2:length(aggdf)])
}

## Aggregate data frame: 'Forget' changepoint indices -> Frequencies of locations
locFreq=function(df) {
    aggregate(df[, "x"], by=list(df[,"Group.1"]), FUN=function(x) sum(x))
}

## data.frame(t(sapply(unique(ip201All$Group.1), function(x) list(x,sum(ip201All$x[ip201All$Group.1 == x])))))

## aggregate(ip201All, by=list(ip201All$Group.1), FUN=function(x) sum(x))
## aggregate(ip201All[,"x"], by=list(ip201All[,"Group.1"]), FUN=function(x) sum(x))


## Calculates MLE for data x and weights w. Here, x must not contain
## multiple elements! x and w must have the same lengths
wtd.mle=function(x,w) {
    x[which.max(w)]
}

## Conditionally selects from y elements for which x==MLE
conditionOnMLE=function(y,x,mle) {
    y[x==mle]
}

## Process RJMCMC results: k, p as read from files, nK should have
## counts added with addCounts
binK=function(k, nK) {
    ## Remove burn-in
    burnIn=nK$Iterations[1]
    k=k[ k$Iterations >= burnIn, ]
    k$Iterations=NULL
    myNK=ncol(k)
    counts=nK$counts[nK$n_K==myNK]
    aggregateData(k,counts)
}

k2locFreq=function(k, nK) {
    k=binK(k, nK)
    k=collapseAggregateData(k)
    locFreq(k)
}

## Transforms lists of changepoints to 
kLists2locFreq=function(kList, nK) {
    locFreqList=lapply(kList, function (k) k2locFreq(k, nK) )
    locFreqs=do.call(rbind, locFreqList)
    locFreq(locFreqs)
}

## Process RJMCMC results: k, p as read from files, nK should have
## counts added with addCounts
rjmcmcReport=function(k,p,nK) {
    ## Counts
    ##    nK$counts=c(sapply(list(nK$Iterations),diff),0)
    myNK=ncol(k)-1##nrow(k)
    counts=nK$counts[nK$n_K==myNK]
    ##nKsamples=nK[nK$n_K==myNK,]
    k$Iterations=NULL
    p$Iterations=NULL
    kBinned=aggregateData(k,counts)
    results=data.frame(kMLE=sapply(kBinned,function(col) wtd.mle(col[,1],col[,2])))
    ## Condition on MLE for k: For those p_(i+1) that are sampled
    ## together with k_i find the one that occurs most frequently
    pMLE=sapply(seq(1:(length(results$kMLE))),
        function(i) which.max(table(conditionOnMLE(p[[i+1]],k[[i]],results$kMLE[i]))))
    ##        function(i) which.max(table(conditionOnMLE(p[[i]],k[[i]],results$kMLE[i]))))
    
    results$pMLE_K=as.numeric(names(pMLE))

    ## Still missing: Segment before first changepoint. This will be
    ## added above the already calculated results data frame
    pMLE0 <- which.max(
        table( conditionOnMLE(p[[1]], k[[1]], results$kMLE[[1]] )
              )
    )
    
    pMLE0_K <- as.numeric(names(pMLE0))

    results <- rbind(data.frame(kMLE=0, pMLE_K=pMLE0_K),
                     results
                     )

    results$mode=factor(x=results$pMLE_K<0.35,c(TRUE,FALSE),labels=c("M1","M2"))

    ## Differences
    results$diff=c(diff(results$kMLE),0)
    results
}

selectMode=function(modesList,mode) {
    modesList==mode
}

selectModeIndices=function(modesList,mode) {
    modes=selectMode(modesList,mode)
    seq(along=modes)[modes]
}

selectModeIntervals=function(modesList,mode) {
    indices=selectModeIndices(modesList,mode)
    seqToIntervals(indices)
}

calculateModeIntervals=function(k,modesList,mode) {
    intervals=selectModeIntervals(modesList,mode)
    ## sapply(seq(1:(length(intervals)/2-1)),function(i) k[intervals[i,][[2]]+2][[1]]-k[intervals[i,][[1]]+1][[1]])
    sapply(seq(1:(length(intervals)/2-1)),function(i) k[intervals[i,][[2]]+1][[1]]-k[intervals[i,][[1]]][[1]])
}

interleaveLists=function(m1List, m2List) {
    ## ifelse(length(m2List)>length(m1List),
    ## idx <- order(c(seq_along(m2List), seq_along(m1List)))
    ## (c(m2List,m1List))[idx],
    ## idx <- order(c(seq_along(m1List), seq_along(m2List)))
    ## (c(m1List,m2List))[idx]
    ## )
    if(length(m2List)>length(m1List)) {
        idx <- order(c(seq_along(m2List), seq_along(m1List)))
        (c(m2List,m1List))[idx]
    }
    else {
        idx <- order(c(seq_along(m1List), seq_along(m2List)))
        (c(m1List,m2List))[idx]
    }
}

combineModes=function(k,modesList) {
    M1List=calculateModeIntervals(k, modesList, "M1")
    M2List=calculateModeIntervals(k, modesList, "M2")
    interleaveLists(M1List,M2List)
}

addCounts=function(nK) {
    nK$counts=c(sapply(list(nK$Iterations),diff),0)
    nK
}

reduceReport=function(df) {
    M1intervals=selectModeIntervals(df$mode,"M1")
    M2intervals=selectModeIntervals(df$mode,"M2")
    M1diffs=calculateModeIntervals(df$kMLE,df$mode,"M1")
    M2diffs=calculateModeIntervals(df$kMLE,df$mode,"M2")
    diffs=interleaveLists(M1diffs,M2diffs)
    M1K=head(df$kMLE[M1intervals], length(M1diffs))
    M2K=head(df$kMLE[M2intervals], length(M2diffs))
    results=data.frame(pMLE_k=interleaveLists(rep(0.0,length(M1diffs)),rep(0.75,length(M2diffs))))
    ## if(M1K[1]<M2K[1]) {
    ##results$kMLE=interleaveLists(M1K,M2K)
    ## }
    ## else {
         results$kMLE=interleaveLists(M2K,M1K)
    ## }
    results$diffs=diffs
    results
}

reduceReportNew <- function(df) {
    ## Get run length encoding
    rleModes <- rle(as.character(df$mode))
    ## Generate constant IDs for each run
    rleIDs <- rep(1:length(rleModes$lengths), rleModes$lengths)

    ## Assign zero probability to mode 1 and 0.75 to mode 2
    p <- ifelse(rleModes$values == "M1", 0.0, 0.75)

    ## Choose first index k from each run
    k <- tapply(df$kMLE, rleIDs, FUN=function(x) head(x,1))

    ## Sum differences for runs
    diffAgg <- tapply(df$diff, rleIDs, FUN=sum)

    data.frame(pMLE_k=p, diff=diffAgg, kMLE=k, mode=rleModes$values)
}

## Generate a vector of unique IDs for the runs in vector v
rleUniqueID=function(v) {
    r <- rle(v)
    r$values <- paste0(r$values, sprintf("_%03d", seq_along(r$values)))
    inverse.rle(r)
}

## Aggregate data according to runs in vector v using function FUN
## aggregateRuns=function(name,d,v,f=sum) {
##     r <- rleUniqueID(v)
##     aggregate(d[name] ~ r,data=d,FUN=f)
## }

## nK2Hist=function(nK) {
##     nK=
## }
##print(rjmcmcReport(ip568_10,p569_10,nKrestart_ip568$counts)
