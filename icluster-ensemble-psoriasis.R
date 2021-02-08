######################################################################
# read in the data, save the lncRNA expression proiles as exprs85034 
# and pData as c85034 and generate the response indicator as pasi75. 
######################################################################

exprs.lnc<-exprs85034
time<-ifelse(c85034$timepoint=="NL", 1, ifelse(c85034$timepoint=="LS", 2, ifelse(c85034$timepoint=="WK 1", 3, ifelse(c85034$timepoint=="WK 2", 4,  ifelse(c85034$timepoint=="WK 4", 5, 6)))))

num=10000
pre<-NULL 
junk.mRNA<-NULL
junk.coef<-NULL 
set.seed(5566)
s.num<-20
for (i in 1:num){
	
	junk<-sample(rownames(exprs.lnc), s.num, replace=FALSE)
	junk.mRNA<-rbind(junk.mRNA, junk)
	exprs.junk<-exprs.lnc[junk, ]
	exp.list.junk<-list()
	exp.list.junk[[1]]<-t(exprs.junk[, which(time==2)])
    exp.list.junk[[2]]<-t(exprs.junk[, which(time==3)])
    exp.list.junk[[3]]<-t(exprs.junk[, which(time==4)])
    exp.list.junk[[4]]<-t(exprs.junk[, which(time==5)])
    
    junk.try<-iCluster(exp.list.junk,k=2, lambda=c(0, 0, 0, 0))
    pre<-rbind(pre, junk.try$clusters) 
    junk.coef<-rbind(junk.coef, unlist(junk.try$beta))
}


floor(30*0.75)

index<-apply(apply(pre, 1, function(x){table(x, pasi75)}), 2, function(x){ifelse(x[2]+x[3]>x[1]+x[4], x[2]+x[3], x[1]+x[4])} )

try<-apply(pre[which(index> floor(30*0.75)
),],1, function(x){junk<-table(x,pasi75); 
	                             if(junk[2]+junk[3]>=junk[1]+junk[4]){return(2-x)}
	                             if(junk[2]+junk[3]<junk[1]+junk[4]){return(x-1)}
                                  }
	                              )

junk.coef1<-junk.coef[which(index> floor(30*0.75)), ]
junk.coef2<-abs(junk.coef1[,1:20])+abs(junk.coef1[,21:40])+abs(junk.coef1[,41:60])+abs(junk.coef1[,61:80])

junk.coef2<-as.vector(junk.coef2)

junk.name<-as.vector(junk.mRNA[which(index>22),])
junk.overall<-tapply(junk.coef2, junk.name, sum)

junk.name<-names(junk.overall[order(junk.overall, decreasing=TRUE)][1:20])

