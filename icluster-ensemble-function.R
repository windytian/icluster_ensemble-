
icluster_ensemble<-function(lnc.ada, time, num, pasi.ada, s.num, cut1, cut2){

pre<-NULL 
junk.mRNA<-NULL
junk.coef<-NULL 
#set.seed(5566)

for (i in 1:num){
	
	junk<-sample(rownames(lnc.ada), s.num, replace=FALSE)
	junk.mRNA<-rbind(junk.mRNA, junk)
	exprs.junk<-lnc.ada[junk, ]
	exp.list.junk<-list()
	time.u<-length(unique(time)) 
	
	for (kk in 1:time.u){
	   exp.list.junk[[kk]]<-t(exprs.junk[, which(time==kk)])
       }
       
    junk.try<-iCluster(exp.list.junk,k=2, lambda=rep(0, time.u))
    pre<-rbind(pre, junk.try$clusters) 
    junk.coef<-cbind(junk.coef, unlist(junk.try$W))
                }





   index<-apply(apply(pre, 1, function(x){table(x, pasi.ada)}), 2, function(x){ifelse(x[2]+x[3]>x[1]+x[4], x[2]+x[3], x[1]+x[4])} )


   junk.coef1<-t(junk.coef)[which(index>cut1), ]
   junk.coef2<-NULL 
   
   for (kk in 0:(time.u-1){
    junk.coef2<-junk.coef2+abs(junk.coef1[,(kk*s.num+1):((kk+1)*s.num)])
    }


junk.coef2<-as.vector((junk.coef2)) 

junk.name<-as.vector(junk.mRNA[which(index>cut1),] )
junk.overall<-tapply(junk.coef2, junk.name, sum)

names(junk.overall[order(junk.overall, decreasing=TRUE)][1:cut2])

 } 




#LOO  using the defined function (iCluster-ensemble)
pred<-NULL 

for (j in 1:length(pasi.ada)){
	     
	    lnc.loo<-lnc.ada[, which(id!=id[j]) ]
	    time.loo<-time[which(id!=id[j])]
	    pasi.loo<-pasi.ada[-j]
	    set.seed(5566)
	    sel.lnc<-icluster_ensemble(lnc.loo, time.loo, 10000, pasi.loo, 20, 10,20)

	   
	    x.sel<-(exprs.junk1[sel.lnc,]+exprs.junk2[sel.lnc, ]+exprs.junk3[sel.lnc, ]+exprs.junk4[sel.lnc, ])/4 #the average 
	
	    fit<-svm(t(x.sel[,-j]), pasi.loo)
		pred<-c(pred, predict(fit, t(x.sel), type="class")[j]) 
		
		
		}
		


table(ifelse(pred>0, 1, 0), pasi.ada)















