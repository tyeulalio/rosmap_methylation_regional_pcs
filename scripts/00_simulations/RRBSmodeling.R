library(caTools)

##Functions for modeling, p-value adjustment, and DMR detection

##This function will fit single-site binomial regression models for a simulated dataset
fitmods.sim=function(x=simtab,colstart=5,nsamp=6,group.ind=c(0,0,0,1,1,1)){
    modres=apply(x[,colstart:ncol(x)],1,runmod.sim,group.ind)
    sc.avg=apply(x[,colstart:(colstart+nsamp-1)],1,mean,na.rm=T)
    mean.pm=apply(x[,(colstart+nsamp):ncol(x)],1,mean,na.rm=T)
    sum.not.na=function(x){sum(!is.na(x))}
    n1=apply(x[,(colstart+nsamp):ncol(x)][group.ind==0],1,sum.not.na)
    n2=apply(x[,(colstart+nsamp):ncol(x)][group.ind==1],1,sum.not.na)
    restab=data.frame(t(modres),as.numeric(x[,1]),as.numeric(x[,2]),as.numeric(x[,3]),round(sc.avg),round(mean.pm),n1,n2)
    names(restab)=c("pval","pm.diff","chrpos","diffvec","clust.ind","mean.score","mean.pm","n1","n2")
    return(restab)
}

##Single-site model fitting subroutine
runmod.sim=function(x,group.ind=c(0,0,0,1,1,1))
{
tryCatch(
{
xmat=as.data.frame(matrix(x,ncol=2))
#print(head(xmat))
svec=round(unlist(xmat[,1])*unlist(xmat[,2])/100)
fvec=unlist(xmat[,1])-svec
sfmat=cbind(svec,fvec)
mod=glm(sfmat~group.ind,family="binomial")
pval=anova(mod,test="Chisq")[2,5]
pm.diff=exp(sum(mod$coef))/(1+exp(sum(mod$coef)))-(exp(mod$coef[1])/(1+exp(mod$coef[1])))
return(c(pval,pm.diff))
}
,error=function(e) rep(NA,2))
}

##QUASASS Adjustment
##The SSQA function will do a full simulation.  For faster computation when working with multiple datasets, the SSQA.pgen function can be used
##to generate a simulation object that can be provided here using the "plist" option
SSQA=function(modmat,nrep=999, plist=NULL)
{
sizemat=matrix(0,nrow(modmat),3)
sizemat[,1]=modmat[,"mean.score"]
sizemat[,2]=apply(modmat[,c("n1","n2")],1,min)
sizemat[,3]=apply(modmat[,c("n1","n2")],1,max)
#print(head(sizemat))
if(length(plist)==0){
plist=SSQA.pgen(sizemat,nrep)}
sizemat.lab=paste(sizemat[,1],sizemat[,2],sizemat[,3],sep=":")
plist.lab=paste(plist[[1]][,1],plist[[1]][,2],plist[[1]][,3],sep=":")
score.map=match(sizemat.lab,plist.lab)
#print(head(score.map))
pmat=plist[[2]]
pval.quant=0
for(i in 1:nrow(modmat)){
if(i%%10000==0){cat("working on row ", i, "\n")}
pval.quant[i]=(sum(pmat[,score.map[i]]<modmat$pval[i])+1)/(nrep+1)}
modmat$p.adj=pval.quant
return(modmat)
}

##Function to compute SSQA null distributions
SSQA.pgen=function(sizemat,nrep){
scprmat.red=unique(sizemat)
cat("number of unique cases = ", nrow(scprmat.red), "\n")
pmat=pdist.loop(scprmat=scprmat.red,nrep=nrep)
return(list(scprmat.red,pmat))
}

##Subroutine of SSQA.pgen
pdist.loop=function(scprmat,nrep){
pmat=apply(scprmat,1,pdist.sim,nrep=nrep)
return(pmat)
}

##Subroutine of SSQA.pgen
pdist.sim=function(scprmat.row,nrep){
sampsize=scprmat.row[1]
n1=scprmat.row[2]
n2=scprmat.row[3]
pvec=replicate(n=nrep,bmod.p(sampsize=sampsize,n1=n1,n2=n2))
return(pvec)
}

##Subroutine of SSQA.pgen
bmod.p=function(sampsize=sampsize,n1=n1,n2=n2)
{
prop=runif(1,min=0.01,max=0.99)
s1=rbinom(n1,sampsize,prop)
s2=rbinom(n2,sampsize,prop)
SF=cbind(c(s1,s2),c(sampsize-s1,sampsize-s2))
#print(SF)
gvec=rep(c(0,1),times=c(n1,n2))
mod=glm(SF~gvec,family="binomial")
pval=anova(mod,test="Chisq")[2,5]
return(pval)
}


##DMR detection
##The function dmrID will find DMRs from a fitmods.sim output matrix based on a the specified minimum CpG density using either raw ("orig") or adjusted ("adj") p-values.

dmrID=function(modmat,lencalc=50,min.density=.01,ptype="adj")
{
if(ptype=="orig"){x=modmat[,"pval"]}
else{x=modmat[,"p.adj"]}
x.dms=sign(modmat[,"pm.diff"])
x.pos=modmat[,"chrpos"]
sigsites=which(!is.na(x) & x<.05)
DMRmat=as.data.frame(matrix(0,length(sigsites),ncol=lencalc))
row.names(DMRmat)=sigsites
DMRmat[,1]=x.dms[sigsites]
lenseq=2:lencalc
DMRmatfunc=function(seqval,x=x,x.dms=x.dms,x.pos=x.pos,sigsites=sigsites,min.density=min.density)
{dmrvec=ifelse(DMRk.log(x,k=seqval)[,2]<log(.05) & PMsign(x.dms,k=seqval)!=0 & calc.dens(x.pos,k=seqval)> min.density, PMsign(x.dms,k=seqval),0)
dmrvec=dmrvec[sigsites]
}
DMRmat[,2:lencalc]=sapply(lenseq,DMRmatfunc,x=x,x.dms=x.dms,x.pos=x.pos,sigsites=sigsites,min.density=min.density)
#DMRmat.end=DMRmat[DMRmat[,1]!=0,]
cat("Number of possible DMR starts ", nrow(DMRmat), "\n")
lenDMR=dmrdir=0
for(i in 1:(nrow(DMRmat)-1)){
if(i%%1000==0){cat("working on row ", i, "\n")}
diffs=as.numeric(row.names(DMRmat)[-(1:i)])-as.numeric(row.names(DMRmat)[i])+1
#cat("diffs = ")
#print(head(diffs))
runlens=which(DMRmat[i,]!=0)
#cat("runlens =")
#print(head(runlens))
dr.int=intersect(diffs,runlens)
runlen=ifelse(length(dr.int>0),max(dr.int),1)
#print(runlen)
lenDMR[i]=runlen
dmrdir[i]=DMRmat[i,runlen]
#print(lenDMR[i])
}
dmr.se=as.data.frame(cbind(starts=as.numeric(row.names(DMRmat)[1:length(lenDMR)]),ends=as.numeric(row.names(DMRmat)[1:length(lenDMR)])+lenDMR-1,direction=dmrdir))
cut1=dmr.se[dmr.se$starts!=dmr.se$ends,]
cut1$pval=apply(cut1,1,calc.logp.dmr,modmat=modmat,ptype=ptype)
keeprow=0
#return(cut1)}
for(i in 1:nrow(cut1)){
keeprow[i]=ifelse(length(which(cut1[,1] <= cut1[i,1] & cut1[,2] >= cut1[i,2]))==1, 1,0)
}
cut2=cut1[keeprow==1,]
cut2$pval=apply(cut2,1,calc.logp.dmr,modmat=modmat,ptype=ptype)
#return(cut2)}
currstart=cut2[1,1]
currend=cut2[1,2]
currdir=cut2[1,3]
cut3=as.data.frame(matrix(0,1,3))
for(i in 2:nrow(cut2)){
#cat("i = ", i, "\n")
if(cut2[i,1] <= currend & cut2[i,3]==currdir){
currdmr=c(currstart,cut2[i,2],currdir)
#cat("case 1: currdmr = \n")
#print(currdmr)
currstart=currdmr[1]
currend=currdmr[2]
currdir=currdmr[3]}
else{
currdmr=c(currstart,currend,currdir)
#cat("case 2: currdmr = \n")
#print(currdmr)
cut3=rbind(cut3,c(currdmr))
currstart=cut2[i,1]
currend=cut2[i,2]
currdir=cut2[i,3]}
}
cut3=rbind(cut3,c(currstart,currend,currdir))
cut3=cut3[-1,]
cut3=unique(cut3)
cut3$dmr.len=modmat[cut3[,2],"chrpos"]-modmat[cut3[,1],"chrpos"]+1
cut3$nsites=cut3[,2]-cut3[,1]+1
site.density=cut3$nsites/cut3$dmr.len
cut3$pval=apply(cut3,1,calc.logp.dmr,modmat=modmat,ptype=ptype)
cut3$mean.diff=apply(cut3,1,calc.mean.diff,modmat=modmat)
cut3$numNA=apply(cut3,1,sum.na,modmat=modmat,ptype=ptype)
#print(site.density)
cut3=cut3[site.density >= min.density,]
dmrtab=as.data.frame(cbind(chr=as.character(modmat[cut3[,1],"chr"]),start=modmat[cut3[,1],"chrpos"],end=modmat[cut3[,2],"chrpos"],len=cut3$dmr.len,nsites=cut3$nsites,NAsites=cut3$numNA, dir=cut3[,3],mean.pm.diff=cut3$mean.diff, dmr.log.pval=cut3$pval, start.id=cut3[,1],end.id=cut3[,2]))
for(i in 2:ncol(dmrtab)){dmrtab[,i]=as.numeric(as.character(dmrtab[,i]))}
return(dmrtab)
}

##Subroutines of dmrID
runm=caTools::runmean

DMRk.log=function(x,k){
x[is.na(x)==TRUE]=1
pval.sum=runm(log(x),k=k,endrule="NA",align="left")*k
log.pval.seq=sapply(pval.sum,log.punif.prod,n=k)
return(cbind(pval.sum,log.pval.seq))
}

PMsign=function(x,k)
{
pm.sign=k*(runm(x,k=k,endrule="NA",align="left"))
sign.dir=ifelse(abs(pm.sign)==k,sign(pm.sign),0)
return(sign.dir)
}

calc.dens=function(x.pos,k)
{
len=diff(x.pos,lag=(k-1))
dens=k/len
dens=append(dens,rep(NA,k-1))
}

calc.logp.dmr=function(vec,modmat,ptype="adj"){
if(ptype=="orig"){pvec=modmat[vec[1]:vec[2],"pval"]}
else{pvec=modmat[vec[1]:vec[2],"p.adj"]}
p.logsum=sum(log(pvec),na.rm=TRUE)
log.p.dmr=log.punif.prod(p.logsum,length(pvec))
}

log.punif.prod=function(x,n){
k=c(0:(n-1))
kterms=(-1)^{n-1-k}*factorial(n-1)/factorial(k)*(sum(x))^(k)
log.prob=sum(x)+log((-1)^(n-1)/factorial(n-1)*sum(kterms))
return(log.prob)}

calc.mean.diff=function(vec,modmat){
pm.vec=modmat[vec[1]:vec[2],"pm.diff"]
mean.diff=mean(pm.vec,na.rm=TRUE)
}

sum.na=function(vec,modmat,ptype="adj"){
if(ptype=="orig"){pvec=modmat[vec[1]:vec[2],"pval"]}
else{pvec=modmat[vec[1]:vec[2],"p.adj"]}
sum(is.na(pvec))
}
























