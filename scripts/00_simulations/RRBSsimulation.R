library(copula)
library(HiddenMarkov)
library(RandomFields)

##sim.beg.vg is the main function to generate a simulation using default parameters.
##The following line will generate 6 samples of 1000 sites divided into two groups with an overall DM rate of 10% and DMRs spanning 250bp with 
##average PMD = 0.10.  When DMRs are included, output is a list of length 2 with the DMR locations included in the second entry.
##sim.n3.d10.l250=sim.bed.vg(nsites=1000,nsamp=6,dm.len=250,prop.diff=0.10,pm.diff=0.10)

nsites = 20
nsamp = 50
dm.len = 10
prop.diff = 0
pm.diff = 0.1
sim.bed.vg=function(nsites=1000,nsamp=6,dm.len=2000,prop.diff=0,pm.diff=0.25)
{
    lagvec=lag.sim.hmm(len=(nsites-1))
    chrpos=c(1,cumsum(lagvec))
    cat("Done generating positions \n")
    
    
    totlen=max(chrpos)
    nscores=length(unique(floor(chrpos/50)))
    scoremat=score.sim.mix(len=nscores,nsamp=nsamp,cor=rep(0.8,choose(nsamp,2)))
    scorepos=as.numeric(as.factor(rank((floor(chrpos/50)),ties.method="average")))
    scoretab=as.data.frame(matrix(0,length(chrpos),nsamp))
    for(i in 1:nscores){
    if(i%%1000==0){cat("Working on score ", i, "\n")}
    for(j in 1:nsamp){
    scoretab[scorepos==i,j]=scoremat[i,j]}}
    names(scoretab)[1:nsamp]=paste("SC",1:nsamp,sep=".")
    cat("Done generating score matrix \n")
    
    
    score.med=round(apply(scoretab,1,median))
    mprob=sim.corr.pm(samp.pos=chrpos)
    dmsites=rep(0,nsites)
    if(prop.diff>0){
        diff.regs=dm.regs(chrpos, dm.len=dm.len, min.sites=dm.len/50)
        avg.nsites.dm=mean(diff.regs[,3])
        num.dm.regs=nsites*prop.diff/avg.nsites.dm
        if(num.dm.regs<nrow(diff.regs)){
            dm.rows=sample(1:nrow(diff.regs),num.dm.regs)}
        else{dm.rows=1:nrow(diff.regs)
            num.dm.regs=nrow(diff.regs)}
        dmsites=rep(0,nsites)
        diff.regs.samp=diff.regs[dm.rows,,drop=FALSE]
        diffdir=sample(c(1,-1),num.dm.regs,replace=T)
        if (num.dm.regs >= 1){
            for(i in 1:num.dm.regs){
                print(diff.regs.samp)
                mprob.med=median(mprob[diff.regs.samp[i,1]:(diff.regs.samp[i,1]+diff.regs.samp[i,3]-1)])
                if((mprob.med-pm.diff)<0){diffdir[i]= 1}
                if((mprob.med+pm.diff)>1){diffdir[i]= -1}
                dmsites[diff.regs.samp[i,1]:(diff.regs.samp[i,1]+diff.regs.samp[i,3]-1)]=diffdir[i]}
        } else{
            mprob.med=0
        }
    }
    cat("Done generating DM regions \n")
    
    
    mprob.diff=mprob+dmsites*pm.diff
    mprob.diff[mprob.diff <0]=0.01
    mprob.diff[mprob.diff >1]=0.99
    pmtab=as.data.frame(matrix(0,length(chrpos),nsamp))
    for(j in 1:(nsamp/2))
    {pmtab[,j]=round(rbinom(nsites,size=scoretab[,j],prob=mprob)*(100/scoretab[,j]))}
    for(j in (nsamp/2+1):nsamp)
    {pmtab[,j]=round(rbinom(nsites,size=scoretab[,j],prob=mprob.diff)*(100/scoretab[,j]))}
    names(pmtab)[1:nsamp]=paste("PM",1:nsamp,sep=".")
    cat("Done generating PM matrix \n")
    
    
    restab=cbind(chrpos,dmsites,mprob,mprob.diff, scoretab,pmtab)
    if(prop.diff>0){
    return(list(restab,diff.regs.samp[order(diff.regs.samp[,1]),]))}
    else{return(restab)}
}


##HMM parameters for lag.sim.hmm
lag.hmm.pars=list()
lag.hmm.pars$Pi=matrix(c(.85,.45,.15,.55),2,2)
lag.hmm.pars$delta=c(0.75,.25)
lag.hmm.pars$pm$meanlog=c(2.08,6.12)
lag.hmm.pars$pm$sdlog=c(0.90,2.84)

##Subroutine to generate CpG site locations
lag.sim.hmm=function(len=1000,hmm.pars=lag.hmm.pars)
{
dx=dthmm(NULL,Pi=hmm.pars$Pi,delta=hmm.pars$delta,distn="lnorm",pm=hmm.pars$pm)
sim.hmm=simulate(dx,nsim=len)
lagvec=round(sim.hmm$x)
hv=which(log(lagvec)>15)
lagvec[hv]=lagvec[hv]/10
lv=which(lagvec<2)
lagvec[lv]=rep(2,length(lv))
return(lagvec)
}

##Subroutine to generate coverage scores using normal copula
score.sim.mix=function(len=500,param=c(0.34,1.62,1.03,13.08,3.70),nsamp=2,cor=c(0.8))
{
stype=rep(0,len)
ntype1=round(len*param[1])
dtype=sample(1:len,ntype1)
stype[dtype]=1
mv1<- mvdc(normalCopula(param=cor,dim=nsamp, dispstr="un"), margins="gamma",
              paramMargins=list(list(shape=param[2],rate=param[3])) ,marginsIdentical=TRUE)
mv2<- mvdc(normalCopula(param=cor,dim=nsamp, dispstr="un"), margins="gamma",
              paramMargins=list(list(shape=param[4],rate=param[5])) ,marginsIdentical=TRUE)
mvsamp1 <- rMvdc(ntype1, mv1)
mvsamp2 <- rMvdc(len-ntype1,mv2)
mvsamp=matrix(0,len,nsamp)
mvsamp[stype==1,]=round(exp(mvsamp1))
mvsamp[stype==0,]=round(exp(mvsamp2))
return(mvsamp)
}


##Gaussian variogram model with estimated parameters
covmodel=RMplus(RMgauss(var=0.151,scale=1594.341),RMnugget(var=0.014),RMtrend(mean=0.68))

##Subroutine to generate locally correlated PM levels
sim.corr.pm=function(samp.pos,shape1=0.14, shape2=0.18)
{
sim.len=length(samp.pos)
pm.ind=rbeta(sim.len,shape1,shape2)
pm.ind[pm.ind<1e-15]=1e-15
pm.ind[pm.ind>1-(1e-15)]=1-(1e-15)
pm.corr=0
pos.vec=samp.pos
pos.breaks=c(0,which(diff(pos.vec)>5000),length(pos.vec))
for(i in 1:(length(pos.breaks)-1))
{
pm.curr=pm.ind[(pos.breaks[i]+1):(pos.breaks[i+1])]
#cat("pm.curr = \n")
#print(pm.curr)
pos.curr=pos.vec[(pos.breaks[i]+1):(pos.breaks[i+1])]
covmat.curr=as.matrix(RFcovmatrix(x=pos.curr,model=covmodel))
#cat("covmat.curr = \n")
#print(covmat.curr)
corrmat.curr=cov2cor(covmat.curr)
chol.curr=t(chol(corrmat.curr))
pm.trans=1-pnorm(chol.curr%*%qnorm(1-pm.curr))
pm.corr[(pos.breaks[i]+1):(pos.breaks[i+1])]=pm.trans
}
return(pm.corr)
}


##Subroutine to identify potential locations for DMRs
dm.regs=function(x=chrpos,dm.len=2000,min.sites=5){
nsites=length(x)
dmstarts=dmlen=sitenum=0
siteend=-1
for(i in 1:(nsites-1)){
if(i>(siteend+1) & (x[i+1]-x[i]) <dm.len){
curr.int=paste("[",x[i], ",",x[i]+dm.len-1,"]")
curr.sites=which(x %in% interval(curr.int))
curr.len=length(curr.sites)
if(curr.len >= min.sites){
dmstarts=c(dmstarts,x[i])
sitenum=c(sitenum,i)
dmlen=c(dmlen,curr.len)
siteend=i+curr.len-1
}
}
}
return(cbind(sitenum[-1],dmstarts[-1],dmlen[-1]))
}

