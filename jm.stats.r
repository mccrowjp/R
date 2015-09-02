require(MASS)

# Standard error for means
stderr = function(x){ sqrt(var(x,na.rm=TRUE)/length(na.omit(x))) }

# Geometric mean
gmean=function(x) { exp(mean(log(x))) }

# Power-law distribution functions
dpareto=function(x, s=1, l=1) s * l^s / x^(s + 1);
ppareto=function(x, s=1, l=1) (x >= l) * (1 - (l / x)^s);
qpareto=function(u, s=1, l=1) l/(1 - u)^(1 / s);
rpareto=function(n, s=1, l=1) qpareto(runif(n), s, l);

# Chi-square goodness of fit, d.f.=(discrete bins) - ((num parameters used in fit)+1)
goodfit.chisq.test=function(x,dfunc,params) {
    n=length(x)
    p=length(params)+1
    xvals=floor(min(x)):(max(x))
    ofset=min(x)-1
    counts=rep(0,length(xvals))
    for(i in 1:n) { counts[x[i]-ofset]=counts[x[i]-ofset]+1 }
    ob=counts[counts>5]
    df=length(ob)-p
    cl=call(dfunc, xvals[counts>5], params)
    if(is.call(cl)) {
        ex=eval(cl)
        ex=round(ex*sum(ob),0)
        ctmat=as.matrix(cbind(ob,ex))
        statistic=chisq.test(ctmat, correct=F)$statistic
        p.value=1-pchisq(statistic,df)
        list(statistic=statistic, df=df, p.value=p.value)
    } else {
        warning("unrecognized function")
    }
}

# List of distributions fit to empirical data x, overlaid on plot
normhist=function(x,col=8,normcol=4,meancol=4,sdcol=5,border=F,xlab=deparse(substitute(x)),...) {
    truehist(x,col=col,border=border,xlab=xlab,...)
    lim1=par("usr")[1]
    lim2=par("usr")[2]
    nm=mean(x)
    nsd=sd(x)
    nx=seq(lim1,lim2,length=1000)
    ny=dnorm(nx,mean=nm,sd=nsd)
    if(!is.na(normcol)){ lines(nx,ny,col=normcol) }
    if(!is.na(meancol)){ lines(c(nm,nm),c(0,dnorm(nm,mean=nm,sd=nsd)),col=meancol) }
    if(!is.na(sdcol)){ lines(c(nm+nsd,nm+nsd),c(0,dnorm(nm+nsd,mean=nm,sd=nsd)),col=sdcol) }
    if(!is.na(sdcol)){ lines(c(nm-nsd,nm-nsd),c(0,dnorm(nm-nsd,mean=nm,sd=nsd)),col=sdcol) }
    p=ks.test(x,"pnorm",nm,nsd)$p
    legend("topright",inset=0.05,legend=c(paste("Normal( ",round(nm,5),", ",round(nsd,5)," )"), "KS-test", paste("p =",signif(p,5))),bty="n")
}

lognormhist=function(x,col=8,lnormcol=4,meancol=4,sdcol=5,border=F,xlab=deparse(substitute(x)),...) {
    truehist(x,col=col,border=border,xlab=xlab,...)
    lim1=par("usr")[1]
    lim2=par("usr")[2]
    nm=gmean(x)
    nsd=sd(x)
    lnm=mean(log(x))
    lnsd=sd(log(x))
    nx=seq(lim1,lim2,length=1000)
    ny=dlnorm(nx,meanlog=lnm,sdlog=lnsd)
    if(!is.na(lnormcol)){ lines(nx,ny,col=lnormcol) }
    if(!is.na(meancol)){ lines(c(nm,nm),c(0,dlnorm(nm,meanlog=lnm,sdlog=lnsd)),col=meancol) }
    if(!is.na(sdcol)){ lines(c(nm+nsd,nm+nsd),c(0,dlnorm(nm+nsd,meanlog=lnm,sdlog=lnsd)),col=sdcol) }
    p=ks.test(x,"plnorm",lnm,lnsd)$p
    legend("topright",inset=0.05,legend=c(paste("Lognormal( ",round(lnm,5),", ",round(lnsd,5)," )"), "KS-test", paste("p =",signif(p,5))),bty="n")
}

gammahist=function(x,shape,legtxt="",col=8,gamcol=4,meancol=4,sdcol=NA,border=F,xlab=deparse(substitute(x)),...) {
    truehist(x,col=col,border=border,xlab=xlab,...)
    lim1=0
    lim2=par("usr")[2]
    nm=mean(x)
    nsd=sd(x)
    nx=seq(lim1,lim2,length=1000)
    ny=dgamma(nx,shape=shape,rate=1/nm)
    if(!is.na(gamcol)){ lines(nx,ny,col=gamcol) }
    if(!is.na(meancol)){ lines(c(nm,nm),c(0,dgamma(nm,shape=shape,rate=1/nm)),col=meancol) }
    if(!is.na(sdcol)){ lines(c(nm+nsd,nm+nsd),c(0,dgamma(nm+nsd,shape=shape,rate=1/nm)),col=sdcol) }
    p=ks.test(x,"pgamma",shape,1/nm)$p
    leg=c("KS-test",paste("p =",signif(p,5)))
    if(nchar(legtxt)>0) {
        leg=c(legtxt, leg)
    } else {
        leg=c(paste("Gamma( ",round(shape,5),", 1/",round(nm,5)," )",sep=""), leg)
    }
    legend("topright",inset=0.05,legend=leg,bty="n")
}

exphist=function(x,col=8,expcol=4,meancol=4,sdcol=NA,border=F,xlab="",...) {
    gammahist(x,shape=1,paste("Exponential( 1, 1/",round(mean(x),5)," )",sep=""), col=col,gamcol=expcol,meancol=meancol,sdcol=sdcol,border=border,xlab=xlab,...)
}

discretehist=function(x,dfunc,param,legtxt="",col=8,denscol=4,meancol=4,sdcol=NA,border=F,xlab=deparse(substitute(x)),...) {
    truehist(x,h=1,col=col,border=border,xlab=xlab,...)
    lim1=0
    lim2=par("usr")[2]
    nm=mean(x)
    nsd=sd(x)
    nx=(lim1):(lim2)
    cl=call(dfunc, nx, param)
    cl2=call(dfunc, floor(nm), param)
    cl3=call(dfunc, floor(nm+nsd), param)
    if(is.call(cl)) {
        ny=eval(cl)
        if(!is.na(denscol)){ lines(nx,ny,col=denscol,type="s") }
        if(!is.na(meancol)){ lines(c(nm,nm),c(0,eval(cl2)),col=meancol) }
        if(!is.na(sdcol)){ lines(c(nm+nsd,nm+nsd),c(0,eval(cl3)),col=sdcol) }
        p=goodfit.chisq.test(x,dfunc,param)$p.value
        leg=c(expression(paste("Pearson ",chi^2,"-test")),paste("p =",signif(p,5)))
        if(nchar(legtxt)>0){ leg=c(legtxt,leg) }
        legend("topright",inset=0.05,legend=leg,bty="n")
    } else {
        warning("not a recognized function")
    }
}

geomhist=function(x,col=8,geomcol=4,meancol=4,sdcol=NA,border=F,xlab=deparse(substitute(x)),...) {
    param=1/mean(x)
    discretehist(x,"dgeom",param,paste("Geometric( 1/",round(mean(x),5)," )",sep=""),col,geomcol,meancol,sdcol,border,xlab,...)
}

poishist=function(x,col=8,poiscol=4,meancol=4,sdcol=NA,border=F,xlab=deparse(substitute(x)),...) {
    param=mean(x)
    discretehist(x,"dpois",param,paste("Poisson( ",round(param,5)," )",sep=""),col,poiscol,meancol,sdcol,border,xlab,...)
}

paretohist=function(x,col=8,pcol=4,meancol=4,sdcol=NA,border=F,xlab=deparse(substitute(x)),...) {
    truehist(x,col=col,border=border,xlab=xlab,...)
    lim1=par("usr")[1]
    lim2=par("usr")[2]
    nm=mean(x)
    nsd=sd(x)
    location=min(x)
    nx=seq(lim1,lim2,length=1000)
    
    # estimate shape parameter, by linear regression of head of distribution on log-log
    logx=log(x)
    lxbins=seq(min(logx),max(logx),length=50)
    lxcounts=rep(0,10)
    for(i in 1:10) { lxcounts[i]=length(logx[logx>=lxbins[i] & logx<lxbins[i+1]])}
    lxmids=(lxbins[2:11]+lxbins[1:10])/2
    shape=-lm(log(lxcounts)~lxmids)$coefficients[2]
    
    ny=dpareto(nx,shape,location)
    if(!is.na(pcol)){ lines(nx,ny,col=pcol) }
    if(!is.na(meancol)){ lines(c(nm,nm),c(0,dpareto(nm,shape,location)),col=meancol) }
    if(!is.na(sdcol)){ lines(c(nm+nsd,nm+nsd),c(0,dpareto(nm+nsd,shape,location)),col=sdcol) }
    p=ks.test(x,"ppareto",shape,location)$p
    legend("topright",inset=0.05,legend=c(paste("Pareto( ",round(shape,5),", ",round(location,5)," )",sep=""),"KS-test",paste("p =",signif(p,5))),bty="n")
}

powerlawhist=paretohist

# Plot 6 distributions for comparison of fit
allhist=function(x,...) {
    par(mfrow=c(2,3))
    options(warn=-1)
    try(normhist(x,nbins=50,...),T)
    try(lognormhist(x,nbins=50,...),T)
    try(exphist(x,nbins=50,...),T)
    try(geomhist(x,nbins=50,...),T)
    try(poishist(x,nbins=50,...),T)
    try(paretohist(x,nbins=50,...),T)
    options(warn=1)
}

