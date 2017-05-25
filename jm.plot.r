## Color palettes
# Better version of rainbow() color palette for real data
RGBcolorpal1=hsv(1/(1+exp(seq(-3.8,3.8,length=1000))),1,1)[300:1000]

# Original VGA 16 colors
VGA16colors=c("#000000","#000080","#008000","#008080","#800000","#800080","#808000","#C0C0C0","#808080","#0000FF","#00FF00","#00FFFF","#FF0000","#FF00FF","#FFFF00","#FFFFFF")

# Default hue color palette for ggplot2
ggplotColourHuePallette = function(n=10, h=c(0, 360)+15, c=100, l=65) {
    if((diff(h)%%360) < 1) h[2]=h[2]-360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = c, l = l)
}

# Draw color gradient of x-values
colored_column = function(x, mark=rep(FALSE,length(x)), minx=min(x), maxx=max(x), colmin="blue", colmax="yellow", res=100) {
    ramp = colorRamp(c(colmin, colmax))
    image(x=1, y=1:(length(x)), z=matrix(x,ncol=length(x),nrow=1), axes=F, xlab="", ylab="", col=rgb(ramp(seq(0,1,length=res)),max=255), zlim=c(minx, maxx))
    for(i in 1:length(x)) {
        if(mark[i]) {
            rect(0.6, i-0.5, 1.4, i+0.5, lwd=2)
        }
    }
}

# Scatter plot with marginal histograms
plotxy_marginhists = function(x, y, margins=c(5,5), border=1, h.axis=TRUE, h1.axis=h.axis, h2.axis=h.axis, h.b=20, h1.b=h.b, h2.b=h.b, h.col=8, h1.col=h.col, h2.col=h.col, log="", ...) {
    lmat = matrix(c(2,4,1,3), ncol=2)
    layout(lmat, widths=c(4/5, 1/5), heights=c(1/5, 4/5))
    
    par(mar=c(0,0,0,0))
    plot(1:1,axes=F,type="n",xlab="",ylab="")
    
    df = data.frame(x,y)
    if(grepl("x", log)) { df = df[df[,1]>0,] }
    if(grepl("y", log)) { df = df[df[,2]>0,] }
    x = df[,1]
    y = df[,2]
    
    hx = hist(x, plot=F, breaks=h1.b)
    hy = hist(y, plot=F, breaks=h2.b)
    if(grepl("x", log)) { hx = hist(log(x), plot=F, breaks=h1.b) }
    if(grepl("y", log)) { hy = hist(log(y), plot=F, breaks=h2.b) }
    
    par(mar=c(0,margins[2],1,1))
    barplot(hx$counts, axes=F, ylim=c(0, max(hx$counts)), space=0, col=h1.col, border=border)
    if(h1.axis) { axis(2) }
    
    par(mar=c(margins[1],0,1,1))
    barplot(hy$counts, axes=F, xlim=c(0, max(hy$counts)), space=0, col=h2.col, border=border, horiz=TRUE)
    if(h2.axis) { axis(1) }
    
    par(mar=c(margins[1],margins[2],1,1))
    plot(x, y, log=log, ...)
}

# MA-plot
plot_MA=function(x, y, sig=rep(FALSE,length(x)), ...) {
    x[x<1]=NA
    y[y<1]=NA
    dat=na.omit(cbind(x,y))
    M=log(dat[,1], 2)-log(dat[,2], 2)
    A=0.5*(log(dat[,1],2) + log(dat[,2], 2))
    sigclr=rep(NA,length(x))
    sigclr[sig==TRUE & M>0]=3
    sigclr[sig==TRUE & M<0]=2
    plot(A,M,cex=0.5,pch=1,col=1, ...)
    points(A,M,cex=0.5,pch=20,col=sigclr)
    abline(h=0,col=4)
}

# Density histogram plot
# x2 is optional weights
# multiple lines in same plot with newplot=F
fuzzyhist=function(x,x2=rep(1,length(x)),d=1,n=100,newplot=TRUE,ylab="",...) {
    y=rep(0,length=n)
    z=seq(min(x,na.rm=T),max(x,na.rm=T),length=n)
    rng=max(x,na.rm=T)-min(x,na.rm=T)
    for(i in 1:n) { y[i]=sum((exp(-abs(z[i]-x)*d))*x2,na.rm=T) }
    bs=n/(rng*sum(y,na.rm=T))
    if(newplot) {
        plot(z,y*bs,type="l",ylab=ylab,...)
    } else {
        lines(z,y*bs,...)
    }
}

# Multi-dimensional scaling plot
plotMDS=function(x, cex=0.5, ...) {
    mds=cmdscale(dist(t(x)),k=2)
    plot(mds[,1],mds[,2],xlab="Dimension 1",ylab="Dimension 2",type="n", ...)
    text(mds,names(x), cex=cex)
}

# PCA plot
plotPCA=function(x, cex=0.5, ...) {
    pcafit=prcomp(t(x))
    v1 = round(summary(pcafit)$importance[2,1],2)*100
    v2 = round(summary(pcafit)$importance[2,2],2)*100
    vc = v1 + v2
    plot(pcafit$x[,1], pcafit$x[,2], xlab=paste("PC1 (",v1,"%)",sep=""), ylab=paste("PC2 (",v2,"%)",sep=""), main=paste("PCA ",title," (",vc,"%)",sep=""), type="n", ...)
    text(pcafit$x[,1], pcafit$x[,2], labels=names(x), cex=cex)
}

# Rarefaction curve plot
rarefaction_curve=function(counts, d=100, iter=10, at=seq(1,sum(counts),length=d)) {
    labels = 1:(length(counts))
    probs = counts / sum(counts)
    rmat = matrix(NA, nrow=iter, ncol=length(at))
    for(i in 1:iter) {
        for(j in 1:length(at)) {
            if(ceiling(at[j]) > sum(counts)+1) {
                rmat[i,j] = NA
            } else {
                s = sample(labels, size=ceiling(at[j]), replace=T, prob=probs)
                rmat[i,j] = length(levels(factor(s)))
            }
        }
    }
    rvec = apply(rmat, 2, mean)
    rvec
}

# 2D Kernal density plot
drawAbundHeatmap = function(fin, fout) {
    res=500
    x=read.table(file=fin, header=F)
    jpeg(filename=fout, width=res, height=res, quality=100)
    par(mar=c(0,0,0,0))
    kd=kde2d(x[,1],30000-x[,2],n=res,lims=c(0,30000,0,30000),h=c(2000,2000))
    image(kd$x,kd$y,sqrt(kd$z),col=RGBcolorpal1[150:700])
    dev.off()
}

plotCrossCorr = function(x1, x2, title=NULL, method="kendall") {
    require(ggplot2)
    df = NULL
    for(i in 1:ncol(x1)) {
        for(j in 1:ncol(x2)) {
            ijrow = c(colnames(x1)[i], colnames(x2)[j], cor(x1[,i], x2[,j], method=method), cor.test(x1[,i], x2[,j], method="kendall")$p.value, 1, NA)
            if(is.null(df)){
                df = ijrow
            } else {
                df = rbind(df, ijrow)
            }  
        }
    }
    df = data.frame(row.names=NULL, df)
    colnames(df) = c("Env","Taxa","Correlation","Pvalue","AdjPvalue","Significance")
    df$Pvalue = as.numeric(as.character(df$Pvalue))
    df$AdjPvalue = as.numeric(as.character(df$AdjPvalue))
    df$Correlation = as.numeric(as.character(df$Correlation))
    df$AdjPvalue = p.adjust(df$Pvalue, method="BH")
    df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
    
    p = ggplot(aes(x=Env, y=Taxa, fill=Correlation), data=df) +
        geom_tile() + scale_fill_gradient2(limits=c(-1,1), low="#2C7BB6", mid="white", high="#D7191C") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
        geom_text(aes(label=Significance), color="black", size=3) + 
        labs(y=NULL, x=NULL, fill=method)
    if(!is.null(title)) {
        p = p + ggtitle(title)
    }
    print(p)
}

## The following plots assume: df, kfit, pcafit, apres
#
# read in the data into df, where the first column is the labels, choose a number of clusters k (here k=4)
# kfit = kmeans(df[,-1], 4)
# pcafit = prcomp(df[,-1])
# plotPcaKmeans(pcafit, kfit, df[,1], "PCA, kmeans(4)")
#
# library(apcluster)
# apres = apcluster(negDistMat(r=2), df[,-1], details=TRUE, q=0)
# plotPcaAPcluster(pcafit, apres, df[,1])

plotPcaAPcluster = function(pcafit, apres, lab=NA, title="") {
    v1 = round(summary(pcafit)$importance[2,1],2)*100
    v2 = round(summary(pcafit)$importance[2,2],2)*100
    vc = v1 + v2
    plot(apres, cbind(pcafit$x[,1],pcafit$x[,2]), xlab=paste("PC1 (",v1,"%)",sep=""), ylab=paste("PC2 (",v2,"%)",sep=""), main=paste("PCA ",title," (",vc,"%)",sep=""))
    text(pcafit$x[,1], pcafit$x[,2], labels=lab, cex=0.5, col=1)
}

plotPcaKmeans = function(pcafit, kfit, lab, title="") {
    clrpal = c(2,3,4,6)
    v1 = round(summary(pcafit)$importance[2,1],2)*100
    v2 = round(summary(pcafit)$importance[2,2],2)*100
    vc = v1 + v2
    plot(pcafit$x[,1],pcafit$x[,2],xlab=paste("PC1 (",v1,"%)",sep=""),ylab=paste("PC2 (",v2,"%)",sep=""),main=paste("PCA ",title," (",vc,"%)",sep=""),pch=5,col=8)
    text(pcafit$x[,1],pcafit$x[,2],labels=lab,cex=0.5,col=clrpal[kfit$cluster])
}

plotPCAfit = function(pcafit, lab, title="", col=8, pch=5, ...) {
    v1 = round(summary(pcafit)$importance[2,1],2)*100
    v2 = round(summary(pcafit)$importance[2,2],2)*100
    vc = v1 + v2
    plot(pcafit$x[,1], pcafit$x[,2], xlab=paste("PC1 (",v1,"%)",sep=""), ylab=paste("PC2 (",v2,"%)",sep=""), main=paste("PCA ",title," (",vc,"%)",sep=""), pch=pch, col=col)
    text(pcafit$x[,1], pcafit$x[,2], labels=lab, cex=0.5, col=col)
}

plotPcaBiplot = function(pcafit, lab, vec.lines=T, title="") {
    labvar=rownames(pcafit$r)
    v1 = round(summary(pcafit)$importance[2,1],2)*100
    v2 = round(summary(pcafit)$importance[2,2],2)*100
    vc = v1 + v2
    plot(pcafit$x[,1],pcafit$x[,2],xlab=paste("PC1 (",v1,"%)",sep=""),ylab=paste("PC2 (",v2,"%)",sep=""),main=paste("PCA ",title," (",vc,"%)",sep=""),type="n",xlim=c(min(pcafit$x[,1]),max(pcafit$x[,1]))*1.1,ylim=c(min(pcafit$x[,2]),max(pcafit$x[,2]))*1.1)
    text(pcafit$x[,1],pcafit$x[,2],labels=lab,cex=0.5,col=1)
    par(new=T)
    plot(pcafit$r[,1],pcafit$r[,2],type="n",axes=F,xlab="",ylab="",xlim=c(min(pcafit$r[,1]),max(pcafit$r[,1]))*1.4,ylim=c(min(pcafit$r[,2]),max(pcafit$r[,2]))*1.4)
    axis(3)
    axis(4)
    if(vec.lines) {
        for(i in 1:(length(pcafit$r[,1]))) {
            lines(c(0,pcafit$r[i,1]), c(0,pcafit$r[i,2]), col=8)
        }
    }
    text(pcafit$r[,1]*1.1,pcafit$r[,2]*1.1,labels=labvar,cex=0.5,col=2)
}

##Draws a log10/log10 plot of data with proper grid lines
#Minexp and Maxexp set plot limits as round exponents of power 10, default to range of x,y data given
#Can set any other plotting parameters: (col, pch, xlab, ylab, main, etc.), but defaults to pch=20 filled circles
#Set gridlinecol to primary and secondary colors, defaults to 1=black, 8=gray
#Set gridlinetype to primary and secondary line types, defaults to 1=solid, 1=solid
#    (You might like gridlinetype=c(1,2) or c(1,3) for dashed or dotted secondary lines

loglog10plot=function(x,y,minexp=NA,maxexp=NA,gridlinecol=c(1,8),gridlinetype=c(1,1),pch=20,...) {
    if(is.na(minexp)) { minexp=round(log(min(min(x),min(y)),10))-1 }
    if(is.na(maxexp)) { maxexp=round(log(max(max(x),max(y)),10))+1 }
    logx=log(x,10)
    logy=log(y,10)
    
    if(maxexp<minexp+1) { maxexp=minexp+1 }
    axispos=10^minexp
    axislab=10^minexp
    axistick=1
    tickcol=gridlinecol[1]
    ticktype=gridlinetype[1]
    
    for(i in minexp:(maxexp-1)) {
        axispos=c(axispos,seq(10^i,10^(i+1),length=10)[2:10])
        axislab=c(axislab,rep(NA,8),10^(i+1))
        axistick=c(axistick,rep(2,8),1)
        tickcol=c(tickcol,rep(gridlinecol[2],8),gridlinecol[1])
        ticktype=c(ticktype,rep(gridlinetype[2],8),gridlinetype[1])
    }
    
    plot(logx,logy,axes=F,type="n",xlim=c(minexp,maxexp),ylim=c(minexp,maxexp),...)
    axis(1,at=log(axispos,10),labels=axislab)
    axis(2,at=log(axispos,10),labels=axislab)
    
    # Secondary grid lines
    for(i in 1:(length(axispos))) {
        if(axistick[i]==2) {
            lines(c(minexp,maxexp),c(log(axispos[i],10),log(axispos[i],10)),col=tickcol[i],lty=ticktype[i])
            lines(c(log(axispos[i],10),log(axispos[i],10)),c(minexp,maxexp),col=tickcol[i],lty=ticktype[i])
        }
    }
    # Primary grid lines
    for(i in 1:(length(axispos))) {
        if(axistick[i]==1) {
            lines(c(minexp,maxexp),c(log(axispos[i],10),log(axispos[i],10)),col=tickcol[i],lty=ticktype[i])
            lines(c(log(axispos[i],10),log(axispos[i],10)),c(minexp,maxexp),col=tickcol[i],lty=ticktype[i])
        }
    }
    # Draw points on top of grid
    points(logx,logy,pch=pch,...)
}

#Example:
# (data must be positive to be logged)
#
# x=5^rnorm(1000)
# y=x*2^rnorm(1000)
#
# loglog10plot(x,y,xlab="Something",ylab="Other Stuff",main="Other Stuff vs. Something",col=2,pch=1)
#
# loglog10plot(x,y,minexp=-3,maxexp=3,gridlinecol=c(1,8),gridlinetype=c(1,3))
#


## Arranges the rows in metrix mat into a 2-dimensional space by correlation of vectors
vectorSort2D = function(mat, xrange=c(floor(sqrt(nrow(mat)))-1, ceiling(sqrt(nrow(mat)))+1)) {
    require(som)
    cortot.best = NA
    dfasgn.best = NA
    for(xdim in seq(xrange[1], xrange[2], by=1)) {
        msize = nrow(mat)
        ydim = ceiling(msize / xdim)
        xsize = xdim * ydim
        
        xsom = som(mat, xdim=xdim, ydim=ydim, topol="rect")
        
        cmat=matrix(nrow=msize*xsize,ncol=3)
        v=0
        for(i in 1:msize) { for(j in 1:xsize) { v=v+1; cmat[v,1]=i; cmat[v,2]=j; cmat[v,3]=cor(xsom$data[i,], xsom$code[j,]) }}
        
        cmat.s = cmat[order(cmat[,3], decreasing=T),]
        masgn = rep(0, msize)
        xasgn = rep(0, xsize)
        cortot = 0
        for(i in 1:(nrow(cmat.s))) {
            if(xasgn[cmat.s[i,2]] == 0 & masgn[cmat.s[i,1]] == 0) {
                xasgn[cmat.s[i,2]] = cmat.s[i,1]
                masgn[cmat.s[i,1]] = cmat.s[i,2]
                cortot = cortot + cmat.s[i,3]
            }
        }
        if(is.na(cortot.best) | cortot > cortot.best) {
            cortot.best = cortot
            dfasgn.best = data.frame(xasgn, xsom$code.sum$x, xsom$code.sum$y)
            names(dfasgn.best) = c("row","x","y")
        }
    }
    dfasgn.best
}

vectorSort2D.plot = function(mat, dfasgn) {
    dfasgn.s = dfasgn[order(dfasgn$y, dfasgn$x),]
    xdim = max(dfasgn$x)+1
    ydim = max(dfasgn$y)+1
    par(mfrow=c(ydim,xdim), mar=c(0.2,0.2,2.0,0.2))
    for(i in 1:(nrow(dfasgn.s))) {
        if(dfasgn.s$row[i] > 0) {
            plot(mat[dfasgn.s$row[i],], type="n", axes=F, xlab="", ylab="", main=rownames(mat)[dfasgn.s$row[i]])
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="#F0F0F0", border=NA)
            lines(mat[dfasgn.s$row[i],])
        } else {
            plot(1,1,type="n",axes=F,xlab="",ylab="")
        }
    }
}
