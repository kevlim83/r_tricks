#generate some nice figures

#1 bar plot with error bars
#data.mean
#data.sd
#data.N
library(plyr)
myrange<-max(data.mean+1.96*data.sd/sqrt(data.N))
b<-barplot(data.mean,las=1,ylim=c(0,myrange)
l_ply(seq_along(b),function(x) arrows(x0=b[x],y0=data.mean[x],x1=b[x],y1=data.mean[x]+1.96*data.sd[x]/sqrt(data.N[x]), code=2, length=0.1, angle=90, lwd=1.5))
