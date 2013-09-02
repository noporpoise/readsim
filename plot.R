args <- commandArgs(trailingOnly = TRUE)
print(args)

x=read.csv(file=args[1],header=F)
x=as.numeric(x)
pdf(file=args[2],width=8,height=8)
plot(1:length(x),x,type='l',ylim=c(0,0.05),
     xlab='ReadBase',ylab='ErrorRate',main='Seq Error Rate')
dev.off()

meanPercent=mean(x)*100
minPercent=min(x)*100
maxPercent=max(x)*100
print(meanPercent)
print(minPercent)
print(maxPercent)
