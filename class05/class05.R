#line plot
weight<- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)
plot(weight$Age, weight$Weight, 
     type="o",pch=15, col="blue", cex=1.5, lwd=2, ylim=c(2,10),
     xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight with Age")

#bar plot
mouse<- read.delim("bimm143_05_rstats/feature_counts.txt", header=TRUE)
par(mar=c(3,11,3,2))
barplot(mouse$Count,horiz=TRUE,ylab="", main="Features in mouse GRCm38 Genome", 
        names.arg = mouse$Feature, las=1, xlim=c(0,80000))

#histogram
x<- c(rnorm(10000),rnorm(10000)+4)
hist(x,breaks=100)

#barplot
gender<- read.delim("bimm143_05_rstats/male_female_counts.txt", header=TRUE)
barplot(gender$Count,col=c("blue","red"), names.arg=gender$Sample, las=2, ylab="Counts")

#plot by color
exp<- read.delim("bimm143_05_rstats/up_down_expression.txt", header=TRUE)
palette(c("blue","grey","red"))
plot(exp$Condition1,exp$Condition2,col=exp$State)

#color density
meth<-read.delim("bimm143_05_rstats/expression_methylation.txt", header=TRUE)
densC<- densCols(meth$gene.meth,meth$expression)
plot(meth$gene.meth,meth$expression,col=densC, pch=20)

#exclude 0
exclude <- meth$expression>0
exdensC <- densCols(meth$gene.meth[exclude],meth$expression[exclude])
plot(meth$gene.meth[exclude],meth$expression[exclude], col=exdensC, pch=20)

#color ramp

c_exdensC<- densCols(meth$gene.meth[exclude],meth$expression[exclude], 
                     colramp= colorRampPalette(c("blue","green","red","yellow")))
plot(meth$gene.meth[exclude],meth$expression[exclude], col=c_exdensC,pch=20)


