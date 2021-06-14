ls()
getwd()
list.files()

GSE4183_data<-read.csv("GSE4183_data.csv")
GSE4183_meta_data<-read.csv("GSE4183_meta_data.csv")

GSE4183_data0<-as.matrix(GSE4183_data[,-1])
rownames(GSE4183_data0)<-as.character(GSE4183_data[,1])

GSE4183_meta_data<-as.matrix(GSE4183_meta_data)

GSE4183_meta_data[17,agrep("cancer",GSE4183_meta_data[17,])]
GSE4183_meta_data[17,agrep("health",GSE4183_meta_data[17,])]

healthy_id<-agrep("health",GSE4183_meta_data[17,])
IBD_id<-agrep("inflamm",GSE4183_meta_data[17,])
Adenoma_id<-agrep("adenoma",GSE4183_meta_data[17,])
cancer_id<-agrep("cancer",GSE4183_meta_data[17,])


data0<-GSE4183_data0[1:500,c(cancer_id,healthy_id)]
colnames(data0)<-c(rep("C",length(cancer_id)),rep("H",length(healthy_id)))
heatmap(data0)

################################################

a<-matrix(c(1:20),4,5)

sum(a)
apply(a,1,sum)
apply(a,2,sum)
mean(a)

data01<-log(data0)
min(data0)
max(data0)
hist(apply(data0,1,mean),breaks=30)
hist(apply(data01,1,mean),breaks=30)

hmcols2 <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#91BFDB", "#4575B4")))(4)
plot(1:4,pch=16,cex=1.5,col=hmcols2 )

sample_info<-cbind(colnames(GSE4183_data0),colnames(GSE4183_data0),colnames(GSE4183_data0))
sample_info[healthy_id,2]<-"healthy"
sample_info[IBD_id,2]<-"IBD"
sample_info[Adenoma_id,2]<-"Adenoma"
sample_info[cancer_id,2]<-"cancer"
sample_info[healthy_id,3]<-hmcols2[1]
sample_info[IBD_id,3]<-hmcols2[2]
sample_info[Adenoma_id,3]<-hmcols2[3]
sample_info[cancer_id,3]<-hmcols2[4]

sample_info<-as.data.frame(sample_info)
colnames(sample_info)<-c("ID","Condition","Color")

library(limma)

a<-plotMDS(GSE4183_data0,col=sample_info$Color,cex=1,ndim=2)

Dimension reduction (Principal componental analysis)

Dendrogram

GSE4183_data01<-GSE4183_data0
colnames(GSE4183_data01)<-sample_info[,2]

plot(hclust(dist(t(GSE4183_data01))))


#############
#Random sampling and random generation of classes
########

xa<-rnorm(100,0,1)
xb<-xa+rnorm(100,0,0.3)
plot(xa,xb)

score<-xa+xb
c1<-which(score<mean(score))
c2<-which(score>=mean(score))
points(xa[c1],xb[c1],col="red")
points(xa[c2],xb[c2],col="blue")

range(score)

par(mfcol=c(2,5))
for(i in 1:10)
{
temp_p<-(score+4.3)/10
rp<-runif(100,0,1)

class<-(temp_p<rp)*1
c1<-which(class==0)
c2<-which(class==1)
plot(xa,xb)
points(xa[c1],xb[c1],col="red",pch=16)
points(xa[c2],xb[c2],col="blue",pch=16)
}


#############
#histogram
#######

mean(score[c1])
mean(score[c2])

p1<-hist(score[c1])
p2<-hist(score[c2])

plot( p1, col=rgb(1,0,0,1/2), xlim=c(-5,6))  # first histogram
plot( p2, col=rgb(0,0,1,1/2), xlim=c(-5,6), add=T)  # second


#############
#t test
###############
t.test(score[c1],score[c2])

#p-value = 2.669e-06

#p-value the probabilty of the observation when (H0:null hypothesis) is true
#H0 : boys' height = girl's height

#######################
#permutation test
################
observed_boys_height<-score[c1]
observed_girls_height<-score[c2]
hd0<-mean(observed_boys_height)-mean(observed_girls_height)


pooled_height<-c(observed_boys_height,observed_girls_height)

hdr_stat<-c()
for(i in 1:10000)
{
	boys_r<-sample(pooled_height,50)
	girls_r<-sample(pooled_height,50)
	hdr<-mean(boys_r)-mean(girls_r)
	hdr_stat<-c(hdr_stat,hdr)
}
hist(hdr_stat,col="lightblue",xlim=c(-2,2))
abline(v=hd0,col=2,lwd=2)


sum(hdr_stat>1)/10000












