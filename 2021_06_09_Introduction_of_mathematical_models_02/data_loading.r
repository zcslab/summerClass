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

cancer_id<-agrep("cancer",GSE4183_meta_data[17,])
healthy_id<-agrep("health",GSE4183_meta_data[17,])


data0<-GSE4183_data0[1:500,c(cancer_id,healthy_id)]
colnames(data0)<-c(rep("C",length(cancer_id)),rep("H",length(healthy_id)))
heatmap(data0)