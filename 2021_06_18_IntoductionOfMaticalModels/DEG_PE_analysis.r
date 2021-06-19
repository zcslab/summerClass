
load("TCGA-COAD.RData")
data_cancer<-log(d.matrix+1)
load("TCGA-COAD_N.RData")
data_normal<-log(d.matrix+1)
#stat0<-wilcox_test_all(data_cancer,data_normal)

load("MsigDB_gene_sets_v6.RData")

cc_genes<-MsigDB_gene_sets_v6[[2]][["KEGG_CELL_CYCLE"]]
cc_genes<-intersect(rownames(data_cancer),cc_genes)

data_cancer_cc<-data_cancer[cc_genes,]
data_normal_cc<-data_normal[cc_genes,]
data_all_cc<-cbind(data_cancer_cc,data_normal_cc)

colname_temp<-c(rep("Cancer",ncol(data_cancer_cc)),rep("N",ncol(data_normal_cc)))
colnames(data_all_cc)<-colname_temp



library(gplots)
colors = c(-200:200)/100

aa<--10:10/5
bb<-rbind(aa,-aa)



h<-heatmap.2(bb,Rowv=F,Colv =F,scale="none",main="",
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(5,5),cexRow=1,cexCol=1)


h<-heatmap.2(data_all_cc[1:10,1:10],Rowv=T,Colv =T,scale="row",main="",
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(5,5),cexRow=1,cexCol=1)




data_all_cc[1:10,1:10]

a<-matrix(,)

my_palette <- colorRampPalette(c("red","white", "blue"))(n =400)

col0<-c(rep("purple",ncol(data_cancer_cc)),rep("yellow",ncol(data_normal_cc)))

h<-heatmap.2(data_all_cc,Rowv=T,Colv =T,scale="row",main="",ColSideColors=col0,
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(5,5),cexRow=1,cexCol=1)


data_all_cc0<-data_all_cc
for(i in 1:nrow(data_all_cc0))
{
	data_all_cc0[i,]<-sample(data_all_cc0[i,])
}
colnames(data_all_cc0)<-colnames(data_all_cc)
for(i in 1:ncol(data_all_cc0))
{
	data_all_cc0[,i]<-sample(data_all_cc0[,i])
}
rownames(data_all_cc0)<-rownames(data_all_cc)


h<-heatmap.2(data_all_cc0,Rowv=T,Colv =T,scale="row",main="",ColSideColors=col0,
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(5,5),cexRow=1,cexCol=1)


plot(hclust(dist(t(data_all_cc))),cex=0.4)
plot(hclust(dist(data_all_cc)),cex=0.4)


h<-hclust(dist(t(data_all_cc)))
clusters<-cutree(h,4)



boxplot(data_cancer[2,],data_normal[2,])



p_all<-c()
mean_diff_all<-c()
for(i in 1:nrow(data_cancer))
{
	p_all<-c(p_all,t.test(data_cancer[i,],data_normal[i,])$p.value)
	mean_diff_all<-c(mean_diff_all,mean(data_cancer[i,])-mean(data_normal[i,]))
}
results<-cbind(p_all,mean_diff_all)
rownames(results)<-rownames(data_cancer)


up_results<-results[which(results[,2]>0),]
up_list<-rownames(up_results[order(up_results[,1])[1:500],])

down_results<-results[which(results[,2]<0),]
down_list<-rownames(down_results[order(down_results[,1])[1:500],])

all_pathways<-MsigDB_gene_sets_v6[[2]]
all_genes<-rownames(results)
p_all<-c()
sign_all<-c()
for(i in 1:length(all_pathways))
{
	stat_table_c<-matrix(0,3,3) #create a 3*3 matrix
	rownames(stat_table_c)<-c("G_pathway","G_n_pathway","Total")
	colnames(stat_table_c)<-c("DEG","nDEG","Total")
	pathway_c<-all_pathways[[i]]
	stat_table_c[3,1]<-length(up_list)
	stat_table_c[1,3]<-length(pathway_c)
	stat_table_c[3,3]<-length(all_genes)
	stat_table_c[3,2]<-stat_table_c[3,3]-stat_table_c[3,1]
	stat_table_c[2,3]<-stat_table_c[3,3]-stat_table_c[1,3] 
	stat_table_c[1,1]<-length(intersect(pathway_c,up_list))
	stat_table_c[1,2]<-stat_table_c[1,3]-stat_table_c[1,1]
	stat_table_c[2,1]<-stat_table_c[3,1]-stat_table_c[1,1]
	stat_table_c[2,2]<-stat_table_c[3,2]-stat_table_c[1,2]
	p_c<-1-phyper(stat_table_c[1,1], m=stat_table_c[3,1], n=stat_table_c[3,3], k=stat_table_c[1,3])
	p_all<-c(p_all,p_c)
	s_c<-sign(stat_table_c[1,1]/stat_table_c[2,1]-stat_table_c[1,2]/stat_table_c[2,2])
	sign_all<-c(sign_all,s_c)
	if(i%%100==1)
	{
		print(i)
	}
}

results<-cbind(sign_all,p_all)
rownames(results)<-names(all_pathways)
PE<-results[which(((results[,1]==1)&(results[,2]<0.05))),]
PE<-PE[order(PE[,2]),]
