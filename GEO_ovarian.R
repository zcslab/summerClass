setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_10_18_Edward")

# BiocManager::install("GEOmetadb")

# install package 
# 1. cran default 
# 2. github
# 3. other 

library(tools)
library(GEOmetadb)
con = dbConnect(SQLite(), "GEOmetadb.sqlite")

library(plyr)
library(dplyr)
con
geo_tables = dbListTables(con)
print(geo_tables)
dbListFields(con,'gse')
dbListFields(con,'gse_gpl')
dbListFields(con,'gds')

gse_database0 = dbGetQuery(con, 'select * from gse limit 10')
gse_database = gse_database0[,-10]	#delete column of "contributor"
#table(gse_database$type)
#length(sort(table(gse_database$type),decreasing=T))
sort(table(gse_database$type),decreasing=T)[1:10]

#counts = apply(gse_database, MARGIN=1, function(x)(grep("expression profiling by high throughput sequencing",tolower(x))))
counts = apply(gse_database, MARGIN=1, function(x)(grep("Expression profiling by high throughput sequencing", x)))
#length(which(sapply(counts, length)>0))
gse_database_RNAseq = gse_database[which(sapply(counts, length)>0),]
length(table(gse_database_RNAseq$type))
#here, we got 17108 high throughput sequencing
sort(table(gse_database_RNAseq$type),decreasing=T)[1:10]
gse_database_RNAseq_unique =  gse_database_RNAseq[which(gse_database_RNAseq[,9]=="Expression profiling by high throughput sequencing"),]
dim(gse_database_RNAseq_unique)
gse_database_RNAseq_unique.gse =  gse_database_RNAseq_unique[,3]

####ovarian	
counts = apply(gse_database_RNAseq_unique , MARGIN=1, function(x)(grep("ovarian",tolower(x))))
length(which(sapply(counts, length)>0))
ovarian = gse_database_RNAseq_unique[which(sapply(counts, length)>0),]
length(ovarian$gse)
#this is GSE search from keyword "ovrian"

####tumor	
counts = apply(gse_database_RNAseq_unique , MARGIN=1, function(x)(grep("tumor",tolower(x))))
length(which(sapply(counts, length)>0))
tumor = gse_database_RNAseq_unique[which(sapply(counts, length)>0),]
length(tumor$gse)
#this is GSE search from keyword "tumor"

#print intersect of "ovarian"&"tumor"
#intersect(ovarian$gse,tumor$gse)

###cancer
counts = apply(gse_database_RNAseq_unique , MARGIN=1, function(x)(grep("cancer",tolower(x))))
length(which(sapply(counts, length)>0))
cancer = gse_database_RNAseq_unique[which(sapply(counts, length)>0),]
length(cancer$gse)
#this is GSE search from keyword "cancer"


###clear
counts = apply(gse_database_RNAseq_unique , MARGIN=1, function(x)(grep("clear",tolower(x))))
length(which(sapply(counts, length)>0))
clear = gse_database_RNAseq_unique[which(sapply(counts, length)>0),]
length(clear$gse)
#this is GSE search from keyword "clear"


aa = Reduce(intersect, list(ovarian$gse,cancer$gse,clear$gse))
bb = Reduce(intersect, list(ovarian$gse,tumor$gse,clear$gse))

unique(union(aa,bb))


##########test
#ss = as.data.frame(matrix(1:12,3,4))
#ss[1,1] = "this"
#ss[1,2] = "is"
#s[1,3] = "wennan"
#apply(ss, MARGIN=1, function(x)(grep("wennan",tolower(x))))
############
