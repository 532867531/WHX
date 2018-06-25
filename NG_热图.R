dir="G:/KJ-DB-PM-BJ161062-07-2018-P11776_20180525_PM-BJ161062-07中国农业大学4个人有参转录组建库测序分析任务单/upload/4_de/2_singleDE_table/"
(files=list.files(path = dir,pattern = ".*\\.csv",full.names = T))
library(NGCHM)
thefile1=read.csv(file = files[1],sep = "\t",stringsAsFactors = F)
colnames(thefile1)[which(grepl(pattern = "ENSG",x = thefile1[10,],perl = T))]="Gene"
thefile1=thefile1[,which(grepl(x=colnames(thefile1),pattern = "^Gene$|^GeneName$|_count|_normalize",perl = TRUE))]
colnames(thefile1)[which(colnames(thefile1)!="Gene")]=paste(colnames(thefile1)[which(colnames(thefile1)!="Gene")],"thefile1",sep = "_")
thefile2=read.csv(file = files[2],sep = ",",stringsAsFactors = F)
colnames(thefile2)[which(grepl(pattern = "ENSG",x = thefile1[10,],perl = T))]="Gene"
thefile2=thefile2[,which(grepl(x=colnames(thefile2),pattern = "^Gene$|^GeneName$|_count|_normalize",perl = TRUE))]
colnames(thefile2)[which(colnames(thefile2)!="Gene")]=paste(colnames(thefile2)[which(colnames(thefile2)!="Gene")],"thefile2",sep = "_")
thefile3=read.csv(file = files[3],sep = ",",stringsAsFactors = F)
colnames(thefile3)[which(grepl(pattern = "ENSG",x = thefile1[10,],perl = T))]="Gene"
thefile3=thefile3[,which(grepl(x=colnames(thefile3),pattern = "^Gene$|^GeneName$|_count|_normalize",perl = TRUE))]
colnames(thefile3)[which(colnames(thefile3)!="Gene")]=paste(colnames(thefile3)[which(colnames(thefile3)!="Gene")],"thefile3",sep = "_")
thefile4=read.csv(file = files[4],sep = ",",stringsAsFactors = F)
colnames(thefile4)[which(grepl(pattern = "ENSG",x = thefile1[10,],perl = T))]="Gene"
thefile4=thefile4[,which(grepl(x=colnames(thefile4),pattern = "^Gene$|^GeneName$|_count|_normalize",perl = TRUE))]
colnames(thefile4)[which(colnames(thefile4)!="Gene")]=paste(colnames(thefile4)[which(colnames(thefile4)!="Gene")],"thefile4",sep = "_")
thefile5=read.csv(file = files[5],sep = ",",stringsAsFactors = F)
colnames(thefile5)[which(grepl(pattern = "ENSG",x = thefile1[10,],perl = T))]="Gene"
thefile5=thefile5[,which(grepl(x=colnames(thefile5),pattern = "^Gene$|^GeneName$|_count|_normalize",perl = TRUE))]
colnames(thefile5)[which(colnames(thefile5)!="Gene")]=paste(colnames(thefile5)[which(colnames(thefile5)!="Gene")],"thefile5",sep = "_")
###merge
thefile=merge(x=thefile1,thefile2,by.x = "Gene",by.y = "Gene",all = TRUE)
thefile=merge(x=thefile,thefile3,by.x = "Gene",by.y = "Gene",all = TRUE)
thefile=merge(x=thefile,thefile4,by.x = "Gene",by.y = "Gene",all = TRUE)
thefile=merge(x=thefile,thefile5,by.x = "Gene",by.y = "Gene",all = TRUE)
rownames(thefile)=thefile$Gene
##排序对比
thefile=thefile[,order(colnames(thefile))]
thefile_trim1=thefile[,which(grepl(pattern ="GeneName_thefile1|normalize" ,x = colnames(thefile),perl = T))]
A=stringr::str_split(string = colnames(thefile_trim1),pattern = "_",simplify = TRUE)
B=paste(A[,1],A[,2],sep  = "_")
C=which(!duplicated(B))
thefile_trim2=thefile_trim1[,C]
##开始准备画热图
rownames(thefile_trim2)=paste(rownames(thefile_trim2),thefile_trim2$GeneName_thefile1,sep = "_")
thefile_trim3=thefile_trim2[,-which(grepl(x = colnames(thefile_trim2),pattern = "GeneName"))]
##查看值的分布
library(ggpubr)
ggpubr::ggdensity(data=unlist(thefile_trim3))
(quantile(unlist(unlist(thefile_trim3)),na.rm = TRUE))
thefile_trim4=thefile_trim3[which(apply(X = thefile_trim3,MARGIN = 1,FUN = function(x){sum(is.na(x))<4})),]
thefile_trim5=thefile_trim4[apply(X = thefile_trim4,MARGIN = 1,FUN = function(x){sum(x<100|is.na(x))<4}),]

######开始画热图
library(NGCHM)
mat=as.matrix(thefile_trim5)
#saveRDS(object = mat,file = "mat.RDS")
mat=readRDS(file = "mat.RDS")
##adjust for the color
PERCENTILE=0.20;lowQ=as.numeric(quantile(unlist(mat),PERCENTILE,na.rm = TRUE));highQ=as.numeric(quantile(unlist(mat),1-PERCENTILE,na.rm = TRUE))
BREAKS=c(min(mat)-20,seq(lowQ,highQ,10),max(mat)+20)
col=colorRampPalette(c("green", "black", "red"))(length(BREAKS))

Agglomerations=c("ward.D", "ward.D2", "single", "complete", "average" , "mcquitty", "median" ,"centroid")
Dists=c( "euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
(Combination_Agglomerations=merge(x = Agglomerations,y = Agglomerations))
(Combination_Dists=merge(x = Dists,y = Dists))


apply(X = Combination_Agglomerations,MARGIN = 1 ,FUN = function(agglos){
  apply(X = Combination_Dists,MARGIN = 1,FUN = function(dists){
    thename=paste("Row",agglos[1],dists[1],"Column",agglos[2],dists[2],sep = "_")    
    thename=stringi::stri_replace_all(str = thename,replacement = "_",regex = "\\.")
    message(thename)
    cmap1 <- chmNewColorMap (BREAKS,col)
    layer1 <- chmNewDataLayer ('layer.name', mat, cmap1)
    hm <- chmNew (paste('WHX',thename,sep = "_"),layer1,rowAgglom = agglos[1],rowDist = dists[1],colAgglom = agglos[2],colDist = dists[2])
    # chmExportToFile(hm,paste('WHX',thename,'.ngchm',sep = "_"),shaidyMapGen = "./ShaidyMapGen.jar",overwrite = T)
    # chmExportToPDF(hm,filename = paste('WHX',thename,'.pdf',sep = "_"),shaidyMapGen = "./ShaidyMapGen.jar",overwrite = T)
    chmExportToFile(hm,paste('WHX',thename,'.ngchm',sep = "_"),shaidyMapGen = "./ShaidyMapGen.jar",overwrite = T)
    chmExportToPDF(hm,filename = paste('WHX',thename,'.pdf',sep = "_"),shaidyMapGen = "./ShaidyMapGen.jar",overwrite = T)
  })
})




    
