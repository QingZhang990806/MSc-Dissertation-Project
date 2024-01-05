########## Preparation ##########
get.adjacency.biogrid <- function(tmp.biogrid, names.genes = NULL){
  if(is.null(names.genes)){
    names.genes <- sort(union(unique(tmp.biogrid[,"Official.Symbol.Interactor.A"]),
                              unique(tmp.biogrid[,"Official.Symbol.Interactor.B"])))
    ind <- seq(1,nrow(tmp.biogrid))}
  else{
    ind.A <- which(tmp.biogrid[,"Official.Symbol.Interactor.A"]%in%names.genes)
    ind.B <- which(tmp.biogrid[,"Official.Symbol.Interactor.B"]%in%names.genes)
    ind <- intersect(ind.A, ind.B)
  }
  mat.biogrid <- matrix(0, nrow=length(names.genes), ncol=length(names.genes), dimnames=list(names.genes, names.genes))
  for(i in ind){
    mat.biogrid[tmp.biogrid[i,"Official.Symbol.Interactor.A"], tmp.biogrid[i,"Official.Symbol.Interactor.B"]] <- mat.biogrid[tmp.biogrid[i,"Official.Symbol.Interactor.B"], tmp.biogrid[i,"Official.Symbol.Interactor.A"]] <- 1
  }
  diag(mat.biogrid) <- 0
  return(mat.biogrid)
}

## GRN Reconstruction Validation
### Read BioGrid database
## Stark, C., Breitkreutz, B. J., Reguly, T., Boucher, L., Breitkreutz, A., & Tyers, M. (2006). BioGRID: a general repository for interaction datasets. Nucleic acids research, 34(Database issue), D535–D539. https://doi.org/10.1093/nar/gkj109
### Check last version in https://thebiogrid.org/download.php 
## Latest Version: 4.4.212
file <- "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-4.4.212/BIOGRID-ALL-4.4.212.tab2.zip"
if(!file.exists(gsub("zip","txt",basename(file)))){
  downloader::download(file,basename(file))
  unzip(basename(file),junkpaths =TRUE)
}

tmp.biogrid <- read.csv(gsub("zip","txt",basename(file)), header=TRUE, sep="\t", stringsAsFactors=FALSE)

mycols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')



########## GRN for Microarray ##########
rt=read.table("diffcamcapExp.txt",sep="\t",header=T,row.names=1,check.names=F)
names.genes.de <- rownames(rt) # 63 DEGs

net.biogrid.de <- get.adjacency.biogrid(tmp.biogrid, names.genes.de) # Unweighted Adjacency Matrix: 63*63
a <- which(net.biogrid.de == 1)
## 12

rt1=as.matrix(rt) 
mydata <- rt1  ## DEG_GEM

### GRNs
t.mydata <- t(mydata)
net.aracne <- minet(t.mydata, method = "aracne")
net.mrnet <- minet(t.mydata, method = "mrnet")
net.clr <- minet(t.mydata, method = "clr")
net.c3net <- c3net(mydata)
library(Rgraphviz)
plot( as( net.aracne ,"graphNEL") )

### Validation
tmp.val <- list(
  validate(net.aracne, net.biogrid.de), 
  validate(net.mrnet, net.biogrid.de),
  validate(net.clr, net.biogrid.de), 
  validate(net.c3net, net.biogrid.de)
)

### ROC Curve
dev1 <- show.roc(tmp.val[[1]],cex=0.3,col=mycols[1],type="l")
res.auc <- auc.roc(tmp.val[[1]])
for(count in 2:length(tmp.val)){
  show.roc(tmp.val[[count]],device=dev1,cex=0.3,col=mycols[count],type="l")
  res.auc <- c(res.auc, auc.roc(tmp.val[[count]]))
}

legend("bottomright", legend=paste(c("aracne","mrnet","clr","c3net"), signif(res.auc,4), sep=": "),
       col=mycols[1:length(tmp.val)],lty=1, bty="n" )




########## GRN for RNA-seq ##########
diffExp=rbind(id=c(Vennplot@IntersectionSets$`11`),TCGA_Pseudo[c(Vennplot@IntersectionSets$`11`),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F) # Expression Matrix for DEGs
rt=read.table("diffmRNAExp.txt",sep="\t",header=T,row.names=1,check.names=F)
names = rt[,1] # 1187 DEGs
names.genes.de <- names # 1847 DEGs
net.biogrid.de <- get.adjacency.biogrid(tmp.biogrid, names.genes.de) # Unweighted Adjacency Matrix: 1847*1847
a <- which(net.biogrid.de == 1)
# There are 6,494 unique interactions between the 1,847 DEGs on RNA-seq data.

rt = rt[,-1]
rt1=as.matrix(rt)
colnames(rt1) <- colnames(TCGA)
mydata <- rt1  ###  DEG_GEM


### GRNs
t.mydata <- t(mydata)
net.aracne <- minet(t.mydata, method = "aracne")
net.mrnet <- minet(t.mydata, method = "mrnet")
net.clr <- minet(t.mydata, method = "clr")
net.c3net <- c3net(mydata)
library(Rgraphviz)
plot( as( net.aracne ,"graphNEL") )

### Validation compared to biogrid network
tmp.val <- list(
  validate(net.aracne, net.biogrid.de), 
  validate(net.mrnet, net.biogrid.de),
  validate(net.clr, net.biogrid.de), 
  validate(net.c3net, net.biogrid.de)
)

### ROC Curve
dev1 <- show.roc(tmp.val[[1]],cex=0.3,col=mycols[1],type="l")
res.auc <- auc.roc(tmp.val[[1]])
for(count in 2:length(tmp.val)){
  show.roc(tmp.val[[count]],device=dev1,cex=0.3,col=mycols[count],type="l")
  res.auc <- c(res.auc, auc.roc(tmp.val[[count]]))
}

legend("bottomright", legend=paste(c("aracne","mrnet","clr","c3net"), signif(res.auc,4), sep=": "),
       col=mycols[1:length(tmp.val)],lty=1, bty="n" )








########## GRN for Protein ##########
names.genes.de <- rownames(Protein2_GEM) # 879 DEGs

net.biogrid.de <- get.adjacency.biogrid(tmp.biogrid, names.genes.de) # Unweighted Adjacency Matrix: 879*879
a <- which(net.biogrid.de == 1)
## 28468

mydata <- Protein2_GEM  ### DEG_GEM

### GRNs
t.mydata <- t(mydata)
net.aracne <- minet(t.mydata, method = "aracne")
net.mrnet <- minet(t.mydata, method = "mrnet")
net.clr <- minet(t.mydata, method = "clr")
net.c3net <- c3net(mydata)
library(Rgraphviz)
plot( as( net.aracne ,"graphNEL") )

### Validation
tmp.val <- list(
  validate(net.aracne, net.biogrid.de), 
  validate(net.mrnet, net.biogrid.de),
  validate(net.clr, net.biogrid.de), 
  validate(net.c3net, net.biogrid.de)
)

### ROC Curve
dev1 <- show.roc(tmp.val[[1]],cex=0.3,col=mycols[1],type="l")
res.auc <- auc.roc(tmp.val[[1]])
for(count in 2:length(tmp.val)){
  show.roc(tmp.val[[count]],device=dev1,cex=0.3,col=mycols[count],type="l")
  res.auc <- c(res.auc, auc.roc(tmp.val[[count]]))
}

legend("bottomright", legend=paste(c("aracne","mrnet","clr","c3net"), signif(res.auc,4), sep=": "),
       col=mycols[1:length(tmp.val)],lty=1, bty="n" )
