library(getopt)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(gridExtra)

command <- matrix(
    c("sig","I",1,"character",
      "Groups","G",1,"character",
      "prefix","P",1,"character",
      "outdir","O",1,"character")
      ,byrow=T,ncol=4
)
opt <- getopt(command)

if(!dir.exists(opt$outdir)){
    dir.create(opt$outdir)
}

da <- read.table(opt$sig,sep="\t",header=F,stringsAsFactors=F,check.names=F,row.names=1)
samples <- rownames(da)

group_info <- read.table(opt$Groups,header=T,sep="\t",stringsAsFactors=F,row.names=1)
groups <- group_info[samples,]
da <- as.data.frame(cbind(da,groups))
colnames(da)[1] <- opt$prefix
colnames(da)[6] <- "MetStatus"
da$MetStatus <- ifelse(da$MetStatus==1,"metastatic","non-metastatic")
#da$MetStatus <- factor(da$MetStatus,levels = c("metastatic","non-metastatic"))
da1 <- da[!is.na(da$MetStatus),]
da1$MetStatus <- factor(da1$MetStatus,levels = c("metastatic","non-metastatic"))
da$Subgroup <- factor(da$Subgroup,levels = c("SHH","Group3","Group4","WNT"))
da1$Subgroup <- factor(da1$Subgroup,levels = c("SHH","Group3","Group4","WNT"))

#fac <- levels(factor(da$groups))
#fac_num <- nlevels(factor(da$groups))
#combines <- combn(fac_num,2)
#my_comps <- lapply(as.data.frame(combines),function(x) return(c(fac[x[1]],fac[x[2]])))
#p1 <- ggboxplot(da, x="groups", y="values", color = "groups",add = "jitter", shape="groups",xlab=opt$Xlab,ylab=opt$Ylab) + stat_compare_means(comparisons = my_comps,method="wilcox.test")
p1 <- ggboxplot(da, x="Subgroup", y=opt$prefix,title=opt$prefix, color = "Subgroup",palette = "jco",add = "jitter",add.params = list(size = 0.8, jitter = 0.1),width=0.5) + stat_compare_means()
p2 <- ggboxplot(da1, x="Subgroup", y=opt$prefix,title=opt$prefix, color = "MetStatus",palette = "jco",add = "jitter",add.params = list(size = 0.8, jitter = 0.1),width=0.5)+ stat_compare_means(aes(group=MetStatus))
p3 <- ggboxplot(da1,x="MetStatus",y=opt$prefix,title=opt$prefix,color="MetStatus",palette="jco",add="jitter",add.params = list(size = 0.8, jitter = 0.1),width=0.4) + stat_compare_means()
p4 <- ggboxplot(da1,x="MetStatus",y=opt$prefix,title=opt$prefix,color="Subgroup",palette="jco",add="jitter",add.params = list(size = 0.8, jitter = 0.1)) + stat_compare_means(aes(group=Subgroup))

ggsave(paste0(opt$outdir,'/',opt$prefix,'_Group_all_compare_boxPlot.pdf'),p1)
ggsave(paste0(opt$outdir,'/',opt$prefix,'_Group_Met_compare_boxPlot.pdf'),p2)
ggsave(paste0(opt$outdir,'/',opt$prefix,'_Met_all_compare_boxPlot.pdf'),p3)
ggsave(paste0(opt$outdir,'/',opt$prefix,'_Met_group_compare_boxPlot.pdf'),p4)
