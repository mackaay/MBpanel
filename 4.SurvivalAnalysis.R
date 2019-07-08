library(factoextra)
library(NbClust)
library(cluster)
library(survival)
library(survminer)
library(forestplot)
library(openxlsx)
library(survivalROC)

command <- matrix(
    c("mRNAsi","R",1,"character",
      "mDNAsi","D",1,"character",
      "Groups","G",1,"character",
      "immuneInfo","I",1,"character",
      "cutoff","C",1,"integer",
      "optimal.cut",1,"logical",
      "outdir","O",1,"character")
      ,byrow=T,ncol=4
)
opt <- getopt(command)

da_mRNA <- read.table(opt$mRNAsi,sep="\t",header=F,stringsAsFactors=F,check.names=F,row.names=1)
da_mDNA <- read.table(opt$mDNAsi,sep="\t",header=F,stringsAsFactors=F,check.names=F,row.names=1)
mat.3 <- read.table(opt$immuneInfo,sep="\t",header=T,row.names=1,stringsAsFactors=F,check.names=F,skip = 2)
samples <- rownames(da_mRNA)

group_info <- read.table(opt$Groups,header=T,sep="\t",stringsAsFactors=F,row.names=1)
groups <- group_info[samples,]
mat.3 <- as.data.frame(t(mat.3[,samples]))
da <- as.data.frame(cbind(da_mRNA,da_mDNA,mat.3$ImmuneScore,groups))
colnames(da)[1] <- "mRNAsi"
colnames(da)[2] <- "mDNAsi"
colnames(da)[3] <- "ImmuneScore"
colnames(da)[8] <- "MetStatus"
colnames(da)[10] <- "OS.years"
da <- da[!(is.na(da$Dead) | is.na(da$OS.years)),]
SHH <- subset(da,Subgroup=="SHH")
Group3 <- subset(da,Subgroup=="Group3")
Group4 <- subset(da,Subgroup=="Group4")
WNT <- subset(da,Subgroup=="WNT")

plot.shiny.km <- function(time, death, x, cut, title = "", ids = NULL,
                          subset = rep(TRUE, length(time)),
                          col = NULL,  xlab = NULL, ylab = NULL, hr.inverse = FALSE,
                          no.plot = FALSE, optimal.cut = FALSE, assign.val=NULL,...) {
  ## filter out missing data ##
  subset = subset & !is.na(time) & !is.na(death) & !is.na(x) & !is.nan(x)
  x = x[subset]; time = time[subset]; death = death[subset];
  if (!is.null(ids)) ids = ids[subset]
  if (length(x) ==0 | length(time) == 0 | length(death) == 0
      | length(unique(death)) < 2) {
    return(invisible(NULL))
  }
  if (is.null(ids)) ids = 1:length(x)
  km.data = data.frame(ID = ids, X = x, Group = NA, time = time, event = death)
  ## settings for median cutoff ##
  #cut = median(x)
  p.adj = NULL
  if(is.numeric(x)){
    #upper = "upper 50%"; lower = "lower 50%"
    #upper = paste0("upper ", 100-cut, "%")
    #lower = paste0("lower ", cut, "%")
    upper = "high risk group"
    lower = "low risk group"
    cut = quantile(x,cut/100)
    ## find optimal cutoff if specified ##
    if (optimal.cut) {
      mod <- coxph(Surv(time, event) ~ X, data = km.data)
      cc = try(cutp(mod), silent = TRUE)
      if (class(cc) %in% "try-error") {
        return(invisible(NULL))
      }
      cut = cc$X$X[1]
      # assign(assign.val,cut,envir = .GlobalEnv)
      p.adj = cc$X$p[1]
      percentile = round(sum(x>=cut) / length(x) * 100,2)
      #upper = paste0("upper ", percentile, "%")
      #lower = paste0("lower ", 100-percentile, "%")
    }
    
    if(!is.null(assign.val)){
      assign(assign.val,cut,envir = .GlobalEnv)
      ## split into high and low groups using appropriate cutoff ##
    }
    newX = x
    newX[x > cut] = upper
    newX[x <= cut] = lower
    expression = factor(newX,levels=c("low risk group","high risk group"))
    km.data$Group = expression
  }else{
    km.data$Group = as.factor(x)
  }
  if (no.plot) return(invisible(km.data))
  n = length(levels(km.data$Group))
  km.group = try(survdiff(Surv(time, death) ~ as.numeric(km.data$Group)), silent = TRUE)
  km.group1 = try(coxph(Surv(time, death) ~ as.numeric(km.data$Group)), silent = TRUE)
  if (class(km.group) %in% "try-error") return(invisible(NULL))
  if (class(km.group1) %in% "try-error") return(invisible(NULL))
  p.km = 1 - pchisq(km.group$chisq,n-1)
  hr = exp(km.group1$coefficients)
  conf = paste0("(",signif(summary(km.group1)$conf.int[,"lower .95"],2),"-",signif(summary(km.group1)$conf.int[,"upper .95"],2),")")
  n = km.group1$n
  
  if (hr.inverse) hr = 1/hr
  
  hr.str = paste0("HR=", round(hr,2),conf,", ")
  p.str = paste0("P=", round(p.km,4))
  if (!is.null(p.adj)) {
    p.str = paste0(p.str, ", P(cutoff) = ", round(p.adj,4))
  }
  
  if (title=="") {
    title = paste(hr.str, p.str, sep = "")
  } else {
    title = paste0(title, "\n",hr.str, p.str)
  }
  
  if (is.null(xlab)) xlab = "Time"
  if (is.null(ylab)) ylab = "Survival"
  
  ## plot graph ### ggplot2/GGally form
  km.group1 = survfit(Surv(time, death) ~ Group,data=km.data)
  km.group = ggsurv(km.group1,
                    main = title,
                    xlab = xlab, ylab = ylab,...) +
    ggplot2::coord_cartesian(ylim = c(0, 1))
  #    theme(plot.title = element_text(vjust = 0.5, colour = "black"))
  
  plot(km.group)
  return(invisible(km.data))
}

for(i in levels(factor(da$Subgroup))){
  dat <- get(i)
  pdf(paste0(opt$outdir,'/',i,"_mRNAsi_survival.pdf"))
  plot.shiny.km(dat$OS.years,dat$Dead,dat$mRNAsi,cut=50,optimal.cut = opt$optimal.cut,assign.val=paste0(i,"_cut"),title = paste0(i,"_mRNAsi"),xlab="time (years)",ylab="OS probability")
  dev.off()
}

SHH_mRNAsi <- ifelse(SHH$mRNAsi>SHH_cut,"SHH_upper","SHH_lower")
Group3_mRNAsi <- ifelse(Group3$mRNAsi>Group3_cut,"Group3_upper","Group3_lower")
Group4_mRNAsi <- ifelse(Group4$mRNAsi>Group4_cut,"Group4_upper","Group4_lower")
WNT_mRNAsi <- ifelse(WNT$mRNAsi>WNT_cut,"WNT_upper","WNT_lower")
SHH <- cbind(SHH,SHH_mRNAsi)
Group3 <- cbind(Group3,Group3_mRNAsi)
Group4 <- cbind(Group4,Group4_mRNAsi)
WNT <- cbind(WNT,WNT_mRNAsi)

diff_group_info <- data.frame(samples=rownames(SHH),Groups=as.character(SHH$SHH_mRNAsi))
sample_info <- cbind(rep('defaultExp',length(diff_group_info$Groups)),rep('RNA-seq',length(diff_group_info$Groups)),rep('',length(diff_group_info$Groups)),as.character(diff_group_info$samples),as.character(diff_group_info$Groups),rep('',length(diff_group_info$Groups)),rep('',length(diff_group_info$Groups)),rep('',length(diff_group_info$Groups)),rep('',length(diff_group_info$Groups)))
colnames(sample_info) <- c("实验名称","数据类型","实验描述","样本名称","分组","描述","文件名","md5","保存路径")
write.xlsx(sample_info,file=paste(opt$outdir,'/','diff_group_info.xlsx',sep=''))

cut_mRNAsi <- cut(da$mRNAsi,quantile(da$mRNAsi,probs=c(0,0.25,0.5,0.75,1)))
cut_mDNAsi <- cut(da$mDNAsi,quantile(da$mDNAsi,probs=c(0,0.25,0.5,0.75,1)))
da <- cbind(da,cut_mRNAsi,cut_mDNAsi)
da1 <- da[!is.na(da$MetStatus),]
da1$MetStatus <- ifelse(da1$MetStatus==1,"metastatic","non-metastatic")

write2table <- function(x){
  df <- get(x)
  df <- as.data.frame(cbind(rownames(df),df))
  colnames(df)[1] <- "features"
  write.table(df,paste0(x,".txt"),row.names = F,col.names = T,quote=F,sep="\t")
}

StemScore_clinical_info <- da
write2table("StemScore_clinical_info")

#cox regresssion model
plotHR <- function(text,res,prefix){
  text[1,] <- gsub("univ_|multi_","",colnames(text))
  png(paste0(opt$outdir,'/',prefix,'_forestplot.png'),width=960, height=(1+dim(res)[1])*50)
  colinfo <- lapply(1:nrow(res),function(x){ifelse(x %% 2==1,
                                                   return(gpar(col="gray",lineend="butt",columns=c(2:6))),return(gpar(col='gray',lineend="butt",columns=c(2:6))))})
  names(colinfo) <- 1:nrow(res)+1
  forestplot(as.matrix(text),mean=c(NA,res$mean),lower = c(NA,res$lower),upper = c(NA,res$upper),
             graph.pos = 4,graphwidth = unit(60,'mm'),lineheight = unit(12,'mm'),line.margin = unit(5,'mm'),
             colgap = unit(4,'mm'),zero = 1,cex=0.9,col=fpColors(box="black",lines="black",zero="gray50"),
             txt_gp = fpTxtGp(label=gpar(cex=1.25),ticks=gpar(cex=1.1),xlab=gpar(cex=1.2),title = gpar(cex=1.2)),
             lwd.ci=1.5,ci.vertices = T,title="Hazard Ratio",boxsize = 0.3,hrzl_lines = colinfo)
  dev.off()
}

univ_analysis <- function(features,da){
  univ_formulas <- sapply(features,function(x) as.formula(paste('Surv(OS.years, Dead)~', paste0('`',x,'`'))))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = da)})
  univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"])
                         beta<-signif(x$coef[1]);#coeficient beta
                         HR <-signif(x$coef[2]);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"])
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"])
                         CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                         res<-c(beta, HR, CI, p.value)
                         names(res)<-c("univ_beta", "univ_HR", "univ_95% CI for HR", "univ_p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
  res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  res$univ_beta <- sapply(res$univ_beta,function(x) as.numeric(as.character(x)))
  res$univ_HR <- sapply(res$univ_HR,function(x) as.numeric(as.character(x)))
  res$`univ_95% CI for HR` <- sapply(res$`univ_95% CI for HR`,as.character)
  res$`univ_p.value` <- sapply(res$`univ_p.value`,function(x) as.numeric(as.character(x)))
  return(res)
}

multi_cox <- function(sig_features,data){
  sig_features_formula <- unlist(sapply(sig_features,function(x){paste0('`',x,'`')}))
  multi_formula <- as.formula(paste('Surv(OS.years, Dead) ~',paste(sig_features_formula,collapse='+')))
  multi_models <- coxph(multi_formula, data = data)
  suma <- summary(multi_models)
  rownames(suma$coef) <- sig_features
  rownames(suma$conf.int) <- sig_features
  multi_res <- sapply(sig_features,function(x){
    p.value <- signif(suma$coef[x,'Pr(>|z|)'])
    beta <- signif(suma$coef[x,'coef'])
    HR <- signif(suma$coef[x,'exp(coef)'])
    HR.confint.lower <- signif(suma$conf.int[x,"lower .95"])
    HR.confint.upper <- signif(suma$conf.int[x,"upper .95"])
    CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
    res<-c(beta, HR, CI, p.value)
    names(res)<-c("multi_beta", "multi_HR", "multi_95% CI for HR", "multi_p.value")
    return(res)
  })
  multi_res <- as.data.frame(t(multi_res))
  multi_res$multi_beta <- sapply(multi_res$multi_beta,function(x) as.numeric(as.character(x)))
  multi_res$multi_HR <- sapply(multi_res$multi_HR,function(x) as.numeric(as.character(x)))
  multi_res$`multi_95% CI for HR` <- sapply(multi_res$`multi_95% CI for HR`,as.character)
  multi_res$`multi_p.value` <- sapply(multi_res$`multi_p.value`,function(x) as.numeric(as.character(x)))
  return(multi_res)
}

forest_plot <- function(res,prefix){
  result <- cbind(rownames(res),res)
  colnames(result)[1] <- "features"
  write.table(result,paste0(opt$outdir,'/',prefix,'_cox_univ_res.txt'),sep="\t",quote=F,row.names=F)
  res$univ_beta <- signif(res$univ_beta,2)
  res$univ_HR <- signif(res$univ_HR,2)
  res$`univ_95% CI for HR` <- sapply(res$`univ_95% CI for HR`,function(x){paste0(signif(as.numeric(unlist(strsplit(x,'-'))),2),collapse = "-")})
  res$univ_p.value <- signif(res$univ_p.value,4)
  res$mean <- as.numeric(as.character(res$univ_HR))
  res$lower <- as.numeric(as.character(sapply(res$`univ_95% CI for HR`,function(x){unlist(strsplit(as.character(x),'-'))[1]})))
  res$upper <- as.numeric(as.character(sapply(res$`univ_95% CI for HR`,function(x){unlist(strsplit(as.character(x),'-'))[2]})))
  text <- rbind(colnames(res)[1:4],as.matrix(res[1:4]))
  rownames(text)[1] <- 'features'
  text <- cbind(rownames(text),text)
  plotHR(text,res,prefix)
}

forest_plot_multi <- function(res,prefix){
  result <- cbind(rownames(res),res)
  colnames(result)[1] <- "features"
  write.table(result,paste0(opt$outdir,'/',prefix,'_cox_multi_res.txt'),sep="\t",quote=F,row.names=F)
  res$multi_beta <- signif(res$multi_beta,2)
  res$multi_HR <- signif(res$multi_HR,2)
  res$`multi_95% CI for HR` <- sapply(res$`multi_95% CI for HR`,function(x){paste0(signif(as.numeric(unlist(strsplit(x,'-'))),2),collapse = "-")})
  res$multi_p.value <- signif(res$multi_p.value,4)  
  res$mean <- as.numeric(as.character(res$multi_HR))
  res$lower <- as.numeric(as.character(sapply(res$`multi_95% CI for HR`,function(x){unlist(strsplit(as.character(x),'-'))[1]})))
  res$upper <- as.numeric(as.character(sapply(res$`multi_95% CI for HR`,function(x){unlist(strsplit(as.character(x),'-'))[2]})))
  text <- rbind(colnames(res)[1:4],as.matrix(res[1:4]))
  rownames(text)[1] <- 'features'
  text <- cbind(rownames(text),text)
  plotHR(text,res,prefix)
}

features <- c("mRNAsi","mDNAsi","ImmuneScore","Age","Gender","MetStatus","Subgroup","Histology")
dat <- da
dat$Gender <- as.numeric(factor(dat$Gender))
dat$MetStatus <- as.numeric(factor(dat$MetStatus))
dat$Subgroup <- as.numeric(factor(dat$Subgroup))
dat$Histology <- as.numeric(factor(dat$Histology))
clinical_univ_res <- univ_analysis(features,dat)
forest_plot(clinical_univ_res,"clinical_info")

SHH_diff_genes <- read.table("SHH_upperVSSHH_lower.result.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
expMat <- read.table("Medulloblastoma_data/GSE85217_rna_exp_matrix.txt",sep="\t",header=T,row.names = 1)
SHH_diff_genes <- rownames(SHH_diff_genes)[SHH_diff_genes$FDR<0.05]
SHH_mRNA_genes <- read.table("Medulloblastoma_data/GSE85217_mRNA.txt",sep="\t",header=T,row.names = 1)
SHH_mRNA_diff_genes <- intersect(rownames(SHH_mRNA_genes),SHH_diff_genes)
SHH_miRNA_genes <- read.table("Medulloblastoma_data/GSE85217_miRNA.txt",sep="\t",header=T,row.names = 1)
SHH_miRNA_diff_genes <- intersect(rownames(SHH_miRNA_genes),SHH_diff_genes)
SHH_lncRNA_genes <- read.table("Medulloblastoma_data/GSE85217_lncRNA.txt",sep="\t",header=T,row.names = 1)
SHH_lncRNA_diff_genes <- intersect(rownames(SHH_lncRNA_genes),SHH_diff_genes)
SHH_mRNA_diff_da <- t(expMat[SHH_mRNA_diff_genes,rownames(SHH)])
SHH_miRNA_diff_da <- t(expMat[SHH_miRNA_diff_genes,rownames(SHH)])
SHH_lncRNA_diff_da <- t(expMat[SHH_lncRNA_diff_genes,rownames(SHH)])
SHH_mRNA_diff_da <- as.data.frame(cbind(SHH$OS.years,SHH$Dead,SHH_mRNA_diff_da))
SHH_miRNA_diff_da <- as.data.frame(cbind(SHH$OS.years,SHH$Dead,SHH_miRNA_diff_da))
SHH_lncRNA_diff_da <- as.data.frame(cbind(SHH$OS.years,SHH$Dead,SHH_lncRNA_diff_da))
colnames(SHH_mRNA_diff_da)[1:2] <- c("OS.years","Dead")
colnames(SHH_miRNA_diff_da)[1:2] <- c("OS.years","Dead")
colnames(SHH_lncRNA_diff_da)[1:2] <- c("OS.years","Dead")

SHH_mRNA_res <- univ_analysis(features = SHH_mRNA_diff_genes,da = SHH_mRNA_diff_da)
SHH_miRNA_res <- univ_analysis(features = SHH_miRNA_diff_genes,da = SHH_miRNA_diff_da)
SHH_lncRNA_res <- univ_analysis(features = SHH_lncRNA_diff_genes,da = SHH_lncRNA_diff_da)

mRNA_cor <- apply(SHH_mRNA_diff_da[,3:dim(SHH_mRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$estimate)
mRNA_p <- apply(SHH_mRNA_diff_da[,3:dim(SHH_mRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$p.value)
miRNA_cor <- apply(SHH_miRNA_diff_da[,3:dim(SHH_miRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$estimate)
miRNA_p <- apply(SHH_miRNA_diff_da[,3:dim(SHH_miRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$p.value)
lncRNA_cor <- apply(SHH_lncRNA_diff_da[,3:dim(SHH_lncRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$estimate)
lncRNA_p <- apply(SHH_lncRNA_diff_da[,3:dim(SHH_lncRNA_diff_da)[2]],2,function(x) cor.test(x,SHH$mRNAsi)$p.value)

SHH_mRNA_res <- as.data.frame(cbind(SHH_mRNA_res,mRNA_cor,mRNA_p))
SHH_miRNA_res <- as.data.frame(cbind(SHH_miRNA_res,miRNA_cor,miRNA_p))
SHH_lncRNA_res <- as.data.frame(cbind(SHH_lncRNA_res,lncRNA_cor,lncRNA_p))
write2table("SHH_mRNA_res")
write2table("SHH_lncRNA_res")
write2table("SHH_miRNA_res")

SHH_mRNA_panel <- rownames(SHH_mRNA_res)[SHH_mRNA_res$univ_p.value < 0.05 & abs(SHH_mRNA_res$mRNA_cor) > 0.6 & SHH_mRNA_res$mRNA_p < 0.05]
SHH_mRNA_panel_r1 <- rownames(SHH_mRNA_res)[SHH_mRNA_res$univ_p.value < 0.05 & abs(SHH_mRNA_res$mRNA_cor) > 0.65 & SHH_mRNA_res$mRNA_p < 0.05]
SHH_mRNA_panel_multi_res <- multi_cox(SHH_mRNA_panel,data = SHH_mRNA_diff_da)
SHH_mRNA_panel_multi_res_r1 <- multi_cox(SHH_mRNA_panel_r1,data = SHH_mRNA_diff_da)
forest_plot_multi(SHH_mRNA_panel_multi_res,"SHH_mRNA_panel")
forest_plot_multi(SHH_mRNA_panel_multi_res_r1,"SHH_mRNA_panel_r1")
SHH_mRNA_risk_score <- NULL
SHH_mRNA_risk_score_r1 <- NULL
for(i in 1:dim(SHH_mRNA_diff_da)[1]){
  risk_score <- sum(SHH_mRNA_panel_multi_res$multi_beta * SHH_mRNA_diff_da[i,rownames(SHH_mRNA_panel_multi_res)])
  names(risk_score) <- rownames(SHH_mRNA_diff_da)[i]
  SHH_mRNA_risk_score <- c(SHH_mRNA_risk_score,risk_score)
}

for(i in 1:dim(SHH_mRNA_diff_da)[1]){
  risk_score <- sum(SHH_mRNA_panel_multi_res_r1$multi_beta * SHH_mRNA_diff_da[i,rownames(SHH_mRNA_panel_multi_res_r1)])
  names(risk_score) <- rownames(SHH_mRNA_diff_da)[i]
  SHH_mRNA_risk_score_r1 <- c(SHH_mRNA_risk_score_r1,risk_score)
}

pdf(paste0(opt$outdir,'/SHH_mRNA_risk_kmplot.pdf'))
plot.shiny.km(SHH$OS.years,SHH$Dead,SHH_mRNA_risk_score,cut=50,optimal.cut = T ,assign.val=NULL,title = "7 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
dev.off()

SHH_lncRNA_panel <- rownames(SHH_lncRNA_res)[SHH_lncRNA_res$univ_p.value<0.05]
SHH_lncRNA_panel_multi_res <- multi_cox(SHH_lncRNA_panel,data = SHH_lncRNA_diff_da)
forest_plot_multi(SHH_lncRNA_panel_multi_res,"SHH_lncRNA_panel")
SHH_lncRNA_risk_score <- NULL
for(i in 1:dim(SHH_lncRNA_diff_da)[1]){
  risk <- sum(SHH_lncRNA_panel_multi_res$multi_beta * SHH_lncRNA_diff_da[i,rownames(SHH_lncRNA_panel_multi_res)])
  names(risk) <- rownames(SHH_lncRNA_diff_da)[i]
  SHH_lncRNA_risk_score <- c(SHH_lncRNA_risk_score,risk)
}

pdf(paste0(opt$outdir,'/SHH_lncRNA_risk_kmplot.pdf'))
plot.shiny.km(SHH$OS.years,SHH$Dead,SHH_lncRNA_risk_score,cut=50,optimal.cut = T,assign.val=NULL,title = "3 lncRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
dev.off()

pdf(paste0(opt$outdir,'/SHH_miRNA_exp_kmplot.pdf'))
plot.shiny.km(SHH$OS.years,SHH$Dead,SHH_miRNA_diff_da[,'MIR4453'],cut=50,optimal.cut = T,assign.val=NULL,title = "SHH_miRNA_expression",xlab="time (years)",ylab="OS probability")
dev.off()

SHH_new <- cbind(SHH,SHH_mRNA_risk_score,SHH_lncRNA_risk_score)
fe <- c("MetStatus","SHH_mRNA_risk_score","SHH_lncRNA_risk_score")
SHH_new_multi_res <- multi_cox(fe,data = SHH_new)
write2table("SHH_new_multi_res")
forest_plot_multi(SHH_new_multi_res,"SHH_risk_clinical_multi")

#scatter plot
scatter_plot <- function(gene,da,value,res,dir,cor_column=1,p_column=2){
  x <- unlist(da[,gene])
  y <- value
  lg <- vector('expression',2)
  lg[1] = substitute(expression(italic(R) == MYVALUE),list(MYVALUE = format(res[gene,cor_column],dig=3)))[2]
  lg[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(res[gene,p_column], digits = 2)))[2]
  pdf(paste0(dir,'/',gene,'_scatter_plot.pdf'))
  plot(x,y,type="p",col=densCols(x,y),pch=20,xlab=paste0(gene,"_expression"),ylab="mRNAsi",main=gene)
  abline(lm(y~x),lty=1,lwd=2,col="grey40")
  legend("topright",lg,bty="n")
  dev.off()
}

for(i in SHH_mRNA_panel){
  scatter_plot(i,SHH_mRNA_diff_da,SHH$mRNAsi,SHH_mRNA_res,"mRNA_scatter_V2",cor_column = 5,p_column = 6)
}

for(i in SHH_lncRNA_panel){
  scatter_plot(i,SHH_lncRNA_diff_da,SHH$mRNAsi,SHH_lncRNA_res,"lncRNA_scatter_V2",cor_column = 5,p_column = 6)
}

scatter_plot("MIR4453",SHH_miRNA_diff_da,SHH$mRNAsi,SHH_miRNA_res,"miRNA_scatter_V2",cor_column = 5,p_column = 6)

#ROC curve
mRNA_risk_1_year <- survivalROC(Stime = SHH_mRNA_diff_da$Time,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score,predict.time = 1,span = 0.25*length(SHH_mRNA_risk_score)^(-0.20))
mRNA_risk_3_year <- survivalROC(Stime = SHH_mRNA_diff_da$Time,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score,predict.time = 3,span = 0.25*length(SHH_mRNA_risk_score)^(-0.20))
mRNA_risk_5_year <- survivalROC(Stime = SHH_mRNA_diff_da$Time,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score,predict.time = 5,span = 0.25*length(SHH_mRNA_risk_score)^(-0.20))

mRNA_risk_1_year_r1 <- survivalROC(Stime = SHH_mRNA_diff_da$OS.years,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score_r1,predict.time = 1,span = 0.25*length(SHH_mRNA_risk_score_r1)^(-0.20))
mRNA_risk_3_year_r1 <- survivalROC(Stime = SHH_mRNA_diff_da$OS.years,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score_r1,predict.time = 3,span = 0.25*length(SHH_mRNA_risk_score_r1)^(-0.20))
mRNA_risk_5_year_r1 <- survivalROC(Stime = SHH_mRNA_diff_da$OS.years,status = SHH_mRNA_diff_da$Dead,marker = SHH_mRNA_risk_score_r1,predict.time = 5,span = 0.25*length(SHH_mRNA_risk_score_r1)^(-0.20))

lncRNA_risk_1_year <- survivalROC(Stime = SHH_lncRNA_diff_da$OS.years,status = SHH_lncRNA_diff_da$Dead,marker = SHH_lncRNA_risk_score,predict.time = 1,span = 0.25*length(SHH_lncRNA_risk_score)^(-0.20))
lncRNA_risk_3_year <- survivalROC(Stime = SHH_lncRNA_diff_da$OS.years,status = SHH_lncRNA_diff_da$Dead,marker = SHH_lncRNA_risk_score,predict.time = 3,span = 0.25*length(SHH_lncRNA_risk_score)^(-0.20))
lncRNA_risk_5_year <- survivalROC(Stime = SHH_lncRNA_diff_da$OS.years,status = SHH_lncRNA_diff_da$Dead,marker = SHH_lncRNA_risk_score,predict.time = 5,span = 0.25*length(SHH_lncRNA_risk_score)^(-0.20))

pdf("mRNA_panel_risk_score_roc.pdf")
plot(mRNA_risk_1_year$FP, mRNA_risk_1_year$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="False positive rate",ylab="True positive rate",main="mRNA panel",col="red",lwd=2,lty=1)
lines(mRNA_risk_3_year$FP, mRNA_risk_3_year$TP,type="l",col="blue",lwd=2,lty=2)
lines(mRNA_risk_5_year$FP, mRNA_risk_5_year$TP,type="l",col="black",lwd=2,lty=5)
lg1 = paste0("one year AUC","=",format(mRNA_risk_1_year$AUC,dig=3))
lg2 = paste0("three year AUC","=",format(mRNA_risk_3_year$AUC,dig=3))
lg3 = paste0("five year AUC","=",format(mRNA_risk_5_year$AUC,dig=3))
legend("topleft",legend=c(lg1,lg2,lg3),bty="n",col=c("red","blue","black"),lty=c(1,2,5))
abline(0,1,lwd=1)
dev.off()

pdf("mRNA_panel_risk_score_roc_r1.pdf")
plot(mRNA_risk_1_year_r1$FP, mRNA_risk_1_year_r1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="False positive rate",ylab="True positive rate",main="mRNA panel",col="red",lwd=2,lty=1)
lines(mRNA_risk_3_year_r1$FP, mRNA_risk_3_year_r1$TP,type="l",col="blue",lwd=2,lty=2)
lines(mRNA_risk_5_year_r1$FP, mRNA_risk_5_year_r1$TP,type="l",col="black",lwd=2,lty=5)
lg1 = paste0("one year auc","=",format(mRNA_risk_1_year_r1$AUC,dig=3))
lg2 = paste0("three year auc","=",format(mRNA_risk_3_year_r1$AUC,dig=3))
lg3 = paste0("five year auc","=",format(mRNA_risk_5_year_r1$AUC,dig=3))
legend("topleft",legend=c(lg1,lg2,lg3),bty="n",col=c("red","blue","black"),lty=c(1,2,5))
abline(0,1,lwd=1)
dev.off()

pdf("lncRNA_panel_risk_score_roc.pdf")
plot(lncRNA_risk_1_year$FP, lncRNA_risk_1_year$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="False positive rate",ylab="True positive rate",main="lncRNA panel",col="red",lwd=2,lty=1)
lines(lncRNA_risk_3_year$FP, lncRNA_risk_3_year$TP,type="l",col="blue",lwd=2,lty=2)
lines(lncRNA_risk_5_year$FP, lncRNA_risk_5_year$TP,type="l",col="black",lwd=2,lty=5)
lg1 = paste0("one year auc","=",format(lncRNA_risk_1_year$AUC,dig=3))
lg2 = paste0("three year auc","=",format(lncRNA_risk_3_year$AUC,dig=3))
lg3 = paste0("five year auc","=",format(lncRNA_risk_5_year$AUC,dig=3))
legend("topleft",legend=c(lg1,lg2,lg3),bty="n",col=c("red","blue","black"),lty=c(1,2,5))
abline(0,1,lwd=1)
dev.off()

auc_optimize <- function(panel,data){
  multi_res <- multi_cox(panel,data)
  risk_score <- NULL
  for(i in 1:dim(data)[1]){
    risk <- sum(multi_res$multi_beta * data[i,rownames(multi_res)])
    names(risk) <- rownames(data)[i]
    risk_score <- c(risk_score,risk)
  }
  risk_1_year <- survivalROC(Stime = data$OS.years,status = data$Dead,marker = risk_score,predict.time = 1,span = 0.25*length(risk_score)^(-0.20))
  return(risk_1_year$AUC)
}

test <- SHH_mRNA_res[SHH_mRNA_res$univ_p.value < 0.05 & abs(SHH_mRNA_res$mRNA_cor) > 0.6 & SHH_mRNA_res$mRNA_p < 0.05,]
test <- test[order(test$univ_beta,decreasing = T),]

panel <- NULL
auc <- NULL
for(i in seq(1,27)){
  panel <- c(panel,rownames(test)[i])
  tmp <- auc_optimize(panel,SHH_mRNA_diff_da)
  auc <- c(auc,tmp)
}
