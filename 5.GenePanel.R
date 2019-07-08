# loading packages
library(survival)
library(survminer)
library(forestplot)
library(survivalROC)
library(survMisc)
library(ggplot2)
library(GGally)
library(parallel)
library(glmnet)
library(caret)
library(dplyr)

# ---------------opt input-------------------
opt <- NULL

opt$train_prefix <- "SHH.train"
opt$outdir <- "PanelScreen"
opt$defFunctions <- TRUE
opt$ifuniv <- TRUE

# -----------------def functions-------------------

if(opt$defFunctions){
  get_feature_survival_info <- function(exp,survival,time_column,outcome_column,dead_status,alive_status){
    data.expr=read.delim(exp,row.names=1,stringsAsFactors=T,check.names=F,sep="\t")
    data.p=read.delim(survival,row.names=1,stringsAsFactors=F,check.names=F,sep="\t")
    common = intersect(colnames(data.expr), rownames(data.p))
    m1 = match(common, colnames(data.expr))
    m2 = match(common, rownames(data.p))
    data.expr = data.expr[,m1]
    data.p = data.p[m2,]
    da=data.frame(cbind(data.p[,time_column-1],data.p[,outcome_column-1],t(data.expr)),check.names = F, stringsAsFactors = F)
    colnames(da)[1:2]=c("Time","Dead")
    da$Dead[da$Dead == dead_status]=1
    da$Dead[da$Dead == alive_status]=0
    da[1:dim(da)[2]]=lapply(da[1:dim(da)[2]],FUN=function(x){as.numeric(x)})
    da = da[!(is.na(da$Time) | da$Time<=0 | is.na(da$Dead)),]
    return(da)
  }
  
  univ_analysis <- function(features,da){
    univ_formulas <- sapply(features,function(x) as.formula(paste('Surv(Time, Dead)~', paste0('`',x,'`'))))
    univ_models <- lapply( univ_formulas, function(x){coxres <- try(coxph(x, data = da),silent = T);if(class(coxres) %in% "try-error"){return(NULL)}else{return(coxres)}})
    univ_models <- univ_models[!sapply(univ_models,is.null)]
    univ_results <- lapply(univ_models,
                           function(x){
                             suma <- summary(x)
                             res <- sapply(rownames(suma$coef),function(x){
                               p.value<-signif(suma$coef[x,'Pr(>|z|)'])
                               beta <- signif(suma$coef[x,'coef']);#coeficient beta
                               HR <- signif(suma$coef[x,'exp(coef)']);#exp(beta)
                               HR.confint.lower <- signif(suma$conf.int[x,"lower .95"])
                               HR.confint.upper <- signif(suma$conf.int[x,"upper .95"])
                               CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                               res<-c(beta, HR, CI, p.value)
                               names(res)<-c("univ_beta", "univ_HR", "univ_95% CI for HR", "univ_p.value")
                               return(res)
                             })
                           })
    res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
    res$univ_beta <- sapply(res$univ_beta,function(x) as.numeric(as.character(x)))
    res$univ_HR <- sapply(res$univ_HR,function(x) as.numeric(as.character(x)))
    res$`univ_95% CI for HR` <- sapply(res$`univ_95% CI for HR`,as.character)
    res$`univ_p.value` <- sapply(res$`univ_p.value`,function(x) as.numeric(as.character(x)))
    return(res)
  }
  
  lasso_screen <- function(da){
    survInfo <- as.matrix(da[,1:2])
    colnames(survInfo) <- c("time","status")
    mat <- as.matrix(da[,3:ncol(da)])
    cvfit <- cv.glmnet(x=mat,y=survInfo,family="cox")
    cvfit_info <- coef(cvfit,s=cvfit$lambda.min)
    cvfit_act <- names(cvfit_info[which(cvfit_info!=0),])
    if(length(cvfit_act)>0){
      return(cvfit_act)
    }else{
      return(NULL)
    }
  }
  
  cross_screen <- function(da,split,seed){
    set.seed(seed)
    trainIndex <- createDataPartition(da$Dead, p=split, list=FALSE)
    data_train <- da[trainIndex,]
    res <- lasso_screen(data_train)
    return(res)
  }
  
  plot.shiny.km <- function(time, death, x, cut, title = "", ids = NULL,
                            subset = rep(TRUE, length(time)),
                            col = NULL,  xlab = NULL, ylab = NULL, hr.inverse = FALSE,
                            no.plot = FALSE, optimal.cut = FALSE, assign.val=NULL,dir=NULL,output=NULL,...) {
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
    km.group2 = surv_fit(Surv(time, death) ~ Group,data=km.data)
    #if(F){
    #km.group = ggsurv(km.group2,
    #                  main = title,
    #                  xlab = xlab, ylab = ylab,...) +
    #  ggplot2::coord_cartesian(ylim = c(0, 1))
    #    theme(plot.title = element_text(vjust = 0.5, colour = "black"))
    #}
    km.group = ggsurvplot(km.group2,data=km.data,
                          pval = TRUE, conf.int = F,
                          risk.table = TRUE, # Add risk table
                          risk.table.col = "strata", # Change risk table color by groups
                          linetype = "strata", # Change line type by groups
                          surv.median.line = "none", # Specify median survival
                          ggtheme = theme_bw(), # Change ggplot2 theme
                          palette = c("#E7B800", "#2E9FDF"))
    return(km.group)
    #return(invisible(km.data))
  }
  
  plot.shiny.km1 <- function(time, death, x, cut, title = "", ids = NULL,
                             subset = rep(TRUE, length(time)),
                             col = NULL,  xlab = NULL, ylab = NULL, hr.inverse = FALSE,
                             no.plot = FALSE, optimal.cut = FALSE, assign.val=NULL,dir=NULL,output=NULL,...) {
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
      #cut = quantile(x,cut/100)
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
    km.group2 = surv_fit(Surv(time, death) ~ Group,data=km.data)
    #if(F){
    #km.group = ggsurv(km.group2,
    #                  main = title,
    #                  xlab = xlab, ylab = ylab,...) +
    #  ggplot2::coord_cartesian(ylim = c(0, 1))
    #    theme(plot.title = element_text(vjust = 0.5, colour = "black"))
    #}
    km.group = ggsurvplot(km.group2,data=km.data,
                          pval = TRUE, conf.int = F,
                          risk.table = TRUE, # Add risk table
                          risk.table.col = "strata", # Change risk table color by groups
                          linetype = "strata", # Change line type by groups
                          surv.median.line = "none", # Specify median survival
                          ggtheme = theme_bw(), # Change ggplot2 theme
                          palette = c("#E7B800", "#2E9FDF"))
    #plot(km.group)
    return(km.group)
    #return(invisible(km.data))
  }
  
  multi_cox <- function(sig_features,data){
    #surv_da <- data[,c("Time","Dead")]
    #data <- data[,sig_features]
    #indx <- sapply(data, is.factor)
    #data[indx] <- lapply(data[indx],function(x) as.numeric(x))
    #data <- cbind(surv_da,data)
    sig_features_formula <- unlist(sapply(sig_features,function(x){paste0('`',x,'`')}))
    multi_formula <- as.formula(paste('Surv(Time, Dead) ~',paste(sig_features_formula,collapse='+')))
    multi_models <- coxph(multi_formula, data = data)
    suma <- summary(multi_models)
    #rownames(suma$coef) <- sig_features
    #rownames(suma$conf.int) <- sig_features
    multi_res <- sapply(rownames(suma$coef),function(x){
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
  
  forest_plot_multi <- function(res,prefix){
    result <- cbind(rownames(res),res)
    colnames(result)[1] <- "features"
    write.table(result,paste0(opt$outdir,'/',prefix,'_res.txt'),sep="\t",quote=F,row.names=F)
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
  
  forest_plot <- function(res,prefix){
    result <- cbind(rownames(res),res)
    colnames(result)[1] <- "features"
    write.table(result,paste0(opt$outdir,'/',opt$prefix,'.',prefix,'_res.txt'),sep="\t",quote=F,row.names=F)
    
    res$mean <- as.numeric(as.character(res[,2]))
    res$lower <- as.numeric(as.character(sapply(res[,3],function(x){unlist(strsplit(as.character(x),'-'))[1]})))
    res$upper <- as.numeric(as.character(sapply(res[,3],function(x){unlist(strsplit(as.character(x),'-'))[2]})))
    text <- rbind(colnames(res)[1:4],as.matrix(res[1:4]))
    rownames(text)[1] <- 'features'
    text <- cbind(rownames(text),text)
    plotHR(text,res,prefix)
  }
  
  risk_cal <- function(da,res){
    risk_score <- NULL
    for(i in 1:dim(da)[1]){
      score <- sum(res$multi_beta * da[i,rownames(res)])
      names(score) <- rownames(da)[i]
      risk_score <- c(risk_score,score)
    }
    return(risk_score)
  }
  
  plot_roc <- function(predict_obj,dir,prefix){
    pdf(paste0(opt$outdir,"/",prefix,"_risk_score_roc.pdf"))
    plot(predict_obj[[1]]$FP, predict_obj[[1]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="False positive rate",ylab="True positive rate",main="Time-dependent ROC curve",col="red",lwd=2,lty=1)
    lines(predict_obj[[2]]$FP, predict_obj[[2]]$TP,type="l",col="blue",lwd=2,lty=2)
    lines(predict_obj[[3]]$FP, predict_obj[[3]]$TP,type="l",col="black",lwd=2,lty=5)
    lg1 = paste0("AUC at 1 year","=",format(predict_obj[[1]]$AUC,dig=3))
    lg2 = paste0("AUC at 3 year","=",format(predict_obj[[2]]$AUC,dig=3))
    lg3 = paste0("AUC at 5 year","=",format(predict_obj[[3]]$AUC,dig=3))
    legend("bottomright",legend=c(lg1,lg2,lg3),bty="n",col=c("red","blue","black"),lty=c(1,2,5))
    abline(0,1,lwd=1)
    dev.off()
  }
  
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
  
  cross_validation <- function(data_trian,data_validation,Group3,Group4,WNT,gene_panel,seed){
    
    set.seed(seed)
    train_data <- data_train
    test_data <- data_validation
    train_multi_cox_res <- multi_cox(gene_panel,train_data)
    forest_plot_multi(train_multi_cox_res,"SHH_train_data_multicox")
    train_rs <- risk_cal(train_data,train_multi_cox_res)
    train_auc_1year_r1 <- survivalROC(Stime=train_data$Time,status=train_data$Dead,marker=train_rs,predict.time=1,span = 0.25*length(train_rs)^(-0.20))
    train_auc_3year_r1 <- survivalROC(Stime=train_data$Time,status=train_data$Dead,marker=train_rs,predict.time=3,span = 0.25*length(train_rs)^(-0.20))
    train_auc_5year_r1 <- survivalROC(Stime=train_data$Time,status=train_data$Dead,marker=train_rs,predict.time=5,span = 0.25*length(train_rs)^(-0.20))
    train_auc_obj_r1 <- list(train_auc_1year=train_auc_1year_r1,train_auc_3year=train_auc_3year_r1,train_auc_5year=train_auc_5year_r1)
    plot_roc(train_auc_obj_r1,opt$outdir,"SHH_train_data")
    pdf(paste0(opt$outdir,'/SHH_train_data_kmplot.pdf'))
    train_data_km <- plot.shiny.km(train_data$Time,train_data$Dead,train_rs,cut=50,optimal.cut = T,assign.val="Train_cutoff",title = "17 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
    print(train_data_km,newpage=F)
    dev.off()
    
    test_rs <- risk_cal(test_data,train_multi_cox_res)
    test_auc_1year_r1 <- survivalROC(Stime=test_data$Time,status=test_data$Dead,marker=test_rs,predict.time=1,span = 0.25*length(test_rs)^(-0.20))
    test_auc_3year_r1 <- survivalROC(Stime=test_data$Time,status=test_data$Dead,marker=test_rs,predict.time=3,span = 0.25*length(test_rs)^(-0.20))
    test_auc_5year_r1 <- survivalROC(Stime=test_data$Time,status=test_data$Dead,marker=test_rs,predict.time=5,span = 0.25*length(test_rs)^(-0.20))
    test_auc_obj_r1 <- list(test_auc_1year=test_auc_1year_r1,test_auc_3year=test_auc_3year_r1,test_auc_5year=test_auc_5year_r1)
    plot_roc(test_auc_obj_r1,opt$outdir,"SHH_test_data")
    pdf(paste0(opt$outdir,'/SHH_test_data_kmplot.pdf'))
    test_data_km <- plot.shiny.km1(test_data$Time,test_data$Dead,test_rs,cut=Train_cutoff,optimal.cut = F,assign.val=NULL,title = "17 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
    print(test_data_km,newpage=F)
    dev.off()
    
    Group3_rs <- risk_cal(Group3,train_multi_cox_res)
    Group3_auc_1year_r1 <- survivalROC(Stime=Group3$Time,status=Group3$Dead,marker=Group3_rs,predict.time=1,span = 0.25*length(Group3_rs)^(-0.20))
    Group3_auc_3year_r1 <- survivalROC(Stime=Group3$Time,status=Group3$Dead,marker=Group3_rs,predict.time=3,span = 0.25*length(Group3_rs)^(-0.20))
    Group3_auc_5year_r1 <- survivalROC(Stime=Group3$Time,status=Group3$Dead,marker=Group3_rs,predict.time=5,span = 0.25*length(Group3_rs)^(-0.20))
    Group3_auc_obj_r1 <- list(Group3_auc_1year=Group3_auc_1year_r1,Group3_auc_3year=Group3_auc_3year_r1,Group3_auc_5year=Group3_auc_5year_r1)
    plot_roc(Group3_auc_obj_r1,opt$outdir,"Group3")
    pdf(paste0(opt$outdir,'/Group3_kmplot.pdf'))
    Group3_km <- plot.shiny.km1(Group3$Time,Group3$Dead,Group3_rs,cut=Train_cutoff,optimal.cut = F,assign.val=NULL,title = "17 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
    print(Group3_km,newpage=F)
    dev.off()
    
    Group4_rs <- risk_cal(Group4,train_multi_cox_res)
    Group4_auc_1year_r1 <- survivalROC(Stime=Group4$Time,status=Group4$Dead,marker=Group4_rs,predict.time=1,span = 0.25*length(Group4_rs)^(-0.20))
    Group4_auc_3year_r1 <- survivalROC(Stime=Group4$Time,status=Group4$Dead,marker=Group4_rs,predict.time=3,span = 0.25*length(Group4_rs)^(-0.20))
    Group4_auc_5year_r1 <- survivalROC(Stime=Group4$Time,status=Group4$Dead,marker=Group4_rs,predict.time=5,span = 0.25*length(Group4_rs)^(-0.20))
    Group4_auc_obj_r1 <- list(Group4_auc_1year=Group4_auc_1year_r1,Group4_auc_3year=Group4_auc_3year_r1,Group4_auc_5year=Group4_auc_5year_r1)
    plot_roc(Group4_auc_obj_r1,opt$outdir,"Group4")
    pdf(paste0(opt$outdir,'/Group4_kmplot.pdf'))
    Group4_km <- plot.shiny.km1(Group4$Time,Group4$Dead,Group4_rs,cut=Train_cutoff,optimal.cut = F,assign.val=NULL,title = "17 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
    print(Group4_km,newpage=F)
    dev.off()
    
    WNT_rs <- risk_cal(WNT,train_multi_cox_res)
    WNT_auc_1year <- survivalROC(Stime=WNT$Time,status=WNT$Dead,marker=WNT_rs,predict.time=1,method="KM")
    WNT_auc_3year <- survivalROC(Stime=WNT$Time,status=WNT$Dead,marker=WNT_rs,predict.time=3,method="KM")
    WNT_auc_5year <- survivalROC(Stime=WNT$Time,status=WNT$Dead,marker=WNT_rs,predict.time=5,method="KM")
    #WNT_auc_1year_r1 <- survivalROC(Stime=WNT$OS.years,status=WNT$Dead,marker=WNT_rs,predict.time=1,span = 0.25*length(WNT_rs)^(-0.20))
    #WNT_auc_3year_r1 <- survivalROC(Stime=WNT$OS.years,status=WNT$Dead,marker=WNT_rs,predict.time=3,span = 0.25*length(WNT_rs)^(-0.20))
    #WNT_auc_5year_r1 <- survivalROC(Stime=WNT$OS.years,status=WNT$Dead,marker=WNT_rs,predict.time=5,span = 0.25*length(WNT_rs)^(-0.20))
    WNT_auc_obj <- list(WNT_auc_1year=WNT_auc_1year,WNT_auc_3year=WNT_auc_3year,WNT_auc_5year=WNT_auc_5year)
    #WNT_auc_obj_r1 <- list(WNT_auc_1year=WNT_auc_1year_r1,WNT_auc_3year=WNT_auc_3year_r1,WNT_auc_5year=WNT_auc_5year_r1)
    plot_roc(WNT_auc_obj,opt$outdir,"WNT")
    #plot_roc(WNT_auc_obj_r1,opt$outdir,"WNT")
    print(WNT_rs)
    pdf(paste0(opt$outdir,'/WNT_kmplot.pdf'))
    WNT_km <- plot.shiny.km1(WNT$Time,WNT$Dead,WNT_rs,cut=Train_cutoff,optimal.cut = F,assign.val=NULL,title = "17 mRNA based risk scoring prognostic model",xlab="time (years)",ylab="OS probability")
    print(WNT_km,newpage=F)
    dev.off()
    
    print(Train_cutoff)
  }
  
}

# -------------------create dir--------------------

if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir)
}
if(!dir.exists(paste0(opt$outdir,"/SHH_train_Cox"))){
  dir.create(paste0(opt$outdir,"/SHH_train_Cox"))
}
# --------------------run train--------------------

load("supplementary_append.Rdata")

set.seed(10086)

features <- diflist[,opt$genecol]

if(opt$ifuniv){
  univ_res <- univ_analysis(features,da=da)
  univ_res <- univ_res[!is.na(univ_res$univ_beta),]
  rownames(univ_res) <- gsub("`","",rownames(univ_res))
  univ_sig_res <- univ_res[univ_res$univ_p.value<=0.001,]
  #****************************************************************************
  newDF1 <- cbind(rownames(univ_res),univ_res)
  newDF1 <- as.data.frame(lapply(newDF1,unlist),stringsAsFactors = F,check.names = F)
  colnames(newDF1)[1] <- "features"
  write.table(newDF1,paste0(opt$outdir,"/SHH_train_Cox/",opt$train_prefix,".univ_res.txt"),row.names = F,col.names = T,quote=F,sep="\t")
  
  newDF <- dplyr::filter(newDF1,univ_p.value<=0.001)
  write.table(newDF,paste0(opt$outdir,"/SHH_train_Cox/",opt$train_prefix,".univ_sig_res.txt"),row.names = F,col.names = T,quote=F,sep="\t")
  #****************************************************************************
  univ_sig_genes <- rownames(univ_sig_res)
  features <- univ_sig_genes
}  # univ cox

do_da <- cbind(da[,1:2],da[,features])  # select univ cox sig gene for lesso

lasso_gene_list <- lapply(1:1000,FUN=function(x){print(x);res <- cross_screen(do_da,0.7,x);return(res)})  # 1000 times lesso with 70% sampling

lasso_gene_list <- lasso_gene_list[!sapply(lasso_gene_list,is.null)]

lasso_gene_freq <- table(unlist(lasso_gene_list))

lasso_gene_freq <- lasso_gene_freq[order(lasso_gene_freq,decreasing = T)]

lasso_gene_freq <- lasso_gene_freq[which(lasso_gene_freq>200)]  # select gene by freq

gene_panel <- names(lasso_gene_freq)

# ----------------select gene panel data----------------------

data_train <- dplyr::select(da,c("Time","Dead",gene_panel))

# ------------------------run validation----------------------------

cross_validation(data_train,data_validations[[1]],data_validations[[2]],data_validations[[3]],data_validations[[4]],gene_panel,1)

save(list=ls(),file="supplementary_append.Rdata")
