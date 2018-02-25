library("ggplot2")
BRCA.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
BRCA.pt <- read.table("/Users/shreyashekhar/Downloads/BRCA_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
A = intersect(gsub('[\\.]','-',colnames(BRCA.ge)),rownames(BRCA.pt))
pt = BRCA.pt[A,]
#Her2
BRCA1.ge = read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
her2=pt$lab_proc_her2_neu_immunohistochemistry_receptor_status
her2_df <- data.frame(her2)
row.names(her2_df)<-A
colnames(BRCA1.ge) = gsub('[\\.]','-',colnames(BRCA1.ge))

BRCA.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
BRCA.pt <- read.table("/Users/shreyashekhar/Downloads/BRCA_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
A = intersect(gsub('[\\.]','-',colnames(BRCA.ge)),rownames(BRCA.pt))
pt = BRCA.pt[A,]

#Her2
BRCA1.ge = read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
her2=pt$lab_proc_her2_neu_immunohistochemistry_receptor_status
her2_df <- data.frame(her2)
row.names(her2_df)<-A
colnames(BRCA1.ge) = gsub('[\\.]','-',colnames(BRCA1.ge))

her2_df[her2=="",] <- NA
#Get list of all row names that have value NA
her2_na =rownames(her2_df)[rowSums(is.na(her2_df)) > 0]
#Delete all rows that have NA
her2_df = her2_df[!rownames(her2_df) %in% her2_na, ]
#Delete all columns that have value NA
`%ni%` <- Negate(`%in%`)
BRCA1.ge = subset(BRCA1.ge,select = names(BRCA1.ge) %ni% her2_na)

t.test(BRCA.ge['ERBB2',her2=='Positive'],BRCA.ge['ERBB2',her2=='Negative'])

#ER
BRCA2.ge = read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
er = pt$breast_carcinoma_estrogen_receptor_status
er_df <- data.frame(er)
row.names(er_df)<-A
colnames(BRCA2.ge) = gsub('[\\.]','-',colnames(BRCA2.ge))

er_df[er=="",] <- NA
er_na =rownames(er_df)[rowSums(is.na(er_df)) > 0]
er_df = er_df[!rownames(er_df) %in% er_na, ]
`%ni%` <- Negate(`%in%`)
BRCA2.ge = subset(BRCA2.ge,select = names(BRCA2.ge) %ni% er_na)

t.test(BRCA.ge['ERBB2',er=='Positive'],BRCA.ge['ERBB2',er=='Negative'])

#PR
BRCA3.ge = read.table("/Users/shreyashekhar/Downloads/HiSeqV2_BRCA.txt",row.names = 1,header=TRUE)
pr = pt$breast_carcinoma_progesterone_receptor_status
pr_df <- data.frame(pr)
row.names(pr_df)<-A
colnames(BRCA3.ge) = gsub('[\\.]','-',colnames(BRCA3.ge))

pr_df[pr=="",] <- NA
pr_na =rownames(pr_df)[rowSums(is.na(pr_df)) > 0]
pr_df = pr_df[!rownames(pr_df) %in% pr_na, ]
`%ni%` <- Negate(`%in%`)
BRCA3.ge = subset(BRCA3.ge,select = names(BRCA3.ge) %ni% pr_na)


t.test(BRCA.ge['ERBB2',pr=='Positive'],BRCA.ge['ERBB2',pr=='Negative'])
s = apply(BRCA.ge,1,sd)

#LUNG
LUNG.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_LUNG.txt",row.names = 1,header=TRUE)
LUNG.pt <- read.table("/Users/shreyashekhar/Downloads/LUNG_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
LUNG.md <- read.table("/Users/shreyashekhar/Downloads/LUNG_mutation_matrix-2.csv",row.names = 1,header=TRUE,sep=',')
LKRAS.md <- read.table("/Users/shreyashekhar/Downloads/KRAS_Mutation_LUNG.csv",row.names = 1,header=TRUE,sep=',')
LALK.md <- read.table("/Users/shreyashekhar/Downloads/LUNG_ALK_Mutation.csv",row.names = 1,header=TRUE,sep=',')
LBRAF.md <-read.table("/Users/shreyashekhar/Downloads/BRAF_Mutation_LUNG.csv",row.names = 1,header=TRUE,sep=',')

B = intersect(gsub('[\\.]','-',colnames(LUNG.ge)),rownames(LUNG.pt))
pt1 = LUNG.pt[B,]
EGFR = pt1$EGFR

#LUNG_MUT
#EGFR
colnames(LUNG.ge) = gsub('[\\.]','-',colnames(LUNG.ge))
B1 = intersect(colnames(LUNG.ge),rownames(LUNG.md))
pt11 <- data.frame(LUNG.md[B1,])
row.names(pt11)<-B1
colnames(pt11)<- "EGFR"
EGFR1 = pt11$EGFR
LUNG1.ge <- data.frame(LUNG.ge[,B1])

#KRAS
B2 = intersect(colnames(LUNG.ge),rownames(LKRAS.md))
pt12 <- data.frame(LKRAS.md[B2,])
row.names(pt12)<-B2
colnames(pt12)<- "KRAS"
KRAS1 = pt12$KRAS
LUNG2.ge <- data.frame(LUNG.ge[,B2])

#ALK
B3= intersect(colnames(LUNG.ge),rownames(LALK.md))
pt13 <- data.frame(LALK.md[B3,])
row.names(pt13)<-B3
colnames(pt13)<- "ALK"
ALK = pt13$ALK
LUNG3.ge <- data.frame(LUNG.ge[,B3])

#BRAF
B4= intersect(colnames(LUNG.ge),rownames(LBRAF.md))
pt14 <- data.frame(LBRAF.md[B4,])
row.names(pt14)<-B4
colnames(pt14)<- "BRAF"
LBRAF = pt14$BRAF
LUNG4.ge <- data.frame(LUNG.ge[,B4])

#LAML 
LAML.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_LAML.txt",row.names = 1,header=TRUE)
LAML.pt <- read.table("/Users/shreyashekhar/Downloads/LAML_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
LKIT.md <-read.table("/Users/shreyashekhar/Downloads/KIT_Mutation_LAML.csv",row.names = 1,header=TRUE, sep=',')

colnames(LAML.ge) = gsub('[\\.]','-',colnames(LAML.ge))
C = intersect(colnames(LAML.ge),rownames(LAML.pt))
pt2 = LAML.pt[C,]
CD  = pt2$immunophenotype_cytochemistry_testing_result
#PH1 = pt2$cytogenetic_abnormality
PH2=pt2$cytogenetic_abnormality_other
PML1 = pt2$cytogenetic_abnormality
PML2 = pt2$cytogenetic_abnormality_other

#C-KIT Mutation
C1 = intersect(colnames(LAML.ge),rownames(LKIT.md))
pt21 <- data.frame(LKIT.md[C1,])
row.names(pt21)<-C1
colnames(pt21)<- "KIT"
KIT = pt21$KIT
LAML1.ge <- data.frame(LAML.ge[,C1])

LAML_fish = pt2$FISH_test_component
LAML_fish_df <-data.frame(LAML_fish)

# PML-RARA Clinical in Fish Test Component
pml_rara_df <-data.frame(LAML_fish)
colnames(pml_rara_df) = 'pml_rara'
pml_rara_df$y[grepl("PML-RAR", pml_rara_df$pml_rara)] = '1'
pml_rara_df$y[!grepl("PML-RAR", pml_rara_df$pml_rara)] = '0'

#Philadelphia chromosome Positive(BCR-ABL in Fish Test Component)
bcr_abl_df <-data.frame(LAML_fish)
colnames(bcr_abl_df) = 'bcr_abl'
bcr_abl_df$y[grepl("BCR-ABL", bcr_abl_df$bcr_abl)] = '1'
bcr_abl_df$y[!grepl("BCR-ABL", bcr_abl_df$bcr_abl)] = '0'

#Philadelphia chromosome Negative(BCR-ABL Missing in Fish Test Component)
bcr_abl_df$x[grepl("BCR-ABL", bcr_abl_df$bcr_abl)] = '0'
bcr_abl_df$x[!grepl("BCR-ABL", bcr_abl_df$bcr_abl)] = '1'

#COADREAD- Mutation and Gene expression data don't intersect from TCGA!
COADREAD.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_COADREAD.txt",row.names = 1,header=TRUE)
COADREAD.pt <- read.table("/Users/shreyashekhar/Downloads/COADREAD_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
CKRAS.md <- read.table("/Users/shreyashekhar/Downloads/KRAS_mutation_COADREAD.csv",row.names = 1,header=TRUE,sep=',')
BRAF1.md <- read.table("/Users/shreyashekhar/Downloads/BRAF_mutation_COADREAD.csv",row.names = 1,header=TRUE,sep=',')
DPYD.md <-read.table("/Users/shreyashekhar/Downloads/DPYD_mutation_COADREAD.csv",row.names = 1,header=TRUE,sep=',')
D = intersect(gsub('[\\.]','-',colnames(COADREAD.ge)),rownames(COADREAD.pt))
pt3 = COADREAD.pt[D,]
KRAS_CODON  = pt3$kras_mutation_codon
KRAS_FOUND = pt3$kras_mutation_found
BRAF = pt3$braf_gene_analysis_result

#COADREAD_MUT
#KRAS
colnames(COADREAD.ge) = gsub('[\\.]','-',colnames(COADREAD.ge))
D1 = intersect(colnames(COADREAD.ge),rownames(CKRAS.md))
pt31 <- data.frame(CKRAS.md[D1,])

#BRAF-COADREAD
D2 = intersect(colnames(COADREAD.ge),rownames(BRAF1.md))
pt32 <- data.frame(BRAF1.md[D2,])

#DYPD Mutation for DPD Deficiency
D3 = intersect(colnames(COADREAD.ge),rownames(DPYD.md))
pt33 <- data.frame(DPYD.md[D3,])
row.names(pt33)<-D3
colnames(pt33)<- "DPYD"
DPYD = pt33$DPYD
COADREAD1.ge <- data.frame(COADREAD.ge[,D3])

#SKCM
SKCM.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_SKCM.txt",row.names = 1,header=TRUE)
BRAF2.md <- read.table("/Users/shreyashekhar/Downloads/BRAF_Mutation_SKCM.csv",row.names = 1,header=TRUE,sep=',')

#BRAF-SKCM
colnames(SKCM.ge) = gsub('[\\.]','-',colnames(SKCM.ge))
E1 = intersect(colnames(SKCM.ge),rownames(BRAF2.md))
pt41 <- data.frame(BRAF2.md[E1,])
row.names(pt41)<-E1
colnames(pt41)<- "BRAF"
BRAF1 = pt41$BRAF
SKCM1.ge <- data.frame(SKCM.ge[,E1])

#THCA
THCA.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_THCA.txt",row.names = 1,header=TRUE)
TBRAF.md <- read.table("/Users/shreyashekhar/Downloads/BRAF_Mutation_THCA.csv",row.names = 1,header=TRUE,sep=',')

#BRAF-SKCM
colnames(THCA.ge) = gsub('[\\.]','-',colnames(THCA.ge))
F1 = intersect(colnames(THCA.ge),rownames(TBRAF.md))
pt51 <- data.frame(TBRAF.md[F1,])
row.names(pt51)<-F1
colnames(pt51)<- "BRAF"
TBRAF = pt51$BRAF
THCA1.ge <- data.frame(THCA.ge[,F1])


echofn <-function(v) {
  deparse(substitute(v))
}

library(plotROC)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(glmnet)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lay <- rbind(c(2,1),
               c(3,4),
               c(5,5))
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

cb <- function(x, print.row.names = TRUE) {
  clip <- pipe("pbcopy", "w");
  if (print.row.names == FALSE) {
    write.table(x,file=clip,sep="\t",row.names=FALSE,quote = FALSE);
  } else {
    write.table(x,file=clip,sep="\t",row.names=TRUE, col.names=NA);
    
  }
  close(clip)
}


#version-2
analyzeFeature2 = function(expr, feat, pos, neg, biostr, ctype) {
  s = apply(expr,1,sd)
  
  t.results = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
                                                    x[grepl(paste(neg, collapse = "|"), feat)]))
  
  t.statistic = unlist(lapply(t.results,FUN=function(x) x$statistic)) 
  p.value = unlist(lapply(t.results,FUN=function(x) x$p.value))
  
  p.adj = p.adjust(p.value, "fdr")
  
  # QQ plot
  pdf(paste0("qqplot_", biostr, "_", ctype, ".pdf"))
  o = -log10(sort(p.value,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(x= e, y= o, xlab = "Expected p-value", 
         ylab = "Observed p-value")
  dev.off()
  
  pl <- list()
  i = 1
  
  t.stats = data.frame(t.statistic)
  p.val = data.frame(p.value)
  p.adj_df = data.frame(p.adj)
  colnames(p.val) = 'p'
  p.val$t = t.stats$t.statistic
  p.val$padj = p.adj_df$p.adj
  p.val = p.val[order(p.val$padj),]
  
  temp_df = head(p.val)
  temp_df$r = rownames(temp_df)
  temp_df$z <- factor(temp_df$r, levels = temp_df$r[order(temp_df$padj)])
  pl[[i]] = ggplot(temp_df,aes(x=z,y=t, fill = r)) +geom_bar(stat='identity') +
    labs(x = "Genes", y = "T-statistic", fill = "Genes") +
    ggtitle("T-test Statistics") +
    theme(axis.text.x=element_blank(), axis.title=element_text(size=9),
          axis.ticks.x=element_blank(),
          plot.title = element_text(size=10),
          legend.position="none")
  #ggsave(paste("ggplot_ttest_", biostr, ".pdf"))
  i = i+1
  
  top_genes = c(temp_df$r)
  
  #pdf(paste("roc_plots_", biostr, ".pdf"))
  
  auc_list <- as.numeric(list())
  auc_glist <- as.numeric(list())
  tmp_df <- data.frame(feat)
  colnames(tmp_df) = 'level'
  tmp_df$level0[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 'Positive'
  tmp_df$level0[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 'Negative'
  names(tmp_df) [ncol(tmp_df)] <- 'lvl0'
  
  samples <- data.frame(tmp_df$lvl0)
  sample_counts <- table(samples, exclude=NULL)
  print(sample_counts)
  samp <- data.frame(sample_counts)
  group.colors <- c("Positive" = "aquamarine2", "Negative" = "deepskyblue1")
  
  pl[[i]] = ggplot(samp, aes(samples, Freq, fill = samples))+ geom_bar(stat="identity")+
    ggtitle(paste(biostr, "- samples")) +
    labs(y = "Num. of samples") +
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.title=element_text(size=9),
          legend.position="bottom",
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(size=10)) + 
          scale_fill_manual(values=group.colors, na.value = "yellow") 
          
  i=i+1
  
  j=1
  
    tmp_df$level1[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 1
  tmp_df$level1[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 0
  
  names(tmp_df) [ncol(tmp_df)] <- 'lvl'
  
  #cb(tmp_df$lvl)
  
  for (g in top_genes) {
    
    tmp_df$gene = as.numeric(expr[g,])
    names(tmp_df) [ncol(tmp_df)] <- g
    tmp_df$tgene = as.numeric(expr[g,])
    
    
    plt <- ggplot(tmp_df, aes(d=lvl, m= tgene)) + geom_roc() 
    auc_list[j] <- round(calc_auc(plt)$AUC, 4)
    auc_glist[g] <- round(calc_auc(plt)$AUC, 4)
    if (auc_list[j] < 0.5) {
        auc_list[j] = 1- auc_list[j]
  
    }
    print(auc_list[j])
    j = j+1
    
  }
  
  longtest <- melt_roc(tmp_df, "lvl", c(temp_df$r))
  for (g in top_genes) {
     if (auc_glist[g] < 0.5) {
       longtest <- within(longtest, D[D == 1 & name %in% g ] <- 2)
       longtest <- within(longtest, D[D == 0 & name %in% g ] <- 1)
       longtest <- within(longtest, D[D == 2 & name %in% g ] <- 0)
     }
    print(auc_glist[g])
  }  
  pl[[i]] = ggplot(longtest, aes(d = D, m = M, color = name)) + 
    #ggplot(longtest, aes(x = D, y = M, label = c, colour = name)) +
    geom_roc(labels = FALSE) + style_roc(theme = theme_gray) +
    theme(axis.text = element_text(colour = "blue"), 
          axis.title=element_text(size=9),
          legend.position="none",
          plot.title = element_text(size=10)) +
    ggtitle(paste(biostr, "- ROC Curves")) + 
    scale_x_continuous(breaks = seq(0, 1, by = .2)) 
  #ggsave(paste("rocCurves_", biostr, ".pdf"))
  i=i+1
  
  
  
  auc_df = data.frame(auc_list)
  #print (auc_list)
  colnames(auc_df) = 'x'
  auc_df$y = c(temp_df$r)
  auc_df = auc_df[order(-auc_df$x),]
  auc_df$z <- factor(auc_df$y, levels = auc_df$y[order(-auc_df$x)])
  print(auc_df)
  
  pl[[i]] = ggplot(auc_df, aes(z, x, fill = z))+ geom_bar(stat="identity")+
    ggtitle(paste(biostr, " AUC values")) + labs(x = "Genes", y = "AUC", fill = "Genes") +
    theme(axis.text.x=element_blank(), #axis.title=element_text(size=9), 
        axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom", 
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(size=10))
  #ggsave(paste("auc_barplot_", biostr, ".pdf"))
  
  i=i+1
  
  temp3_df <- data.frame(head(auc_df, 3))
  print(temp3_df)
  btest <- melt(tmp_df,id.vars = c("lvl0","lvl"), measure.vars = c(temp3_df$y))
  pl[[i]] <- ggplot(btest, aes(x=factor(lvl0),y=value, fill = factor(lvl0))) +
    geom_boxplot(outlier.size=0,show.legend = F,lwd=0.1) + labs(title="Box plots for Top 3 Genes") +facet_wrap(~variable) +
    theme(axis.text.x=element_blank(),
          #axis.title=element_text(size=9),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title = element_text(size=6), 
          legend.text = element_text(size=4),
          plot.title = element_text(size=10),
          legend.position="bottom") +
          scale_fill_manual(values=group.colors, na.value="yellow")
  i=i+1
  
  lay <- rbind(c(2,1),
               c(3,4),
               c(5,5))
  ml = grid.arrange(grobs = pl, layout_matrix = lay, top= paste(biostr, " in ", ctype)) 
                    #legend = legend)
  
  
  ggsave(paste0("plots_", biostr, "_in_", ctype, ".pdf"), plot=ml, width = 6, height = 8, units = "in")
  
} 


#Machine Learning for each Biomarker

mlasso = function(expr, feat, feat_str){
  data  <-expr
  
  data <- data.frame(t(data))
  
  #making data set
  feat_bin <- as.numeric( feat == "1" )
  
  data$feat_str <- feat_bin
  feat_target <- data$feat_str
  
  #training + testing sets
  feat_trsize<- floor(0.75 * nrow(data))
  feat_tr <- sample(seq_len(nrow(data)), size = feat_trsize)
  feat_train <- data[feat_tr, ]
  feat_test <- data[-feat_tr, ]
  
  #training model
  feat_responsecol <- which(colnames(feat_train) == "feat_str")
  feat_trainxdm <- data.matrix(feat_train[, -feat_responsecol])
  feat_lasso_model <- cv.glmnet(x = feat_trainxdm, y = feat_train$feat_str, alpha = 1)
  
  #test set predictions
  feat_testxdm <- data.matrix(feat_test[, -feat_responsecol])
  
  feat_predictions_lasso <- predict(feat_lasso_model, newx = feat_testxdm, type = "response", 
                                    s = "lambda.min")[, 1]
  
  #Model performance (AUC)
  print(auc(feat_test$feat_str, feat_predictions_lasso))
  
}


analyzeFeature2(BRCA1.ge,her2_df,c('Positive','Equivocal'),c('Negative','Indeterminate'), echofn(Her2), echofn(BRCA))
analyzeFeature2(BRCA2.ge,er_df,c('Positive','Equivocal'),c('Negative','Indeterminate'), echofn(ER), echofn(BRCA))
analyzeFeature2(BRCA3.ge,pr_df,c('Positive','Equivocal'),c('Negative','Indeterminate'), echofn(PR), echofn(BRCA))


#EGFR in LUNG
pos_vec <- seq(746,753,1)
pos_v = as.vector(apply(as.matrix(pos_vec), 2, function(x) paste0(x,"del")))
p_vec = append(pos_v, "L858R")
n_vec = as.vector(unique(EGFR[!grepl(paste(p_vec, collapse="|"), EGFR)])) 
n_vec = unlist(lapply(n_vec,FUN=function(x) strsplit(x, split=", ")))
analyzeFeature2(LUNG.ge,EGFR, c(p_vec), c(n_vec), echofn(EGFR), echofn(LUNG))

#Test EGFR w/mutation data

analyzeFeature2(LUNG1.ge,EGFR1, c("1"), c("0"), echofn(EGFR), echofn(LUNG-Mutation))

#KRAS Mutation in LUNG
analyzeFeature2(LUNG2.ge,KRAS1, c("1"), c("0"), echofn(KRAS), echofn(LUNG-Mutation))

#ALK Mutation in LUNG
analyzeFeature2(LUNG3.ge,ALK, c("1"), c("0"), echofn(ALK), echofn(LUNG-Mutation))

#BRAF Mutation in LUNG
analyzeFeature2(LUNG4.ge,LBRAF, c("1"), c("0"), echofn(BRAF), echofn(LUNG-Mutation))


#LAML- CD20

cd_df <- data.frame(CD)
colnames(cd_df) = 'cd20'

cd_df$y[grepl("CD20", cd_df$cd20)] = 'Positive'
cd_df$y[grepl("CD20 Negative", cd_df$cd20)] = 'Negative'
analyzeFeature2(LAML.ge, cd_df$y, c('Positive'), c('Negative'), echofn(CD20), echofn(LAML))


#COADREAD - KRAS
krasc_df <- data.frame(KRAS_CODON)
colnames(krasc_df) = 'codon'
krasf_df <- data.frame(KRAS_FOUND)
colnames(krasf_df) = 'found'

p_vec <- c("12", "13")
krasc_df$y[grepl(paste(p_vec, collapse="|"), krasc_df$codon)] = 'Positive'
krasc_df$y[grepl("NO", krasf_df$found)] = 'Negative'
analyzeFeature2(COADREAD.ge, krasc_df$y, c('Positive'), c('Negative'), echofn(KRAS), echofn(COADREAD))


#COADREAD - BRAF
analyzeFeature2(COADREAD.ge, BRAF, c('Abnormal'), c('Normal'), echofn(BRAF), echofn(COADREAD))

#LAML- Philadelphia
ph2_df <- data.frame(PH2)
colnames(ph2_df) = 'ph2'
ph2_df$y[grepl("9;22", ph2_df$ph2)] = 'Positive'
ph2_df$y[(ph2_df$ph2 == "")] = 0
ph2_df$y[(!grepl("9;22", ph2_df$ph2)) & (!grepl("0", ph2_df$y))] = 'Negative'
ph2_df$y[grepl("0", ph2_df$y)] = NA
analyzeFeature2(LAML.ge, ph2_df$y , c('Positive'), c('Negative'), echofn(PH), echofn(LAML))

#LAML- PML/RARa
pml_df <- data.frame(PML1)
colnames(pml_df) = 'pml1'
pml_df$pml2 = PML2
pml_df$y[grepl("15;17", pml_df$pml1)] = 'Positive'
pml_df$y[grepl("15;17", pml_df$pml2)] = 'Positive'
pml_df$y[(pml_df$pml1 == "") & (pml_df$pml2 == "")] = 0
pml_df$y[(!grepl("15;17", pml_df$pml1))& (!grepl("15;17", pml_df$pml2)) & (!grepl("0", pml_df$y))] = 'Negative'
pml_df$y[grepl("0", pml_df$y)] = NA

analyzeFeature2(LAML.ge, pml_df$y , c('Positive'), c('Negative'), echofn(PML_RARa), echofn(LAML))

# LAML PML-RARA in FISH_TEST
analyzeFeature2(LAML.ge, pml_rara_df$y, c("1"), c("0"), echofn(PML-RARA), echofn(LAML))

#Philadelphia chromosome Positive LAML BCR-ABL in FISH_TEST
analyzeFeature2(LAML.ge, bcr_abl_df$y, c("1"), c("0"), echofn(BCR-ABL), echofn(LAML))

#Philadelphia chromosome Negative BCR-ABL Missing in FISH_TEST 
analyzeFeature2(LAML.ge, bcr_abl_df$x, c("1"), c("0"), echofn(BCR-ABL), echofn(LAML))

#LAML KIT Mutation
analyzeFeature2(LAML1.ge, KIT , c("1"), c("0"), echofn(KIT), echofn(LAML))

#SKCM- BRAF Mutation
analyzeFeature2(SKCM1.ge,BRAF1, c("1"), c("0"), echofn(BRAF), echofn(SKCM-Mutation))

#THCA BRAF Mutation
analyzeFeature2(THCA1.ge,TBRAF, c("1"), c("0"), echofn(BRAF), echofn(THCA-Mutation))
