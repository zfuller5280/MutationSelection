#File paths will need to be changed to run locally!

#Packages
library(tidyr)
library(ggplot2)
library("HistogramTools")
library("ggbeeswarm")
require(gridExtra)
library(dplyr)
library(MASS)
library("tram")

#Read in simulation input data
updated_sim_input<-read.table("../data_files/sim.check.txt", sep="\t", header=T)
an_too_low<-updated_sim_input[(updated_sim_input$flag1=="AN_too_low"),]
updated_sim_input<-updated_sim_input[!(updated_sim_input$flag1=="AN_too_low"),]
mu_na<-updated_sim_input[is.na((updated_sim_input$mu_lof)),]
updated_sim_input<-updated_sim_input[!is.na((updated_sim_input$mu_lof)),]
updated_sim_input<-updated_sim_input[updated_sim_input$mu_lof>0,]
updated_n_occur<-data.frame(table(updated_sim_input$SYMBOL))
duplicated_genes<-updated_n_occur[updated_n_occur$Freq>1,]
updated_sim_input<-updated_sim_input[!(updated_sim_input$SYMBOL %in% updated_n_occur[updated_n_occur$Freq>1,]$Var1),]
head(updated_sim_input)
#Get raw autosome and X Chr input data
autosome_input<-subset(updated_sim_input, !(updated_sim_input$chromosome=="X"))
x_raw_input<-subset(updated_sim_input, (updated_sim_input$chromosome=="X"))

#Read in PAR genes
par_list<-read.table("../data_files/par_hs_input.data.7_7.tsv", header=T)
par_gene_list<-par_list$gene
x_chr_input<-subset(x_raw_input, !(x_raw_input$SYMBOL %in% par_gene_list))

#Determine which genes do not fit the model (i.e., too high of LoF freqs because of annotation)
neutral_auto_sims<-read.table("../data_files/automsome_neutral_summary.txt")
names(neutral_auto_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_auto_sims)
wrongModel_genes<-(neutral_auto_sims[neutral_auto_sims$obs_ptile>=.9, ])
neutral_x_sims<-read.table("../data_files/x_neutral_summary.txt")
names(neutral_x_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_x_sims)
wrongModel_x_genes<-(neutral_x_sims[neutral_x_sims$obs_ptile>=.9, ])
neutral_par_sims<-read.table("../data_files/par_neutral_summary.txt")
names(neutral_par_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_par_sims)
wrongModel_par_genes<-(neutral_par_sims[neutral_par_sims$obs_ptile>=.9, ])
valid_x_input<-subset(x_chr_input, !(x_chr_input$SYMBOL %in% wrongModel_x_genes$gene))
valid_x_input<-subset(valid_x_input, !(valid_x_input$SYMBOL %in% par_gene_list))

#X Chr Posteriors
x_chr_exp<-read.table("../data_files/x_chr_posteriors.8_10.txt")
names(x_chr_exp)<-c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")
x_chr_exp<-subset(x_chr_exp, !(x_chr_exp$Gene %in% par_gene_list))
x_chr_exp<-subset(x_chr_exp, !(x_chr_exp$Gene %in% wrongModel_x_genes$gene))
x_param_posts<-read.table("../data_files/xchr.posteriors.h_and_hs.8_12.tsv")
names(x_param_posts)<-c("Gene", "s_unscaled_hpd_low", "s_unscaled_hpd_high", "s_unscaled_hpd_map", "s_log10_ci_low", "s_log10_ci_high", "s_log10_map","hs_unscaled_hpd_low", "hs_unscaled_hpd_high", "hs_unscaled_hpd_map", "hs_log10_ci_low", "hs_log10_ci_high", "hs_log10_map")
x_param_posts<-subset(x_param_posts, x_param_posts$Gene %in% valid_x_input$SYMBOL)
x_param_posts<-merge(x_param_posts, x_chr_exp[,c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")], by="Gene")

#Autosomal/PAR genes
exp_posteriors<-read.table("../data_files/autosome_posteriors.10_11.txt") #17440 on 8_19
names(exp_posteriors)<-c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")
exp_posteriors<-exp_posteriors[!duplicated(exp_posteriors$Gene),] #Some mistakes in generating posterios, so a small number of genes are duplcited. Remove here
exp_par_posteriors<-subset(exp_posteriors, exp_posteriors$Gene %in% par_gene_list)
exp_posteriors<-subset(exp_posteriors, !(exp_posteriors$Gene %in% wrongModel_genes$gene))
exp_posteriors<-subset(exp_posteriors, !(exp_posteriors$Gene %in% par_gene_list))
subset(exp_posteriors, !(exp_posteriors$Gene %in% autosome_input$SYMBOL))
exp_posteriors<-subset(exp_posteriors, (exp_posteriors$Gene %in% autosome_input$SYMBOL))

#Label genes if they are on the NPX and have a Y homolog
npy_genes<-read.table("../data_files/NPY_genes.txt")
npy_input<-subset(x_chr_input, (x_chr_input$SYMBOL %in% npy_genes$V1))
npy_input$chrom<-24
npx_posts<-read.table("../data_files/npx_posteriors.9_13.txt")
names(npx_posts)<-c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")
head(npx_posts)
npx_xmodel<-subset(x_param_posts, Gene %in%  npx_posts$Gene)

#Label genes that are problematic by their status
problematic_df<-data.frame(gene=wrongModel_genes$gene, flag="wrong_model")
problematic_df<-rbind(problematic_df, data.frame(gene=remaining_autosome_input$SYMBOL, flag="not_finished"))
problematic_df<-rbind(problematic_df, data.frame(gene=an_too_low$SYMBOL, flag="AN_too_low"))
problematic_df<-rbind(problematic_df, data.frame(gene=mu_na$SYMBOL, flag="mu_NA"))
problematic_df<-rbind(problematic_df, data.frame(gene=duplicated_genes$Var1, flag="duplicated"))
problematic_df<-rbind(problematic_df, data.frame(gene=wrongModel_x_genes$gene, flag="wrong_model_x"))
problematic_df<-problematic_df[!duplicated(problematic_df$gene),]

#Get external data for Weghorn et al. estimates
weghorn<-read.table("../data_files/Supplementary_Table_1.txt",header=T)
head(weghorn)
merged_auto_weghorn<-merge(weghorn, exp_posteriors, by.x="Gene", by.y="Gene")
merged_x_weghorn<-merge(weghorn, x_chr_exp, by.x="Gene", by.y="Gene")
cor(merged_auto_weghorn$s_het_drift, 10^merged_auto_weghorn$log10_map, method="spearman")
cor(merged_auto_weghorn$s_het_drift, 10^merged_auto_weghorn$log10_map)^2
cor(merged_x_weghorn$s_het_drift, 10^merged_x_weghorn$log10_map, method="spearman")
cor(merged_x_weghorn$s_het_drift, 10^merged_x_weghorn$log10_map)

#Make supplemental figures comparing shet with hs
plot((merged_auto_weghorn$log10_map), log10(merged_auto_weghorn$s_het_drift),xlim=c(-6,0),ylim=c(-6,0))
auto_weg_comp<-ggplot(merged_auto_weghorn, aes(x=log10_map, y=log10(s_het_drift))) + geom_bin2d(bins=200)+theme_bw() + ylim(c(-6,0)) + xlim(c(-6,0)) + ggtitle("Autosomes") 
x_weg_comp<-ggplot(merged_x_weghorn, aes(x=log10_map, y=log10(s_het_drift))) + geom_bin2d(bins=100)+theme_bw() + ylim(c(-6,0)) + xlim(c(-6,0)) + ggtitle("X Chromosome")
grid.arrange(auto_weg_comp, x_weg_comp, nrow=1)

#Number of genes > -1e-3
length(exp_posteriors[exp_posteriors$log10_map>-3,1])/length(exp_posteriors[,1])
length(exp_posteriors[exp_posteriors$log10_ci_low>-3,1])/length(exp_posteriors[,1])

length(x_chr_exp[x_chr_exp$log10_map>-3,1])/length(x_chr_exp[,1])
length(x_chr_exp[x_chr_exp$log10_ci_low>-3,1])/length(x_chr_exp[,1])

#Get X-inactive status
tukiainen_escape<-read.csv("../data_files/tukiaianen_xci_status.csv")
tukiainen_escape<-tukiainen_escape[tukiainen_escape$Region=="nonPAR",]

tukiainen_escape<-tukiainen_escape[!(tukiainen_escape$Gene.name %in% npy_genes$V1),]
x_inact<-x_chr_exp[x_chr_exp$Gene %in% tukiainen_escape[tukiainen_escape$Reported.XCI.status=="Inactive","Gene.name"],]
x_escape<-x_chr_exp[x_chr_exp$Gene %in% tukiainen_escape[tukiainen_escape$Reported.XCI.status=="Escape","Gene.name"],]
x_variable<-x_chr_exp[x_chr_exp$Gene %in% tukiainen_escape[tukiainen_escape$Reported.XCI.status=="Variable","Gene.name"],]
test_p<-wilcox.test(x_inact$log10_map, x_escape$log10_map)

boxplot(x_inact$log10_map, x_escape$log10_map, names=c("Inactive","Escape"),main=paste0("pvalue: ",format(round(test_p$p.value, 4),nsmall=4)))
wilcox.test(x_variable$log10_map, x_inact$log10_map)

#Make plot of hs distribution for autosomes
legend.col <- function(col, lev){
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  xx <- rep(box.cx, each = 2)
  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}
dev.off()
auto_P<-ecdf(exp_posteriors$log10_map)
exp_posteriors$ecdf<-auto_P(exp_posteriors$log10_map)

#Make plot for autosomes
plot(exp_posteriors$log10_map,exp_posteriors$ecdf,xlim=c(-7,0),cex=0.25,xlab="log10(hs)",ylab="Cumulative Probability",main="Autosomes",ylim=c(0,1))
abline(v=-3,lty=3,col="red")
rbPal <- colorRampPalette(c('blue','cyan'))
exp_posteriors$ci_width<-exp_posteriors$log10_ci_high - exp_posteriors$log10_ci_low
ci_levels<-seq(6,0,-0.1)
exp_posteriors$Col <- rbPal(length(ci_levels))[as.numeric(cut(exp_posteriors$ci_width,breaks=ci_levels))]
segments(exp_posteriors$log10_ci_low,exp_posteriors$ecdf,exp_posteriors$log10_ci_high,lwd=.125,col=exp_posteriors$Col)
legend.col(col = rbPal(length(ci_levels)), lev = ci_levels)

#Distribution of hs for X chromosome
x_P<-ecdf(x_chr_exp$log10_map)
x_chr_exp$ecdf<-x_P(x_chr_exp$log10_map)
min(x_chr_exp$ecdf)
x_chr_exp$ci_width<-x_chr_exp$log10_ci_high - x_chr_exp$log10_ci_low
rbPal <- colorRampPalette(c('blue','cyan'))
x_chr_exp$Col <- rbPal(length(ci_levels))[as.numeric(cut(x_chr_exp$ci_width,breaks=ci_levels))]
#Make plot for X
plot(x_chr_exp$log10_map,x_chr_exp$ecdf,xlim=c(-7,0),cex=0.5,xlab="log10(hs)",ylab="Cumulative Probability",main="X Chromosome",ylim=c(0,1))
abline(v=-3,lty=3,col="red")
segments(x_chr_exp$log10_ci_low,x_chr_exp$ecdf,x_chr_exp$log10_ci_high,lwd=.35,col=x_chr_exp$Col)
legend.col(col = rbPal(length(ci_levels)), lev = ci_levels)

#Distribution of 0.5hs + 0.5s for X chromosome
x_param_posts$avg_sel_map<-(x_param_posts$log10_map * 0.5) + (x_param_posts$s_log10_map * 0.5)
x_param_posts$avg_sel_lo<-(x_param_posts$log10_ci_low * 0.5) + (x_param_posts$s_log10_ci_low * 0.5)
x_param_posts$avg_sel_hi<-(x_param_posts$log10_ci_high * 0.5) + (x_param_posts$s_log10_ci_high * 0.5)
x_P<-ecdf(x_param_posts$avg_sel_map)
x_param_posts$ecdf<-x_P(x_param_posts$avg_sel_map)
x_param_posts$ci_width<-x_param_posts$avg_sel_hi - x_param_posts$avg_sel_lo
rbPal <- colorRampPalette(c('blue','cyan'))
x_param_posts$Col <- rbPal(10)[as.numeric(cut(x_param_posts$ci_width,breaks = 10))]
#Make plot for X avg
plot(x_param_posts$avg_sel_map,x_param_posts$ecdf,xlim=c(-7,0),cex=0.5,xlab="log10((hs+s)/2)",ylab="Cumulative Probability",main="X Chromosome",ylim=c(0,1))
abline(v=-3,lty=3,col="red")
segments(x_param_posts$avg_sel_lo,x_param_posts$ecdf,x_param_posts$avg_sel_hi,lwd=.35,col=x_param_posts$Col)
legend.col(col = rbPal(10), lev = x_param_posts$ci_width)
length(x_param_posts[x_param_posts$avg_sel_lo>-3,1])/length(x_param_posts$Gene)

#Prepare jitter plots for different hs in different compartments
exp_par_posteriors$strong_sel<-ifelse(exp_par_posteriors$log10_ci_low>-3,"strong","other")
x_posteriors_exp_plot<-x_chr_exp[!(x_chr_exp$Gene %in% par_gene_list),]
x_posteriors_exp_plot$strong_sel<-ifelse(x_posteriors_exp_plot$log10_ci_low>-3,"strong","other")
x_posteriors_exp_plot_total<-x_posteriors_exp_plot
x_posteriors_exp_plot_total$Combined_XCI_status<-"XChrom"
exp_posteriors$strong_sel<-ifelse(exp_posteriors$log10_ci_low>-3,"strong","other")
exp_posteriors$Combined_XCI_status<-"Autosomes"
exp_par_posteriors$Combined_XCI_status<-"PAR"
x_posteriors_exp_compare<-rbind(x_posteriors_exp_plot_total[,colnames(exp_par_posteriors)], exp_par_posteriors)
sampled_subset_posteriors<-exp_posteriors[sample(nrow(exp_posteriors), 1000),]
x_posteriors_exp_compare<-rbind(sampled_subset_posteriors[,colnames(x_posteriors_exp_compare)], x_posteriors_exp_compare)
x_posteriors_exp_compare$Combined_XCI_status<-factor(x_posteriors_exp_compare$Combined_XCI_status, levels=c("Autosomes","XChrom","PAR"),ordered = T)
#specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
dev.off()
x_posteriors_exp_compare<-subset(x_posteriors_exp_compare, !(Combined_XCI_status=="PAR"))
#Make jitter plot and show
ggplot(x_posteriors_exp_compare, aes(x=Combined_XCI_status, y=log10_map)) +
  geom_boxplot( outlier.shape=NA) + 
  geom_quasirandom(aes(color=strong_sel),shape=21,bandwidth=0.4,varwidth=TRUE,method="pseudorandom")+
  xlab("Gene Status") + ylab("log10(hs) MAP") + scale_color_brewer(palette="Paired") + theme_bw()

#Not significant between 'escape' and 'autosomes'
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map,exp_posteriors$log10_map)
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='XChrom')$log10_map,subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map)
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='XChrom')$log10_map,subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='Autosomes')$log10_map)
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map,subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='Autosomes')$log10_map)
ks.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map,exp_posteriors$log10_map)

