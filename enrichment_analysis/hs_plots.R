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

dev.off()
#Read in data
#X Chr DDD
xchr_dnms<-read.csv("ddd_xchr_dnms.csv")

#Katsir et al. escape status
#katsir_escape<-read.csv("katsir_escape_status.csv")
#katsir_escape<-merge(katsir_escape,x_posteriors_exp_plot, by.x="Gene.Symbol", by.y="Gene")
#boxplot(log10_map~Annotation,data=katsir_escape)
#table(katsir_escape$Annotation)
#X Chr Posteriors
x_chr_exp<-read.table("/VOLUMES/Seagate/abc_hs/expanded_outfiles/x_chr_posteriors.4_28.txt")
names(x_chr_exp)<-c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")

#All gene input data
all_genes<-read.table("/Users/zachfuller/abc_docker/all.lof_summary.11_16.nfe_canonical.txt", header=T)
colnames(all_genes)[1] <- "Gene"
all_genes$mu_lof<-as.numeric(as.character(all_genes$mu_lof))
all_genes<-all_genes[!(is.na(all_genes$mu_lof)),]
all_genes<-all_genes[all_genes$mu_lof>0,]
all_x_genes<-subset(all_genes, all_genes$chromosome=="X")
neut_x_genes<-subset(all_x_genes, (all_x_genes$Gene %in% invalid_x_exp_gene_list))
x_chr_exp<-subset(x_chr_exp, !(x_chr_exp$Gene %in% invalid_x_exp_gene_list))
neut_x_genes$log10_map<--7
neut_x_genes$log10_ci_low<--7
neut_x_genes$log10_ci_high<--7
neut_x_genes[setdiff(names(x_chr_exp), names(neut_x_genes))] <- NA
gnomad_genes<-read.table("gnomad_mu_xi_par.txt",header=T)
gnomad_x_genes[gnomad_x_genes$mu_lof>0,][!(gnomad_x_genes[gnomad_x_genes$mu_lof>0,]$SYMBOL %in% x_posteriors_exp$Gene),]
table(gnomad_x_genes$SYMBOL)
table(gnomad_x_genes$Combined_XCI_status)
par_gene_list<-subset(gnomad_genes, !(is.na(gnomad_genes$par_status)))
neut_x_genes<-neut_x_genes[,colnames(neut_x_genes) %in% colnames(x_chr_exp)]
x_chr_exp<-rbind(x_chr_exp, neut_x_genes)
x_posteriors_exp<-merge(x_chr_exp, gnomad_x_genes, by.x="Gene", by.y="SYMBOL")
par_genes<-subset(x_posteriors_exp, !(is.na(x_posteriors_exp$par_status))) #These hs estiamtes are under the wrong (X Chr)
n_par_genes<-subset(par_genes, par_genes$log10_map<=-7)
x_chr_exp<-subset(x_chr_exp, !(x_chr_exp$Gene %in% par_gene_list$SYMBOL))
x_posteriors_exp<-subset(x_posteriors_exp, !(x_posteriors_exp$Gene %in% par_gene_list$SYMBOL))
estimated_x<-subset(x_posteriors_exp, (is.na(x_posteriors_exp$par_status)) & (x_posteriors_exp$log10_map>-7))
neutral_x<-subset(x_posteriors_exp, (is.na(x_posteriors_exp$par_status)) & (x_posteriors_exp$log10_map==-7))

#Autosomal genes
exp_posteriors<-read.table("/VOLUMES/Seagate/abc_hs/expanded_outfiles/autosome_posteriors.4_26.txt")
names(exp_posteriors)<-c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")
exp_par_posteriors<-subset(exp_posteriors, exp_posteriors$Gene %in% par_gene_list$SYMBOL)
neut_par_genes<-data.frame(read.table("/VOLUMES/Seagate/abc_hs/expanded_outfiles/par_neutral_gene_list.txt")[,1])
neut_par_genes<-subset(all_genes, all_genes$Gene %in% neut_par_genes[,1])
all_vars<-read.table("all_ddd_vars.csv",sep=",")
valid_sim_input<-read.table("all.lof.11_16.input.tsv")
names(valid_sim_input)<-c("gene","chromosome","mu_lof","AF_nfe","AN_nfe","NFE_k")
all_auto_genes<-subset(all_genes, !(all_genes$chromosome=="X"))
invalid_gene_list[!(invalid_gene_list %in% neut_auto_genes$Gene)]
neut_auto_genes<-subset(all_auto_genes, all_auto_genes$Gene %in% invalid_gene_list)
neut_auto_genes$log10_map<--7
neut_auto_genes$log10_ci_low<--7
neut_auto_genes$log10_ci_high<--7
neut_auto_genes[setdiff(names(exp_posteriors), names(neut_auto_genes))] <- NA
neut_auto_genes<-neut_auto_genes[,colnames(neut_auto_genes) %in% colnames(exp_posteriors)]
neut_par_genes$log10_map<--7
neut_par_genes$log10_ci_low<--7
neut_par_genes$log10_ci_high<--7
neut_par_genes[setdiff(names(exp_posteriors), names(neut_par_genes))] <- NA
neut_par_genes<-neut_par_genes[,colnames(neut_par_genes) %in% colnames(exp_posteriors)]
exp_posteriors<-rbind(exp_posteriors, neut_auto_genes)
exp_posteriors<-merge(exp_posteriors, gnomad_genes, by.x="Gene", by.y="SYMBOL")
exp_par_posteriors<-merge(exp_par_posteriors, gnomad_genes, by.x="Gene", by.y="SYMBOL")
neut_par_genes<-merge(neut_par_genes, gnomad_genes, by.x="Gene", by.y="SYMBOL")
weghorn<-read.table("~/Downloads/msz092_supplementary_data/Supplementary_Table_1.txt",header=T)
gnomad_mu<-read.table("~/Downloads/gnomad.v2.1.1.lof_metrics.by_transcript.txt",header=T,fill=T,sep='\t')
canonical_gnomad<-gnomad_mu[gnomad_mu$canonical=="true",]
exp_posteriors<-subset(exp_posteriors, !(exp_posteriors$Gene %in% exp_par_posteriors$Gene))
x_posteriors_exp<-rbind(x_posteriors_exp, exp_par_posteriors, neut_par_genes)

#EDS scores from Wang and Goldstein
eds<-read.csv("eds_scores.csv")
eds<-merge(canonical_gnomad, eds, by.x="gene_id", by.y="GeneSymbol")
eds_posteriors<-merge(exp_posteriors, eds, by.x="Gene", by.y="gene")

#DDD data from Samocha
denovo_west<-read.table("forZach_synexp_extended_denovoWEST_results.txt",sep="\t",header=T)
head(denovo_west)

#Imprinted genes
imprinted_genes<-read.table("imprinted_genes.txt")
invalid_imprinted_genes<-c("CPA4","CST1","GRB10","INPP5F","KIF25","LPAR6","MAGEL2","MAGI2","NTM","PRSS50","LINC01629","THEGL","UBE3A","UGT2B4","UTS2")
imprinted_genes<-subset(imprinted_genes, !(imprinted_genes$V1 %in% invalid_imprinted_genes))
imprinted_genes_posteriors<-subset(exp_posteriors, exp_posteriors$Gene %in% imprinted_genes[,1])
imprinted_genes_posteriors$Combined_XCI_status<-"Imprinted"
imprinted_genes_posteriors$strong_sel<-ifelse(imprinted_genes_posteriors$log10_ci_low>-3,"strong","other")

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
min(exp_posteriors$ecdf)
#Make plot
plot(exp_posteriors$log10_map,exp_posteriors$ecdf,xlim=c(-7,0),cex=0.25,xlab="log10(hs)",ylab="Cumulative Probability",main="Autosomes",ylim=c(0,1))
abline(v=-3,lty=3,col="red")
rbPal <- colorRampPalette(c('blue','cyan'))
exp_posteriors$ci_width<-exp_posteriors$log10_ci_high - exp_posteriors$log10_ci_low
exp_posteriors$Col <- rbPal(10)[as.numeric(cut(exp_posteriors$ci_width,breaks = 10))]
segments(exp_posteriors$log10_ci_low,exp_posteriors$ecdf,exp_posteriors$log10_ci_high,lwd=.125,col=exp_posteriors$Col)
legend.col(col = rbPal(10), lev = exp_posteriors$ci_width)

head(merged_obs_posteriors)
merged_obs_posteriors<-merge(exp_posteriors, canonical_gnomad, by.x="Gene", by.y="gene",all.x=T)
sum(as.numeric(as.character(merged_obs_posteriors$possible_lof)),na.rm=T)
merged_obs_posteriors[is.na(as.numeric(as.character(merged_obs_posteriors$possible_lof)))==TRUE,]
sum(as.numeric(as.character(merged_obs_posteriors[merged_obs_posteriors$log10_ci_low> -3.30103,"possible_lof"])),na.rm=T)

sum(merged_obs_posteriors[as.numeric(as.character(merged_obs_posteriors$s_het))> 0.9,"possible_lof"],na.rm=T)
plot(merged_obs_posteriors$log10_ci_low, merged_obs_posteriors$pLI)
length(exp_posteriors[exp_posteriors$log10_ci_low>-3,1])/(length(exp_posteriors[,1])+718)


#Distribution of hs for X chromosome
x_P<-ecdf(x_posteriors_exp$log10_map)
x_posteriors_exp$ecdf<-x_P(x_posteriors_exp$log10_map)
min(x_posteriors_exp$ecdf)
dev.off()
x_posteriors_exp$ci_width<-x_posteriors_exp$log10_ci_high - x_posteriors_exp$log10_ci_low
rbPal <- colorRampPalette(c('blue','cyan'))
x_posteriors_exp$Col <- rbPal(10)[as.numeric(cut(x_posteriors_exp$ci_width,breaks = 10))]
#Make plot
plot(x_posteriors_exp$log10_map,x_posteriors_exp$ecdf,xlim=c(-7,0),cex=0.5,xlab="log10(hs)",ylab="Cumulative Probability",main="X Chromosome",ylim=c(0,1))
abline(v=-3,lty=3,col="red")
segments(x_posteriors_exp$log10_ci_low,x_posteriors_exp$ecdf,x_posteriors_exp$log10_ci_high,lwd=.35,col=x_posteriors_exp$Col)
legend.col(col = rbPal(10), lev = x_posteriors_exp$ci_width)

neut_x_genes[neut_x_genes$Gene %in% escape_x$Gene,]

#Jitter plots for different hs in different compartments
x_par_posts<-x_posteriors_exp[!(is.na(x_posteriors_exp$par_status)),]
#write.table(x_par_posts[x_par_posts$log10_map>-7,names(x_chr_exp)],"5_3.par_posteriors.raw.txt",quote=F, row.names = F,col.names = F)
x_par_posts$strong_sel<-ifelse(x_par_posts$log10_ci_low>-3,"strong","other")
x_posteriors_exp_plot<-x_posteriors_exp[(is.na(x_posteriors_exp$par_status)),]
#new_escape_genes<-subset(katsir_escape, katsir_escape$Annotation=="confirmed")

#new_escape_genes<-katsir_escape$Gene.Symbol
#x_posteriors_exp_plot[which(x_posteriors_exp_plot$Combined_XCI_status=="escape"),]$Combined_XCI_status<-NA
table(x_posteriors_exp_plot$Combined_XCI_status)
#x_posteriors_exp_plot[x_posteriors_exp_plot$Gene %in% new_escape_genes,]$Combined_XCI_status<-"escape"
#escape_x<-subset(x_posteriors_exp_plot, (x_posteriors_exp_plot$Combined_XCI_status=="escape") & (x_posteriors_exp_plot$log10_map>-7))
#redo_escape<-subset(all_x_genes, all_x_genes$Gene %in% escape_x$Gene)
#head(redo_escape)
#write.table(redo_escape,"4_27_x_chr.escape.input.txt",row.names = F, col.names=F, quote=F, sep=" ")
x_posteriors_exp_plot$strong_sel<-ifelse(x_posteriors_exp_plot$log10_ci_low>-3,"strong","other")
x_posteriors_exp_plot_total<-x_posteriors_exp_plot
x_posteriors_exp_plot_total$Combined_XCI_status<-"XChrom"
exp_posteriors$strong_sel<-ifelse(exp_posteriors$log10_ci_low>-3,"strong","other")
exp_posteriors$Combined_XCI_status<-"Autosomes"
x_par_posts$Combined_XCI_status<-"PAR"
x_posteriors_exp_compare<-rbind(x_posteriors_exp_plot_total, x_par_posts)
#chr1_posteriors<-subset(exp_posteriors, exp_posteriors$chromosome=='22')
chr1_posteriors<-exp_posteriors[sample(nrow(exp_posteriors), 1000),]
x_posteriors_exp_compare<-rbind(x_posteriors_exp_compare, chr1_posteriors)
x_posteriors_exp_compare<-subset(x_posteriors_exp_compare, !(x_posteriors_exp_compare$Combined_XCI_status=="NA"))
x_posteriors_exp_compare$Combined_XCI_status<-factor(x_posteriors_exp_compare$Combined_XCI_status, levels=c("Autosomes","XChrom","PAR"),ordered = T)
#x_posteriors_exp_compare$Combined_XCI_status<-factor(x_posteriors_exp_compare$Combined_XCI_status, levels=c("Autosomes","PAR","escape","variable","inactive"),ordered = T)

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
auto_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="Autosomes") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="Autosomes"))$Gene)
#escape_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="escape") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="escape"))$Gene)
#variable_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="variable") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="variable"))$Gene)
#inactive_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="inactive") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="inactive"))$Gene)
par_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="PAR") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="PAR"))$Gene)
dev.off()
#xchr_neut_num<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="XChrom") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="XChrom"))$Gene)
#imprint_neut<-length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="Imprinted") & (x_posteriors_exp_compare$log10_map==-7))$Gene)/length(subset(x_posteriors_exp_compare, (x_posteriors_exp_compare$Combined_XCI_status=="Imprinted"))$Gene)
x_posteriors_exp_test<-merge(x_posteriors_exp_compare, canonical_gnomad, by.x="Gene", by.y="gene")
auto_par_colr<-Colr(log10_map ~ Combined_XCI_status + gene_length, data=x_posteriors_exp_test[x_posteriors_exp_test$Combined_XCI_status %in% c("Autosomes", "PAR"),])
summary(auto_par_colr)
auto_par_colr<-Colr(log10_map ~ Combined_XCI_status + gene_length, data=x_posteriors_exp_test[x_posteriors_exp_test$Combined_XCI_status %in% c("Autosomes", "XChrom"),])
plot(auto_par_colr)

mean(subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='PAR')$log10_map)
mean(subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='Autosomes')$log10_map)/mean(subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='PAR')$log10_map)

wilcox.test(subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='PAR')$possible_lof,subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='Autosomes')$possible_lof)
wilcox.test(subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='XChrom')$log10_map,subset(x_posteriors_exp_test,x_posteriors_exp_test$Combined_XCI_status=='Autosomes')$log10_map)

ggplot(x_posteriors_exp_compare, aes(x=Combined_XCI_status, y=log10_map)) +
  geom_boxplot( outlier.shape=NA) + 
  #geom_jitter(position = position_jitter(width = .05, height =0), shape=21, size=1.5,aes(color=strong_sel)) +
  geom_quasirandom(aes(color=strong_sel),shape=21,bandwidth=0.4,varwidth=TRUE,method="pseudorandom")+
  xlab("Gene Status") + ylab("log10(hs) MAP") + scale_color_brewer(palette="Paired") + theme_bw()
#  annotate("text",x=1, y=-7.25, label=specify_decimal(auto_neut, 3)) +
#  annotate("text",x=4, y=-7.25, label=specify_decimal(escape_neut, 3)) +
#  annotate("text",x=5, y=-7.25, label=specify_decimal(variable_neut, 3)) +
#  annotate("text",x=6, y=-7.25, label=specify_decimal(inactive_neut, 3))
#Get MWU p-values for testing between different compartments
#Not significant between 'escape' and 'autosomes'
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map,exp_posteriors$log10_map)
ks.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='PAR')$log10_map,exp_posteriors$log10_map)
#Significant between 'escape' and 'inactive'
wilcox.test(subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='escape')$log10_map,subset(x_posteriors_exp_compare,x_posteriors_exp_compare$Combined_XCI_status=='inactive')$log10_map)

gnomad_par<-canonical_gnomad[canonical_gnomad$gene %in% par_gene_list$SYMBOL,]
head(gnomad_par)
wilcox.test(gnomad_par$pLI, canonical_gnomad$pLI)

#DDD Autosomal Enrichment Plots
#Our hs point estimates
#DDD data from Samocha
denovo_west<-read.table("forZach_synexp_extended_denovoWEST_results.txt",sep="\t",header=T)
head(auto_dnm_counts)
denovo_west<-merge(denovo_west, canonical_gnomad, by.x="symbol", by.y="gene")
auto_dnm_counts<-merge(denovo_west, exp_posteriors, by.x="symbol", by.y="Gene",all.y=T)
auto_dnm_counts$mu_syn.y
#auto_dnm_counts<-merge(auto_dnm_counts, weghorn, by.x="symbol", by.y="Gene",all.x=T)

hs_bins<-c(-8,seq(-5,-.3,.3))
hs_bins<-c(-8, log10(c(1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1)))
hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
hs_labels<-paste(head(hs_bins,-1),tail(hs_bins,-1),sep="-")
auto_dnm_counts$hs_bin<-cut(auto_dnm_counts$log10_map, hs_bins, labels=hs_labels)
auto_dnm_counts[(is.na(auto_dnm_counts$lofexpected)),]
auto_dnm_counts[(is.na(auto_dnm_counts$lofexpected)),"lofexpected"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$synexp)),"synexp"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$lofcount)),"lofcount"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$synonymous_variant)),"synonymous_variant"]<-0
auto_dnm_counts$synexpected<-(auto_dnm_counts$mu_syn.y *2*31058)
hs_mode_agg<-aggregate(list(lof_expect=auto_dnm_counts$lofexpected,lof_obs=auto_dnm_counts$lofcount,syn_expect=auto_dnm_counts$synexpected,syn_obs=auto_dnm_counts$synonymous_variant),by=list(hs_bin=auto_dnm_counts$hs_bin),sum)
hs_mode_agg$lof_enrich<-hs_mode_agg$lof_obs/hs_mode_agg$lof_expect
hs_mode_agg$syn_enrich<-hs_mode_agg$syn_obs/hs_mode_agg$syn_expect
hs_mode_agg
iter_hs_enrich<-list()
iter_hs_syn_enrich<-list()
#Change the number range below for number of bootstraps (100 is for quick, 1000 takes longer but for final figures)
for(i in seq(1:1000)){
  iter_df<-list()
  count=1
  for(j in unique(auto_dnm_counts$hs_bin)){
    binned_hs<-auto_dnm_counts[auto_dnm_counts$hs_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_mode_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexpected,syn_obs=total_iter_df$synonymous_variant),by=list(hs_bin=total_iter_df$hs_bin),sum)
  iter_mode_agg$enrich<-iter_mode_agg$obs/iter_mode_agg$expect
  iter_mode_agg$syn_enrich<-iter_mode_agg$syn_obs/iter_mode_agg$syn_expect
  iter_hs_enrich[[i]]<-iter_mode_agg$enrich
  iter_hs_syn_enrich[[i]]<-iter_mode_agg$syn_enrich
}
iter_hs_enrich<-data.frame(iter_hs_enrich)
iter_hs_syn_enrich<-data.frame(iter_hs_syn_enrich)
names(iter_hs_enrich)<-seq(1:1000)
boot_hs_enrich<-apply(iter_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_enrich)<-c("lof_lo", "lof_hi")
boot_hs_syn_enrich<-apply(iter_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_syn_enrich)<-c("syn_lo", "syn_hi")
hs_enrich_df<-cbind(hs_mode_agg,t(boot_hs_enrich),t(boot_hs_syn_enrich))
colnames(hs_enrich_df)[1]<-"hsbin"
hs_enrich_df <- hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(hs_enrich_df)[2]<-"DNM type"
hs_plot<-(ggplot(hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
          + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
          + ylim(c(0,16)) + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
          + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
          + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3))
#Show plot for hs (point estimate) enrichment
hs_plot

#Enrichment for s_het
shet_bins<-seq(-3,0,0.5)
shet_labels<-paste(head(shet_bins,-1),tail(shet_bins,-1),sep="-")
auto_dnm_counts$shet_bin<-cut(log10(auto_dnm_counts$s_het_det), shet_bins, labels=shet_labels)
shet_bin_agg<-aggregate(list(lof_expect=auto_dnm_counts$lofexpected,lof_obs=auto_dnm_counts$lofcount,syn_expect=auto_dnm_counts$synexpected,syn_obs=auto_dnm_counts$synonymous_variant),by=list(shet_bin=auto_dnm_counts$shet_bin),sum)
shet_bin_agg$lof_enrich<-shet_bin_agg$lof_obs/shet_bin_agg$lof_expect
shet_bin_agg$syn_enrich<-shet_bin_agg$syn_obs/shet_bin_agg$syn_expect
iter_shet_enrich<-list()
iter_shet_syn_enrich<-list()
for(i in seq(1:100)){
  iter_df<-list()
  count=1
  for(j in unique(auto_dnm_counts$shet_bin)){
    binned_hs<-auto_dnm_counts[auto_dnm_counts$shet_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_shet_ci_low_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexpected,syn_obs=total_iter_df$synonymous_variant),by=list(shet_bin=total_iter_df$shet_bin),sum)
  iter_shet_ci_low_agg$enrich<-iter_shet_ci_low_agg$obs/iter_shet_ci_low_agg$expect
  iter_shet_ci_low_agg$syn_enrich<-iter_shet_ci_low_agg$syn_obs/iter_shet_ci_low_agg$syn_expect
  iter_shet_enrich[[i]]<-iter_shet_ci_low_agg$enrich
  iter_shet_syn_enrich[[i]]<-iter_shet_ci_low_agg$syn_enrich
}
iter_shet_enrich<-data.frame(iter_shet_enrich)
iter_shet_syn_enrich<-data.frame(iter_shet_syn_enrich)
names(iter_shet_enrich)<-seq(1:100)
boot_shet_enrich<-apply(iter_shet_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_shet_enrich)<-c("lof_lo", "lof_hi")
boot_shet_syn_enrich<-apply(iter_shet_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_shet_syn_enrich)<-c("syn_lo", "syn_hi")
shet_enrich_df<-cbind(shet_bin_agg,t(boot_shet_enrich),t(boot_shet_syn_enrich))
colnames(shet_enrich_df)[1]<-"shetbin"
shet_enrich_df <- shet_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(shet_enrich_df)[2]<-"DNM type"
shet_plot<-(ggplot(shet_enrich_df, aes(x=shetbin, y=(enrich), colour=`DNM type`)) + geom_point() 
            + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
            + ylim(c(0,12)) + xlab("shet bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
            + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
            + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3))
#Show plot for shet enrichment
shet_plot

#If you want to overlay s_het and hs plots
colnames(shet_enrich_df)[1]<-"hsbin"
shet_enrich_df$measure<-"shet"
hs_enrich_df$measure<-"hs"
enrich_df<-rbind(hs_enrich_df, shet_enrich_df)
shet_hs_plot<-(ggplot(enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
            + geom_line(aes(group=interaction(`DNM type`, measure),linetype=measure)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
            + ylim(c(0,12)) + xlab("shet bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray"))
shet_hs_plot

#pLI enrichment plot
pli_bins<-seq(0,1,0.1)
pli_labels<-paste(head(pli_bins,-1),tail(pli_bins,-1),sep="-")
auto_dnm_counts$pli_bin<-cut((auto_dnm_counts$pLI.x), pli_bins, labels=pli_labels)
pli_bin_agg<-aggregate(list(lof_expect=auto_dnm_counts$lofexpected,lof_obs=auto_dnm_counts$lofcount,syn_expect=auto_dnm_counts$synexpected,syn_obs=auto_dnm_counts$synonymous_variant),by=list(pli_bin=auto_dnm_counts$pli_bin),sum)
pli_bin_agg$lof_enrich<-pli_bin_agg$lof_obs/pli_bin_agg$lof_expect
pli_bin_agg$syn_enrich<-pli_bin_agg$syn_obs/pli_bin_agg$syn_expect

iter_pli_enrich<-list()
iter_pli_syn_enrich<-list()
for(i in seq(1:100)){
  iter_df<-list()
  count=1
  for(j in unique(auto_dnm_counts$pli_bin)){
    binned_hs<-auto_dnm_counts[auto_dnm_counts$pli_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_pli_ci_low_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexpected,syn_obs=total_iter_df$synonymous_variant),by=list(pli_bin=total_iter_df$pli_bin),sum)
  iter_pli_ci_low_agg$enrich<-iter_pli_ci_low_agg$obs/iter_pli_ci_low_agg$expect
  iter_pli_ci_low_agg$syn_enrich<-iter_pli_ci_low_agg$syn_obs/iter_pli_ci_low_agg$syn_expect
  iter_pli_enrich[[i]]<-iter_pli_ci_low_agg$enrich
  iter_pli_syn_enrich[[i]]<-iter_pli_ci_low_agg$syn_enrich
}
iter_pli_enrich<-data.frame(iter_pli_enrich)
iter_pli_syn_enrich<-data.frame(iter_pli_syn_enrich)
names(iter_pli_enrich)<-seq(1:100)
boot_pli_enrich<-apply(iter_pli_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_pli_enrich)<-c("lof_lo", "lof_hi")
boot_pli_syn_enrich<-apply(iter_pli_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_pli_syn_enrich)<-c("syn_lo", "syn_hi")
pli_enrich_df<-cbind(pli_bin_agg,t(boot_pli_enrich),t(boot_pli_syn_enrich))
colnames(pli_enrich_df)[1]<-"shetbin"
pli_enrich_df <- pli_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(pli_enrich_df)[2]<-"DNM type"
pli_plot<-(ggplot(pli_enrich_df, aes(x=shetbin, y=(enrich), colour=`DNM type`)) + geom_point() 
           + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
           + ylim(c(0,12)) + xlab("pLI bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
           + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
           + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3))
#Show plot for pLI enrichment
pli_plot

#LOEUF enrichment plot
luf_bin<-seq(0,2,0.25)

luf_labels<-paste(head(luf_bin,-1),tail(luf_bin,-1),sep="-")
auto_dnm_counts$luf_bin<-cut((auto_dnm_counts$oe_lof_upper), luf_bin, labels=luf_labels)
luf_bin_agg<-aggregate(list(lof_expect=auto_dnm_counts$lofexpected,lof_obs=auto_dnm_counts$lofcount,syn_expect=auto_dnm_counts$synexpected,syn_obs=auto_dnm_counts$synonymous_variant),by=list(luf_bin=auto_dnm_counts$luf_bin),sum)
luf_bin_agg$lof_enrich<-luf_bin_agg$lof_obs/luf_bin_agg$lof_expect
luf_bin_agg$syn_enrich<-luf_bin_agg$syn_obs/luf_bin_agg$syn_expect

iter_luf_enrich<-list()
iter_luf_syn_enrich<-list()
for(i in seq(1:100)){
  iter_df<-list()
  count=1
  for(j in unique(auto_dnm_counts$luf_bin)){
    binned_hs<-auto_dnm_counts[auto_dnm_counts$luf_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_luf_ci_low_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexpected,syn_obs=total_iter_df$synonymous_variant),by=list(luf_bin=total_iter_df$luf_bin),sum)
  iter_luf_ci_low_agg$enrich<-iter_luf_ci_low_agg$obs/iter_luf_ci_low_agg$expect
  iter_luf_ci_low_agg$syn_enrich<-iter_luf_ci_low_agg$syn_obs/iter_luf_ci_low_agg$syn_expect
  iter_luf_enrich[[i]]<-iter_luf_ci_low_agg$enrich
  iter_luf_syn_enrich[[i]]<-iter_luf_ci_low_agg$syn_enrich
}
iter_luf_enrich<-data.frame(iter_luf_enrich)
iter_luf_syn_enrich<-data.frame(iter_luf_syn_enrich)
names(iter_luf_enrich)<-seq(1:100)
boot_luf_enrich<-apply(iter_luf_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_luf_enrich)<-c("lof_lo", "lof_hi")
boot_luf_syn_enrich<-apply(iter_luf_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_luf_syn_enrich)<-c("syn_lo", "syn_hi")
luf_enrich_df<-cbind(luf_bin_agg,t(boot_luf_enrich),t(boot_luf_syn_enrich))
colnames(luf_enrich_df)[1]<-"shetbin"
luf_enrich_df <- luf_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(luf_enrich_df)[2]<-"DNM type"
luf_plot<-(ggplot(luf_enrich_df, aes(x=(shetbin), y=(enrich), colour=`DNM type`)) + geom_point() 
           + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
           + ylim(c(0,12)) + xlab("LOEUF bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
           + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
           + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3) + scale_x_discrete(limits=rev(levels(luf_enrich_df$shetbin))))
luf_plot

#Plot all 4 measures in one figure
grid.arrange(hs_plot, shet_plot,pli_plot,luf_plot)
dev.off()

#X chr enrichment
x_dnm_counts<-merge(denovo_west, x_posteriors_exp, by.x="symbol", by.y="Gene")

x_dnm_counts<-merge(x_dnm_counts, weghorn, by.x="symbol", by.y="Gene",all.x=T)
x_hs_bins<-c(-8,seq(-4,0,0.5))
x_hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
x_hs_labels<-paste(head(x_hs_bins,-1),tail(x_hs_bins,-1),sep="-")
x_dnm_counts$x_hs_bin<-cut(x_dnm_counts$log10_map, x_hs_bins, labels=x_hs_labels)
x_shet_bins<-seq(-3,0,0.5)
x_shet_labels<-paste(head(x_shet_bins,-1),tail(x_shet_bins,-1),sep="-")
x_dnm_counts$x_shet_bin<-cut(log10(x_dnm_counts$s_het_det), x_shet_bins, labels=x_shet_labels)
x_dnm_counts[(is.na(x_dnm_counts$lofexpected)),"lofexpected"]<-0
x_dnm_counts[(is.na(x_dnm_counts$synexp)),"synexp"]<-0
x_hs_ci_low_agg<-aggregate(list(lof_expect=x_dnm_counts$lofexpected,lof_obs=x_dnm_counts$lofcount,syn_expect=x_dnm_counts$synexp,syn_obs=x_dnm_counts$synonymous_variant),by=list(x_hs_bin=x_dnm_counts$x_hs_bin),sum)
x_hs_ci_low_agg$lof_enrich<-x_hs_ci_low_agg$lof_obs/x_hs_ci_low_agg$lof_expect
x_hs_ci_low_agg$syn_enrich<-x_hs_ci_low_agg$syn_obs/x_hs_ci_low_agg$syn_expect

iter_x_hs_enrich<-list()
iter_x_hs_syn_enrich<-list()
for(i in seq(1:1000)){
  iter_df<-list()
  count=1
  for(j in unique(x_dnm_counts$x_hs_bin)){
    binned_hs<-x_dnm_counts[x_dnm_counts$x_hs_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_x_hs_ci_low_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexp,syn_obs=total_iter_df$synonymous_variant),by=list(x_hs_bin=total_iter_df$x_hs_bin),sum)
  iter_x_hs_ci_low_agg$enrich<-iter_x_hs_ci_low_agg$obs/iter_x_hs_ci_low_agg$expect
  iter_x_hs_ci_low_agg$syn_enrich<-iter_x_hs_ci_low_agg$syn_obs/iter_x_hs_ci_low_agg$syn_expect
  iter_x_hs_enrich[[i]]<-iter_x_hs_ci_low_agg$enrich
  iter_x_hs_syn_enrich[[i]]<-iter_x_hs_ci_low_agg$syn_enrich
}

iter_x_hs_enrich<-data.frame(iter_x_hs_enrich)
iter_x_hs_syn_enrich<-data.frame(iter_x_hs_syn_enrich)
names(iter_x_hs_enrich)<-seq(1:1000)
boot_x_hs_enrich<-apply(iter_x_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_hs_enrich)<-c("lof_lo", "lof_hi")
boot_x_hs_syn_enrich<-apply(iter_x_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_hs_syn_enrich)<-c("syn_lo", "syn_hi")
x_hs_enrich_df<-cbind(x_hs_ci_low_agg,t(boot_x_hs_enrich),t(boot_x_hs_syn_enrich))
colnames(x_hs_enrich_df)[1]<-"hsbin"
x_hs_enrich_df <- x_hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
#x_hs_plot<-ggplot(x_hs_enrich_df, aes(x=hsbin, y=(enrich), colour=DNM_type)) + geom_point() + geom_line(aes(group=DNM_type)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() + ylim(c(-1,25))
colnames(x_hs_enrich_df)[2]<-"DNM type"
x_hs_plot<-(ggplot(x_hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
          + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
          + ylim(c(0,16)) + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
          + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
          + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3))
x_hs_plot

x_shet_bin_agg<-aggregate(list(lof_expect=x_dnm_counts$lofexpected,lof_obs=x_dnm_counts$lofcount,syn_expect=x_dnm_counts$synexp,syn_obs=x_dnm_counts$synonymous_variant),by=list(x_shet_bin=x_dnm_counts$x_shet_bin),sum)
x_shet_bin_agg$lof_enrich<-x_shet_bin_agg$lof_obs/x_shet_bin_agg$lof_expect
x_shet_bin_agg$syn_enrich<-x_shet_bin_agg$syn_obs/x_shet_bin_agg$syn_expect
iter_x_shet_enrich<-list()
iter_x_shet_syn_enrich<-list()
for(i in seq(1:1000)){
  iter_df<-list()
  count=1
  for(j in unique(x_dnm_counts$x_shet_bin)){
    binned_hs<-x_dnm_counts[x_dnm_counts$x_shet_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_x_shet_ci_low_agg<-aggregate(list(expect=total_iter_df$lofexpected,obs=total_iter_df$lofcount,syn_expect=total_iter_df$synexp,syn_obs=total_iter_df$synonymous_variant),by=list(x_shet_bin=total_iter_df$x_shet_bin),sum)
  iter_x_shet_ci_low_agg$enrich<-iter_x_shet_ci_low_agg$obs/iter_x_shet_ci_low_agg$expect
  iter_x_shet_ci_low_agg$syn_enrich<-iter_x_shet_ci_low_agg$syn_obs/iter_x_shet_ci_low_agg$syn_expect
  iter_x_shet_enrich[[i]]<-iter_x_shet_ci_low_agg$enrich
  iter_x_shet_syn_enrich[[i]]<-iter_x_shet_ci_low_agg$syn_enrich
}
iter_x_shet_enrich<-data.frame(iter_x_shet_enrich)
iter_x_shet_syn_enrich<-data.frame(iter_x_shet_syn_enrich)
names(iter_x_shet_enrich)<-seq(1:1000)
boot_x_shet_enrich<-apply(iter_x_shet_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_shet_enrich)<-c("lof_lo", "lof_hi")
boot_x_shet_syn_enrich<-apply(iter_x_shet_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_shet_syn_enrich)<-c("syn_lo", "syn_hi")
x_shet_enrich_df<-cbind(x_shet_bin_agg,t(boot_x_shet_enrich),t(boot_x_shet_syn_enrich))
colnames(x_shet_enrich_df)[1]<-"shetbin"
x_shet_enrich_df <- x_shet_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
x_shet_plot<-ggplot(x_shet_enrich_df, aes(x=shetbin, y=(enrich), colour=DNM_type)) + geom_point() + geom_line(aes(group=DNM_type)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() + ylim(c(-1,25))
x_shet_plot
require(gridExtra)
grid.arrange(x_hs_plot, x_shet_plot)

#If you want to overlay s_het and hs plots for the X
colnames(x_hs_enrich_df)[2]<-"DNM type"
colnames(x_shet_enrich_df)[2]<-"DNM type"
colnames(x_shet_enrich_df)[1]<-"hsbin"
x_shet_enrich_df$measure<-"shet"
x_hs_enrich_df$measure<-"hs"
x_enrich_df<-rbind(x_hs_enrich_df, x_shet_enrich_df)
x_shet_hs_plot<-(ggplot(x_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
               + geom_line(aes(group=interaction(`DNM type`, measure),linetype=measure)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
               + ylim(c(0,15)) + xlab("shet bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray"))
x_shet_hs_plot

#EDS scores from Wang and Goldstein
eds<-read.csv("eds_scores.csv")
eds<-merge(canonical_gnomad, eds, by.x="gene_id", by.y="GeneSymbol")
eds_posteriors<-merge(exp_posteriors, eds, by.x="Gene", by.y="gene")

eds_posteriors<-merge(eds_posteriors, weghorn, by.x="Gene", by.y="Gene")
eds_posteriors$hs_bin_decile<-ntile(eds_posteriors$log10_ci_low, 5)
eds_posteriors$pli_bin_decile<-ntile(eds_posteriors$pLI.y, 5)
eds_posteriors$loeuf_bin_decile<-ntile(-1*eds_posteriors$oe_lof_upper, 5)
eds_posteriors$shet_bin_decile<-ntile(eds_posteriors$s_het_det, 5)
#eds_posteriors_df<-eds_posteriors_df[complete.cases(eds_posteriors_df),]
eds_posteriors_df<-eds_posteriors %>% pivot_longer(cols = contains("_bin_decile"), names_to=c("measure",".value"), names_pattern= "(.*?)_(.*)")

ggplot(eds_posteriors_df, aes(x=bin_decile, y=EDS,fill=measure,colour=measure))+ 
  stat_summary(fun=mean, position=position_dodge(0.2), geom="point",size=3)  + 
  stat_summary(fun=mean, position=position_dodge(0.2), geom="line",size=0.5,linetype="dashed") + theme_bw() + 
  stat_summary(fun.data=mean_cl_boot, position=position_dodge(0.2), geom="errorbar",width=0.1) 



dev.off()

head(eds_posteriors_df[eds_posteriors_df$measure=="hs","bin_decile"])



