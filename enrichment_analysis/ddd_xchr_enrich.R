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

#Determine which genes do not fit the model (i.e., too high of LoF freqs because of annotation)
neutral_auto_sims<-read.table("/Volumes/Seagate/abc_hs/automsome_neutral_summary.txt")
names(neutral_auto_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_auto_sims)
wrongModel_genes<-(neutral_auto_sims[neutral_auto_sims$obs_ptile>=.9, ])
neutral_x_sims<-read.table("/Volumes/Seagate/abc_hs/x_neutral_summary.txt")
names(neutral_x_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_x_sims)
wrongModel_x_genes<-(neutral_x_sims[neutral_x_sims$obs_ptile>=.9, ])
neutral_par_sims<-read.table("/Volumes/Seagate/abc_hs/par_neutral_summary.txt")
names(neutral_par_sims)<-c("gene", "mu_lof", "AF_nfe", "AN_nfe", "NFE_k", "mean_sim_AF", "obs_ptile", "mean_cnt", "exact_count", "median", "lower_mid", "upper_mid", "count_tol_2", "count_tol_5", "count_tol_10")
head(neutral_par_sims)
wrongModel_par_genes<-(neutral_par_sims[neutral_par_sims$obs_ptile>=.9, ])
valid_x_input<-subset(x_chr_input, !(x_chr_input$SYMBOL %in% wrongModel_x_genes$gene))
valid_x_input<-subset(valid_x_input, !(valid_x_input$SYMBOL %in% par_gene_list))

#Read in X Chr data
x_param_posts<-read.table("/VOLUMES/Seagate/abc_hs/expanded_outfiles/xchr.posteriors.h_and_hs.8_12.tsv")
names(x_param_posts)<-c("Gene", "s_unscaled_hpd_low", "s_unscaled_hpd_high", "s_unscaled_hpd_map", "s_log10_ci_low", "s_log10_ci_high", "s_log10_map","hs_unscaled_hpd_low", "hs_unscaled_hpd_high", "hs_unscaled_hpd_map", "hs_log10_ci_low", "hs_log10_ci_high", "hs_log10_map")
x_param_posts<-subset(x_param_posts, x_param_posts$Gene %in% valid_x_input$SYMBOL)
x_param_posts<-merge(x_param_posts, x_chr_exp[,c("Gene","unscaled_hpd_low","unscaled_hpd_high","unscaled_hpd_map","log10_ci_low","log10_ci_high","log10_map","2.5th","97.5th","mean","median")], by="Gene")
cat(subset(valid_x_input, !(valid_x_input$SYMBOL %in% x_param_posts$Gene))$SYMBOL,sep="\n")

ddd_vars<-read.table("~/Downloads/DDD_RUMC_GDX_denovos_cadd_shet_wweights_2020_01_17.txt",sep="\t",header=T)
ddd_ids<-read.table("fordist_joint_dnm_ID_sex_2019_08_30.txt",header=T)
head(ddd_vars)
male_ids<-subset(ddd_ids, ddd_ids$sex=="M")$id
female_ids<-subset(ddd_ids, ddd_ids$sex=="F")$id
male_samp_size<-length(male_ids)
female_samp_size<-length(female_ids)
ddd_vars<-merge(ddd_vars, ddd_ids, by="id")

ddd_vars[ddd_vars$symbol=="MECP2",]

male_ddd<-subset(ddd_vars, ddd_vars$sex=="M")
female_ddd<-subset(ddd_vars, ddd_vars$sex=="F")

male_agg<-male_ddd %>% group_by(symbol, consequence) %>% summarize(count=n()) %>% spread(consequence, count) %>% as.data.frame()
male_agg[is.na(male_agg)]<-0
male_agg$lof_obs<-male_agg$stop_gained + male_agg$splice_acceptor_variant + male_agg$splice_donor_variant + male_agg$frameshift_variant 
head(male_agg)
female_agg<-female_ddd %>% group_by(symbol, consequence) %>% summarize(count=n()) %>% spread(consequence, count) %>% as.data.frame()
female_agg[is.na(female_agg)]<-0
female_agg$lof_obs<-female_agg$stop_gained + female_agg$splice_acceptor_variant + female_agg$splice_donor_variant + female_agg$frameshift_variant 

male_x_agg<-subset(male_agg, male_agg$symbol %in% x_param_posts$Gene)
female_x_agg<-subset(female_agg, female_agg$symbol %in% x_param_posts$Gene)

#X chr enrichment
alpha=3.5
x_male_counts<-merge(male_x_agg, x_param_posts, by.x="symbol", by.y="Gene",all.y=T)
x_male_counts<-x_male_counts[!(x_male_counts$symbol %in% par_param_posts$Gene),]
x_male_counts[(is.na(x_male_counts$lof_obs)),"lof_obs"]<-0
x_male_counts[(is.na(x_male_counts$synonymous_variant)),"synonymous_variant"]<-0
head(x_male_counts)
#x_dnm_counts<-merge(x_dnm_counts, weghorn, by.x="symbol", by.y="Gene",all.x=T)
x_hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
x_hs_labels<-paste(head(x_hs_bins,-1),tail(x_hs_bins,-1),sep="-")
x_male_counts$x_hs_bin<-cut(x_male_counts$s_log10_map, x_hs_bins, labels=x_hs_labels)
x_male_counts<-merge(x_male_counts, denovo_west, by.x="symbol", by.y="symbol",all.x=T)
gnomad_x_genes<-canonical_gnomad[canonical_gnomad$chromosome=="X",]
x_male_counts<-merge(x_male_counts, gnomad_x_genes, by.x="symbol", by.y="gene",all.x=T)
x_male_counts$lof_rate<-ifelse(is.na(x_male_counts$p_lof),x_male_counts$mu_lof,x_male_counts$p_lof)
x_male_counts$syn_rate<-ifelse(is.na(x_male_counts$p_syn),x_male_counts$mu_lof,x_male_counts$p_syn)
x_male_counts$lofexpected<-3* x_male_counts$lof_rate * (alpha/(2+alpha)) * male_samp_size
x_male_counts$synexpected<-3* x_male_counts$mu_syn * (alpha/(2+alpha)) * male_samp_size
x_male_counts[x_male_counts$symbol %in% par_param_posts$Gene,]$lofexpected<-x_male_counts[x_male_counts$symbol %in% par_param_posts$Gene,]$p_lof*male_samp_size
x_male_counts[x_male_counts$symbol %in% par_param_posts$Gene,]$synexpected<-x_male_counts[x_male_counts$symbol %in% par_param_posts$Gene,]$p_syn*male_samp_size

x_male_counts[(is.na(x_male_counts$lofexpected)),"lofexpected"]<-0
x_male_counts[(is.na(x_male_counts$synexpected)),"synexpected"]<-0
x_male_counts[x_male_counts$symbol %in% par_param_posts$Gene,]
x_hs_male_agg<-aggregate(list(lof_expect=x_male_counts$lofexpected,lof_obs=x_male_counts$lof_obs,syn_expect=x_male_counts$synexpected,syn_obs=x_male_counts$synonymous_variant.x),by=list(x_hs_bin=x_male_counts$x_hs_bin),sum)

x_female_counts<-merge(female_x_agg, x_param_posts, by.x="symbol", by.y="Gene", all.y=T)
x_female_counts<-x_female_counts[!(x_female_counts$symbol %in% par_param_posts$Gene),]
x_female_counts[(is.na(x_female_counts$lof_obs)),"lof_obs"]<-0
x_female_counts[(is.na(x_female_counts$synonymous_variant)),"synonymous_variant"]<-0

x_hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
x_hs_labels<-paste(head(x_hs_bins,-1),tail(x_hs_bins,-1),sep="-")
x_female_counts$x_hs_bin<-cut(x_female_counts$hs_log10_map, x_hs_bins, labels=x_hs_labels)
x_female_counts<-merge(x_female_counts, denovo_west, by.x="symbol", by.y="symbol",all.x=T)
x_female_counts<-merge(x_female_counts, gnomad_x_genes, by.x="symbol", by.y="gene",all.x=T)
x_female_counts$lof_rate<-ifelse(is.na(x_female_counts$p_lof),x_female_counts$mu_lof,x_female_counts$p_lof)
x_female_counts$syn_rate<-ifelse(is.na(x_female_counts$p_syn),x_female_counts$mu_lof,x_female_counts$p_syn)
x_female_counts$lofexpected<-3* x_female_counts$lof_rate * (1/(2+alpha)) * 2 * female_samp_size
x_female_counts$synexpected<-3* x_female_counts$mu_syn * (1/(2+alpha)) * 2 * female_samp_size

x_female_counts[(is.na(x_female_counts$lofexpected)),"lofexpected"]<-0
x_female_counts[(is.na(x_female_counts$synexpected)),"synexpected"]<-0
x_hs_female_agg<-aggregate(list(lof_expect=x_female_counts$lofexpected,lof_obs=x_female_counts$lof_obs,syn_expect=x_female_counts$synexpected,syn_obs=x_female_counts$synonymous_variant.x),by=list(x_hs_bin=x_female_counts$x_hs_bin),sum)
x_chr_agg<-merge(x_hs_male_agg, x_hs_female_agg, by="x_hs_bin")
x_chr_agg$lof_expect<-x_chr_agg$lof_expect.x + x_chr_agg$lof_expect.y
x_chr_agg$lof_obs<-x_chr_agg$lof_obs.x + x_chr_agg$lof_obs.y
x_chr_agg$syn_expect<-x_chr_agg$syn_expect.x + x_chr_agg$syn_expect.y
x_chr_agg$syn_obs<-x_chr_agg$syn_obs.x + x_chr_agg$syn_obs.y

x_chr_agg$lof_enrich<-x_chr_agg$lof_obs/x_chr_agg$lof_expect
x_chr_agg$syn_enrich<-x_chr_agg$syn_obs/x_chr_agg$syn_expect
x_chr_agg$lof_prob<-(x_chr_agg$lof_obs - x_chr_agg$lof_expect)/x_chr_agg$lof_obs

iter_x_hs_enrich<-list()
iter_x_hs_syn_enrich<-list()
iter_x_hs_prob<-list()
for(i in seq(1:1000)){
  male_iter_df<-list()
  female_iter_df<-list()
  count=1
  for(j in unique(x_male_counts$x_hs_bin)){
    binned_hs<-x_male_counts[x_male_counts$x_hs_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    male_iter_df[[count]]<-resampled_hs
    
    binned_hs<-x_female_counts[x_female_counts$x_hs_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    female_iter_df[[count]]<-resampled_hs
    
    count=count+1
  }
  male_total_iter_df<-do.call(rbind, male_iter_df)
  female_total_iter_df<-do.call(rbind, female_iter_df)
  iter_x_hs_male_agg<-aggregate(list(lof_expect=male_total_iter_df$lofexpected,lof_obs=male_total_iter_df$lof_obs,syn_expect=male_total_iter_df$synexpected,syn_obs=male_total_iter_df$synonymous_variant.x),by=list(x_hs_bin=male_total_iter_df$x_hs_bin),sum)
  iter_x_hs_female_agg<-aggregate(list(lof_expect=female_total_iter_df$lofexpected,lof_obs=female_total_iter_df$lof_obs,syn_expect=female_total_iter_df$synexpected,syn_obs=female_total_iter_df$synonymous_variant.x),by=list(x_hs_bin=female_total_iter_df$x_hs_bin),sum)
  
  iter_x_chr_agg<-merge(iter_x_hs_male_agg, iter_x_hs_female_agg, by="x_hs_bin")
  iter_x_chr_agg$lof_expect<-iter_x_chr_agg$lof_expect.x + iter_x_chr_agg$lof_expect.y
  iter_x_chr_agg$lof_obs<-iter_x_chr_agg$lof_obs.x + iter_x_chr_agg$lof_obs.y
  iter_x_chr_agg$syn_expect<-iter_x_chr_agg$syn_expect.x + iter_x_chr_agg$syn_expect.y
  iter_x_chr_agg$syn_obs<-iter_x_chr_agg$syn_obs.x + iter_x_chr_agg$syn_obs.y
  
  iter_x_chr_agg$lof_enrich<-iter_x_chr_agg$lof_obs/iter_x_chr_agg$lof_expect
  iter_x_chr_agg$syn_enrich<-iter_x_chr_agg$syn_obs/iter_x_chr_agg$syn_expect
  iter_x_chr_agg$prob<-(iter_x_chr_agg$lof_obs - iter_x_chr_agg$lof_expect)/iter_x_chr_agg$lof_obs
  
  iter_x_hs_enrich[[i]]<-iter_x_chr_agg$lof_enrich
  iter_x_hs_syn_enrich[[i]]<-iter_x_chr_agg$syn_enrich
  iter_x_hs_prob[[i]]<-iter_x_chr_agg$prob
}

iter_x_hs_enrich<-data.frame(iter_x_hs_enrich)
iter_x_hs_syn_enrich<-data.frame(iter_x_hs_syn_enrich)
iter_x_hs_prob<-data.frame(iter_x_hs_prob)
names(iter_x_hs_enrich)<-seq(1:1000)
names(iter_x_hs_prob)<-seq(1:1000)
boot_x_hs_enrich<-apply(iter_x_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_hs_enrich)<-c("lof_lo", "lof_hi")
boot_x_hs_syn_enrich<-apply(iter_x_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_hs_syn_enrich)<-c("syn_lo", "syn_hi")
x_hs_enrich_df<-cbind(x_chr_agg,t(boot_x_hs_enrich),t(boot_x_hs_syn_enrich))
colnames(x_hs_enrich_df)[1]<-"hsbin"
boot_x_hs_prob<-apply(iter_x_hs_prob, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_x_hs_prob)<-c("prob_lo", "prob_hi")
x_prob_df<-cbind(x_chr_agg[,c("x_hs_bin","lof_prob")], t(boot_x_hs_prob))
x_hs_enrich_df <- x_hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(x_hs_enrich_df)[2]<-"DNM type"
x_hs_plot<-(ggplot(x_hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
            + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
            + ylim(c(0,16)) + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray")
            + geom_text(aes(label=ifelse(`DNM type`=="lof",obs,''), hjust=1.25,vjust=-2.5),size=3)
            + geom_text(aes(label=ifelse(`DNM type`=="syn",obs,''), hjust=1.25,vjust=2.5),size=3))
x_hs_plot