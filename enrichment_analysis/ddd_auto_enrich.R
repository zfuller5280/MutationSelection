#File paths will need to be changed to run locally!
#Packages
library(tidyr)
library(ggplot2)
library("HistogramTools")
library("ggbeeswarm")
require(gridExtra)
library(dplyr)
library("tidyverse")

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

#Read in our posterior data for hs
exp_posteriors<-read.table("../data_files/autosome_posteriors.8_11.tsv", header=T, sep="\t")
canonical_gnomad<-read.table("../data_files/canonical_gnomad.data.tsv", header=T, sep="\t")

#DDD Autosomal Enrichment Plot
#Our hs point estimates
#DDD data from Samocha
denovo_west<-read.table("./forZach_synexp_extended_denovoWEST_results.txt",sep="\t",header=T)
gnomad_mu<-read.table("../data_files/gnomad.v2.1.1.lof_metrics.by_transcript.txt",header=T,fill=T,sep='\t')
canonical_gnomad<-gnomad_mu[gnomad_mu$canonical=="true",]
denovo_west<-merge(denovo_west, canonical_gnomad, by.x="symbol", by.y="gene")
auto_dnm_counts<-merge(denovo_west, exp_posteriors, by.x="symbol", by.y="Gene")
hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
hs_labels<-paste(head(hs_bins,-1),tail(hs_bins,-1),sep="-")
auto_dnm_counts$hs_bin<-cut(auto_dnm_counts$log10_map, hs_bins, labels=hs_labels)
auto_dnm_counts[(is.na(auto_dnm_counts$lofexpected)),]
auto_dnm_counts[(is.na(auto_dnm_counts$lofexpected)),"lofexpected"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$synexp)),"synexp"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$lofcount)),"lofcount"]<-0
auto_dnm_counts[(is.na(auto_dnm_counts$synonymous_variant)),"synonymous_variant"]<-0
auto_dnm_counts$synexpected<-(auto_dnm_counts$mu_syn *2*31058)
hs_mode_agg<-aggregate(list(lof_expect=auto_dnm_counts$lofexpected,lof_obs=auto_dnm_counts$lofcount,syn_expect=auto_dnm_counts$synexpected,syn_obs=auto_dnm_counts$synonymous_variant),by=list(hs_bin=auto_dnm_counts$hs_bin),sum)
hs_mode_agg$lof_enrich<-hs_mode_agg$lof_obs/hs_mode_agg$lof_expect
hs_mode_agg$syn_enrich<-hs_mode_agg$syn_obs/hs_mode_agg$syn_expect
hs_mode_agg$lof_prob<-(hs_mode_agg$lof_obs - hs_mode_agg$lof_expect)/hs_mode_agg$lof_obs

iter_hs_enrich<-list()
iter_hs_syn_enrich<-list()
iter_hs_prob<-list()
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
  iter_mode_agg$prob<-(iter_mode_agg$obs - iter_mode_agg$expect)/iter_mode_agg$obs
  iter_hs_enrich[[i]]<-iter_mode_agg$enrich
  iter_hs_syn_enrich[[i]]<-iter_mode_agg$syn_enrich
  iter_hs_prob[[i]]<-iter_mode_agg$prob
}
iter_hs_enrich<-data.frame(iter_hs_enrich)
iter_hs_syn_enrich<-data.frame(iter_hs_syn_enrich)
iter_hs_prob<-data.frame(iter_hs_prob)
names(iter_hs_enrich)<-seq(1:1000)
names(iter_hs_prob)<-seq(1:1000)
boot_hs_enrich<-apply(iter_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_enrich)<-c("lof_lo", "lof_hi")
boot_hs_syn_enrich<-apply(iter_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_syn_enrich)<-c("syn_lo", "syn_hi")
hs_enrich_df<-cbind(hs_mode_agg,t(boot_hs_enrich),t(boot_hs_syn_enrich))
colnames(hs_enrich_df)[1]<-"hsbin"
boot_hs_prob<-apply(iter_hs_prob, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_prob)<-c("prob_lo", "prob_hi")
prob_df<-cbind(hs_mode_agg[,c("hs_bin","lof_prob")], t(boot_hs_prob))
hs_enrich_df <- hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(hs_enrich_df)[2]<-"DNM type"
hs_plot<-(ggplot(hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
          + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
          + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray") + ylim(c(0,15)) )

#Show plot for hs (point estimate) enrichment
hs_plot
