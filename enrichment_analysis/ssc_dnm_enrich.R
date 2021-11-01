#File paths will need to be changed to run locally!

#Packages
library(tidyr)
library(ggplot2)
library("HistogramTools")
library("ggbeeswarm")
require(gridExtra)
library(dplyr)
library("tidyverse")

#Raw SSC total DNMs
raw_ssc_dnms<-read.csv("raw_ssc_dnms.csv",skip=1)
head(raw_ssc_dnms)
table(raw_ssc_dnms$Consequence)
table(raw_ssc_dnms$Chr)
tmp_ssc_dnms<-subset(raw_ssc_dnms, (raw_ssc_dnms$GENCODE_LoF==1)|(raw_ssc_dnms$Consequence=="synonymous_variant"))
tmp_ssc_dnms<-subset(tmp_ssc_dnms, !(tmp_ssc_dnms$Consequence=="frameshift_variant"))

ssc_samples<-data.frame(do.call(rbind, strsplit(as.character(tmp_ssc_dnms$SampleID), '\\.')))
names(ssc_samples)<-c("ID","patient")
ssc_dnms<-cbind(tmp_ssc_dnms, ssc_samples)
ssc_dnms<-raw_ssc_dnms[,c("SYMBOL","Chr","Pos","Fam","SampleID","Pheno","Consequence","GENCODE_LoF")]
ssc_samples<-data.frame(do.call(rbind, strsplit(as.character(ssc_dnms$SampleID), '\\.')))
names(ssc_samples)<-c("ID","patient")
ssc_dnms<-cbind(ssc_dnms, ssc_samples)

paired_ssc_dnms <- ssc_dnms %>% group_by(Fam) %>% mutate(unique_patients = n_distinct(patient)) %>% filter(unique_patients>1) %>%ungroup()  %>% group_by(Pheno,ID) %>% arrange(desc(Pheno), ID, patient) %>% filter(case_when(Pheno=="control" ~ patient == first(patient), T ~ unique_patients>1))
head(paired_ssc_dnms)
lof_ssc_dnms<-subset(paired_ssc_dnms, paired_ssc_dnms$Consequence %in% c("synonymous_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant", "stop_lost", "start_lost", "stop_gained&splice_region_variant"))
lof_ssc_dnms<-lof_ssc_dnms[!(duplicated(lof_ssc_dnms[c("SYMBOL","SampleID")])),]
lof_ssc_dnms$obs_lof<-ifelse(lof_ssc_dnms$GENCODE_LoF=="1", 1, 0)
lof_ssc_dnms$obs_lof[is.na(lof_ssc_dnms$obs_lof)]<-0
lof_ssc_dnms$obs_syn<-ifelse(lof_ssc_dnms$Consequence=="synonymous_variant",1, 0)

ssc_cases<-subset(lof_ssc_dnms, lof_ssc_dnms$Pheno=="case")
ssc_conts<-subset(lof_ssc_dnms, lof_ssc_dnms$Pheno=="control")
length(unique(paired_ssc_dnms[paired_ssc_dnms$Pheno=="case",]$SampleID))
length(unique(paired_ssc_dnms[paired_ssc_dnms$Pheno=="control",]$SampleID))
ssc_case_counts<-aggregate(x=ssc_cases[,c("obs_lof","obs_syn")], by=list(ssc_cases$SYMBOL), FUN=sum)
ssc_control_counts<-aggregate(x=ssc_conts[,c("obs_lof","obs_syn")], by=list(ssc_conts$SYMBOL), FUN=sum)
names(ssc_case_counts)<-c("gene","obs_case_lof","obs_case_syn")
ssc_case_counts[,c("obs_case_lof","obs_case_syn")]<-(ssc_case_counts[,c("obs_case_lof","obs_case_syn")])
names(ssc_control_counts)<-c("gene","obs_cont_lof","obs_cont_syn")
ssc_control_counts[,c("obs_cont_lof","obs_cont_syn","obs_cont_mis")]<-(ssc_control_counts[,c("obs_cont_lof","obs_cont_syn")])
ssc_counts<-merge(ssc_case_counts, ssc_control_counts, all=T, by="gene")
ssc_counts[is.na(ssc_counts)]<-0
ssc_counts<-merge(ssc_counts, updated_sim_input, by.x="gene",by.y="SYMBOL")
ssc_counts<-merge(ssc_counts, canonical_gnomad, by="gene")
ssc_counts<-merge(ssc_counts, exp_posteriors, by.x="gene", by.y="Gene")
ssc_counts$exp_syn<-(ssc_counts$mu_syn*2*3689)
ssc_counts$exp_lof<-(ssc_counts$mu_lof.x*2*3689)

hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))

hs_labels<-paste(head(hs_bins,-1),tail(hs_bins,-1),sep="-")
ssc_counts$hs_bin<-cut(ssc_counts$log10_map, hs_bins, labels=hs_labels)

hs_ssc_mode_agg<-aggregate(list(lof_case_obs=ssc_counts$obs_case_lof,lof_cont_obs=ssc_counts$obs_cont_lof,syn_case_obs=ssc_counts$obs_case_syn,syn_cont_obs=ssc_counts$obs_cont_syn,possible_lofs=ssc_counts$possible_lof),by=list(hs_bin=ssc_counts$hs_bin),sum)
ssc_gene_counts<-table(ssc_counts$hs_bin)

hs_ssc_mode_agg$lof_enrich<-(hs_ssc_mode_agg$lof_case_obs/sum(hs_ssc_mode_agg$lof_case_obs))/(hs_ssc_mode_agg$lof_cont_obs/sum(hs_ssc_mode_agg$lof_cont_obs))
hs_ssc_mode_agg$syn_enrich<-(hs_ssc_mode_agg$syn_case_obs/sum(hs_ssc_mode_agg$syn_case_obs))/(hs_ssc_mode_agg$syn_cont_obs/sum(hs_ssc_mode_agg$syn_cont_obs))
hs_ssc_mode_agg$lof_enrich<-(hs_ssc_mode_agg$lof_case_obs)/(hs_ssc_mode_agg$lof_cont_obs)
hs_ssc_mode_agg$lof_prob<-(hs_ssc_mode_agg$lof_case_obs-hs_ssc_mode_agg$lof_cont_obs)/(hs_ssc_mode_agg$lof_case_obs)
hs_ssc_mode_agg$syn_enrich<-(hs_ssc_mode_agg$syn_case_obs)/(hs_ssc_mode_agg$syn_cont_obs)
hs_ssc_mode_agg$syn_prob<-(hs_ssc_mode_agg$syn_case_obs-hs_ssc_mode_agg$syn_cont_obs)/(hs_ssc_mode_agg$syn_case_obs)

hs_ssc_mode_agg$prob_causal<-(hs_ssc_mode_agg$lof_case_obs - hs_ssc_mode_agg$lof_cont_obs)/(hs_ssc_mode_agg$lof_case_obs) * (hs_ssc_mode_agg$syn_case_obs/hs_ssc_mode_agg$syn_cont_obs)

ssc_hs_enrich<-list()
ssc_hs_syn_enrich<-list()
ssc_hs_prob<-list()
#Change the number range below for number of bootstraps (100 is for quick, 1000 takes longer but for final figures)
for(i in seq(1:1000)){
  iter_df<-list()
  count=1
  for(j in unique(ssc_counts$hs_bin)){
    binned_hs<-ssc_counts[ssc_counts$hs_bin==j,]
    resampled_hs<-binned_hs[sample(nrow(binned_hs),nrow(binned_hs),replace=T),]
    iter_df[[count]]<-resampled_hs
    count=count+1
  }
  total_iter_df<-do.call(rbind, iter_df)
  iter_mode_agg<-aggregate(list(lof_case_obs=total_iter_df$obs_case_lof,lof_cont_obs=total_iter_df$obs_cont_lof,syn_case_obs=total_iter_df$obs_case_syn,syn_cont_obs=total_iter_df$obs_cont_syn),by=list(hs_bin=total_iter_df$hs_bin),sum)
  iter_mode_agg$enrich<-((iter_mode_agg$lof_case_obs))/((iter_mode_agg$lof_cont_obs+1))
  iter_mode_agg$syn_enrich<-((iter_mode_agg$syn_case_obs))/((iter_mode_agg$syn_cont_obs+1))
  iter_mode_agg$prob<-(iter_mode_agg$lof_case_obs + 1 - iter_mode_agg$lof_cont_obs)/(iter_mode_agg$lof_case_obs + 1) * (iter_mode_agg$syn_case_obs/iter_mode_agg$syn_cont_obs)
  ssc_hs_enrich[[i]]<-iter_mode_agg$enrich
  ssc_hs_syn_enrich[[i]]<-iter_mode_agg$syn_enrich
  ssc_hs_prob[[i]]<-iter_mode_agg$prob
}
ssc_hs_enrich<-data.frame(ssc_hs_enrich)
ssc_hs_syn_enrich<-data.frame(ssc_hs_syn_enrich)
ssc_hs_prob<-data.frame(ssc_hs_prob)
names(ssc_hs_enrich)<-seq(1:1000)
names(ssc_hs_prob)<-seq(1:1000)
boot_hs_enrich<-apply(ssc_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_enrich)<-c("lof_lo", "lof_hi")
boot_hs_syn_enrich<-apply(ssc_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_syn_enrich)<-c("syn_lo", "syn_hi")
boot_ssc_hs_prob<-apply(ssc_hs_prob, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_ssc_hs_prob)<-c("prob_lo", "prob_hi")
ssc_prob_df<-cbind(hs_ssc_mode_agg[,c("hs_bin","lof_prob")], t(boot_ssc_hs_prob))
ssc_hs_enrich_df<-cbind(hs_ssc_mode_agg[,c("hs_bin","lof_enrich","syn_enrich")],t(boot_hs_enrich),t(boot_hs_syn_enrich))
colnames(ssc_hs_enrich_df)[1]<-"hsbin"
ssc_hs_enrich_df <- ssc_hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(ssc_hs_enrich_df)[2]<-"DNM type"
ssc_hs_enrich_df<-ssc_hs_enrich_df[!(ssc_hs_enrich_df$`DNM type`=='mis'),]

ssc_hs_plot<-(ggplot(ssc_hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
              + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
              + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray") + ylim(c(0,15)))

#Show plot for hs (point estimate) enrichment
ssc_hs_plot

#Show plot for prob of interest
ssc_prob_df[!(is.finite(ssc_prob_df$prob_lo)),]$prob_lo<- -1
ssc_prob_df$label="prob"
ssc_prob_plot<-ggplot(ssc_prob_df, aes(x=hs_bin, y=lof_prob, group=label, colour=label)) + geom_point() + geom_line() + geom_errorbar(aes(ymin=(prob_lo),ymax=(prob_hi)),width=0.1)+ theme_bw() + ylim(c(-0.6,1)) + geom_text(aes(label=ssc_gene_counts, hjust=1.25,vjust=-2.5),size=3)
ssc_prob_plot