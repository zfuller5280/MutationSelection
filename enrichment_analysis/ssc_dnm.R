#File paths will need to be changed to run locally!

#Packages
library(tidyr)
library(ggplot2)
library("HistogramTools")
library("ggbeeswarm")
require(gridExtra)
library(dplyr)
library("tidyverse")

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

#DDD data from Samocha
denovo_west<-read.table("forZach_synexp_extended_denovoWEST_results.txt",sep="\t",header=T)
head(denovo_west)

head(exp_posteriors)
#Read in DNM data
dnms<-read.table("dnm.lof.validation.txt",header=T)
table(dnms$data)
ddd_dnms<-subset(dnms, dnms$data=="ddd")
ssc_dnms<-subset(dnms, dnms$data=="ssc")

#DDD data from Samocha
denovo_west<-read.table("forZach_synexp_extended_denovoWEST_results.txt",sep="\t",header=T)
head(denovo_west)

#Our hs point estimates
table(ddd_dnms$HC_LOF_canonical)
length(unique(ssc_dnms$SampleID))

lof_ddd_dnms<-subset(ddd_dnms,(ddd_dnms$HC_LOF_canonical=="1")|(ddd_dnms$canonical_csq=="synonymous_variant"))
lof_ddd_dnms<-lof_ddd_dnms[!(duplicated(lof_ddd_dnms[c("SYMBOL","SampleID")])),]
lof_ddd_dnms$obs_lof<-ifelse(lof_ddd_dnms$HC_LOF_canonical=="1", 1, 0)
lof_ddd_dnms$obs_lof[is.na(lof_ddd_dnms$obs_lof)]<-0
lof_ddd_dnms$obs_syn<-ifelse(lof_ddd_dnms$canonical_csq=="synonymous_variant",1, 0)
ddd_counts<-aggregate(x=lof_ddd_dnms[,c("obs_lof","obs_syn")], by=list(lof_ddd_dnms$SYMBOL), FUN=sum)
names(ddd_counts)<-c("gene","obs_lof","obs_syn")
ddd_counts<-merge(ddd_counts, valid_sim_input, by="gene")
ddd_counts<-merge(ddd_counts, canonical_gnomad, by="gene")
ddd_counts<-merge(ddd_counts, exp_posteriors, by.x="gene", by.y="Gene")
ddd_counts$exp_syn<-(ddd_counts$mu_syn.x*2*52664)
ddd_counts$exp_lof<-(ddd_counts$mu_lof.x*2*52664)
hs_bins<-c(-8, log10(c(5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1)))
hs_labels<-paste(head(hs_bins,-1),tail(hs_bins,-1),sep="-")
ddd_counts$hs_bin<-cut(ddd_counts$log10_map, hs_bins, labels=hs_labels)
hs_mode_agg<-aggregate(list(lof_expect=ddd_counts$exp_lof,lof_obs=ddd_counts$obs_lof.x,syn_expect=ddd_counts$exp_syn,syn_obs=ddd_counts$obs_syn.x),by=list(hs_bin=ddd_counts$hs_bin),sum)
hs_mode_agg$lof_enrich<-hs_mode_agg$lof_obs/hs_mode_agg$lof_expect
hs_mode_agg$syn_enrich<-hs_mode_agg$syn_obs/hs_mode_agg$syn_expect
hs_mode_agg
head(ddd_counts)
head(ssc_dnms)
ssc_samples<-data.frame(do.call(rbind, strsplit(as.character(ssc_dnms$SampleID), '\\.')))
names(ssc_samples)<-c("ID","patient")
ssc_dnms<-cbind(ssc_dnms, ssc_samples)
#complete_ssc_dnms <- ssc_dnms %>% group_by(ID) %>% mutate(unique_patients = n_distinct(patient)) %>% filter(unique_patients>1) %>% ungroup() 
paired_ssc_dnms <- ssc_dnms %>% group_by(ID) %>% mutate(unique_patients = n_distinct(patient)) %>% filter(unique_patients>1) %>%ungroup()  %>% group_by(Pheno,ID) %>% arrange(desc(Pheno), ID, patient) %>% filter(case_when(Pheno=="control" ~ patient == first(patient), T ~ unique_patients>1))
head(paired_ssc_dnms)
#lof_ssc_dnms<-subset(paired_ssc_dnms,(paired_ssc_dnms$HC_LOF_canonical=="1")|(paired_ssc_dnms$canonical_csq=="synonymous_variant"))
lof_ssc_dnms<-subset(paired_ssc_dnms,(paired_ssc_dnms$HC_LOF_canonical=="1")|(paired_ssc_dnms$canonical_csq=="synonymous_variant")|(paired_ssc_dnms$canonical_csq=="missense_variant"))
lof_ssc_dnms<-lof_ssc_dnms[!(duplicated(lof_ssc_dnms[c("SYMBOL","SampleID")])),]
lof_ssc_dnms$obs_lof<-ifelse(lof_ssc_dnms$HC_LOF_canonical=="1", 1, 0)
lof_ssc_dnms$obs_lof[is.na(lof_ssc_dnms$obs_lof)]<-0
lof_ssc_dnms$obs_syn<-ifelse(lof_ssc_dnms$canonical_csq=="synonymous_variant",1, 0)
lof_ssc_dnms$obs_mis<-ifelse(lof_ssc_dnms$canonical_csq=="missense_variant",1, 0)
ssc_cases<-subset(lof_ssc_dnms, lof_ssc_dnms$Pheno=="case")
ssc_conts<-subset(lof_ssc_dnms, lof_ssc_dnms$Pheno=="control")
length(unique(paired_ssc_dnms[paired_ssc_dnms$Pheno=="case",]$SampleID))
length(unique(paired_ssc_dnms[paired_ssc_dnms$Pheno=="control",]$SampleID))
ssc_case_counts<-aggregate(x=ssc_cases[,c("obs_lof","obs_syn","obs_mis")], by=list(ssc_cases$SYMBOL), FUN=sum)
ssc_control_counts<-aggregate(x=ssc_conts[,c("obs_lof","obs_syn","obs_mis")], by=list(ssc_conts$SYMBOL), FUN=sum)
names(ssc_case_counts)<-c("gene","obs_case_lof","obs_case_syn","obs_case_mis")
ssc_case_counts[,c("obs_case_lof","obs_case_syn","obs_case_mis")]<-(ssc_case_counts[,c("obs_case_lof","obs_case_syn","obs_case_mis")])
names(ssc_control_counts)<-c("gene","obs_cont_lof","obs_cont_syn","obs_cont_mis")
ssc_control_counts[,c("obs_cont_lof","obs_cont_syn","obs_cont_mis")]<-(ssc_control_counts[,c("obs_cont_lof","obs_cont_syn","obs_cont_mis")])
ssc_counts<-merge(ssc_case_counts, ssc_control_counts, all=T, by="gene")
ssc_counts[is.na(ssc_counts)]<-0
ssc_counts<-merge(ssc_counts, sim_input, by="gene")
ssc_counts<-merge(ssc_counts, canonical_gnomad, by="gene")
ssc_counts<-merge(ssc_counts, exp_posteriors, by.x="gene", by.y="Gene")
ssc_counts$exp_syn<-(ssc_counts$mu_syn.x*2*3689)
ssc_counts$exp_lof<-(ssc_counts$mu_lof.x*2*3689)
#hs_bins<-c(-8, log10(c(1e-6,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1)))
#hs_bins<-c(-8, log10(c(1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1)))
hs_bins<-c(-8, log10(c(1e-4,1e-3,1e-2,1e-1,1)))
#hs_bins<-c(-8, -6,  -3, -2, -1, 0)

hs_labels<-paste(head(hs_bins,-1),tail(hs_bins,-1),sep="-")
ssc_counts$hs_bin<-cut(ssc_counts$log10_map, hs_bins, labels=hs_labels)
hs_ssc_mode_agg<-aggregate(list(lof_case_obs=ssc_counts$obs_case_lof,lof_cont_obs=ssc_counts$obs_cont_lof,syn_case_obs=ssc_counts$obs_case_syn,syn_cont_obs=ssc_counts$obs_cont_syn,mis_case_obs=ssc_counts$obs_case_mis,mis_cont_obs=ssc_counts$obs_cont_mis),by=list(hs_bin=ssc_counts$hs_bin),sum)
hs_ssc_mode_agg$lof_enrich<-(hs_ssc_mode_agg$lof_case_obs/sum(hs_ssc_mode_agg$lof_case_obs))/(hs_ssc_mode_agg$lof_cont_obs/sum(hs_ssc_mode_agg$lof_cont_obs))
hs_ssc_mode_agg$syn_enrich<-(hs_ssc_mode_agg$syn_case_obs/sum(hs_ssc_mode_agg$syn_case_obs))/(hs_ssc_mode_agg$syn_cont_obs/sum(hs_ssc_mode_agg$syn_cont_obs))
hs_ssc_mode_agg$lof_enrich<-(hs_ssc_mode_agg$lof_case_obs)/(hs_ssc_mode_agg$lof_cont_obs)
hs_ssc_mode_agg$syn_enrich<-(hs_ssc_mode_agg$syn_case_obs)/(hs_ssc_mode_agg$syn_cont_obs)
hs_ssc_mode_agg$mis_enrich<-(hs_ssc_mode_agg$mis_case_obs)/(hs_ssc_mode_agg$mis_cont_obs)
#hs_ssc_mode_agg$lof_enrich<-(hs_ssc_mode_agg$lof_case_obs/length((ssc_dnms[ssc_dnms$Pheno=="case",]$SampleID)))/(hs_ssc_mode_agg$lof_cont_obs/length((ssc_dnms[ssc_dnms$Pheno=="control",]$SampleID)))
#hs_ssc_mode_agg$syn_enrich<-(hs_ssc_mode_agg$syn_case_obs/length((ssc_dnms[ssc_dnms$Pheno=="case",]$SampleID)))/(hs_ssc_mode_agg$syn_cont_obs/length((ssc_dnms[ssc_dnms$Pheno=="control",]$SampleID)))
hs_ssc_mode_agg
ssc_hs_enrich<-list()
ssc_hs_syn_enrich<-list()
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
  iter_mode_agg<-aggregate(list(lof_case_obs=total_iter_df$obs_case_lof,lof_cont_obs=total_iter_df$obs_cont_lof,syn_case_obs=total_iter_df$obs_case_syn,syn_cont_obs=total_iter_df$obs_cont_syn,mis_case_obs=ssc_counts$obs_case_mis,mis_cont_obs=ssc_counts$obs_cont_mis),by=list(hs_bin=total_iter_df$hs_bin),sum)
  #iter_mode_agg$enrich<-((iter_mode_agg$lof_case_obs+1)/length((ssc_dnms[ssc_dnms$Pheno=="case",]$SampleID)))/((iter_mode_agg$lof_cont_obs+1)/length((ssc_dnms[ssc_dnms$Pheno=="control",]$SampleID)))
  #iter_mode_agg$syn_enrich<-((iter_mode_agg$syn_case_obs+1)/length((ssc_dnms[ssc_dnms$Pheno=="case",]$SampleID)))/((iter_mode_agg$syn_cont_obs+1)/length((ssc_dnms[ssc_dnms$Pheno=="control",]$SampleID)))
  iter_mode_agg$enrich<-((iter_mode_agg$lof_case_obs))/((iter_mode_agg$lof_cont_obs+1))
  iter_mode_agg$syn_enrich<-((iter_mode_agg$syn_case_obs))/((iter_mode_agg$syn_cont_obs+1))
  ssc_hs_enrich[[i]]<-iter_mode_agg$enrich
  ssc_hs_syn_enrich[[i]]<-iter_mode_agg$syn_enrich
}
ssc_hs_enrich<-data.frame(ssc_hs_enrich)
ssc_hs_syn_enrich<-data.frame(ssc_hs_syn_enrich)
names(ssc_hs_enrich)<-seq(1:1000)
boot_hs_enrich<-apply(ssc_hs_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_enrich)<-c("lof_lo", "lof_hi")
boot_hs_syn_enrich<-apply(ssc_hs_syn_enrich, 1, quantile, probs=c(0.025, 0.975),na.rm=T)
rownames(boot_hs_syn_enrich)<-c("syn_lo", "syn_hi")
ssc_hs_enrich_df<-cbind(hs_ssc_mode_agg,t(boot_hs_enrich),t(boot_hs_syn_enrich))
colnames(ssc_hs_enrich_df)[1]<-"hsbin"
ssc_hs_enrich_df <- ssc_hs_enrich_df %>% pivot_longer(cols = contains("_"), names_to=c("DNM_type",".value"), names_pattern= "(.*?)_(.*)")
colnames(ssc_hs_enrich_df)[2]<-"DNM type"
ssc_hs_enrich_df<-ssc_hs_enrich_df[!(ssc_hs_enrich_df$`DNM type`=='mis'),]
ssc_hs_plot<-(ggplot(ssc_hs_enrich_df, aes(x=hsbin, y=(enrich), colour=`DNM type`)) + geom_point() 
          + geom_line(aes(group=`DNM type`)) + geom_errorbar(aes(ymin=(lo),ymax=(hi)),width=0.1)+ theme_bw() 
          + ylim(c(0,16)) + xlab("hs bin") + ylab("Enrichment")+ geom_hline(yintercept=1, linetype="dashed",color="gray"))
#Show plot for hs (point estimate) enrichment
ssc_hs_plot

grid.arrange(hs_plot, x_hs_plot, ssc_hs_plot, nrow=1)

#Burden of DNMs
table(ssc_cases[ssc_cases$obs_lof>0,]$ID)
ssc_cases_hs<-merge(ssc_cases, exp_posteriors, by.x="SYMBOL", by.y="Gene")
ssc_conts_hs<-merge(ssc_conts, exp_posteriors, by.x="SYMBOL", by.y="Gene")
ssc_case_w <- ssc_cases_hs %>% group_by(ID) %>% filter(obs_lof > 0) %>% mutate(fitness = prod(1-10^log10_map))
ssc_conts_w <- ssc_conts_hs %>% group_by(ID) %>% filter(obs_lof > 0) %>% mutate(fitness = prod(1-10^log10_map))
dev.off()
par(mfrow=c(2,1))
library(HistogramTools)
PlotRelativeFrequency(hist(ssc_case_w$fitness,plot=F),main="Cases",ylab="Frequency",xlab="Fitness",xlim=c(0.5,1),ylim=c(0,0.7))
PlotRelativeFrequency(hist(ssc_conts_w$fitness,plot=F),main="Controls",ylab="Frequency",xlab="Fitness",xlim=c(0.5,1),ylim=c(0,0.7))
boxplot(ssc_case_w$fitness, ssc_conts_w$fitness)
mean(ssc_conts_w$log10_map)
ssc_conts_w[order(ssc_conts_w$fitness),]
1-10^-7
