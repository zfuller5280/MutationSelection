#Packages
library(tidyr)
library(ggplot2)
library("HistogramTools")
library("ggbeeswarm")
require(gridExtra)
library(dplyr)
library(MASS)

#X Chr Posteriors
avg_x_param_posts<-read.table("/VOLUMES/Seagate/abc_hs/expanded_outfiles/x_chr_avg_hs.txt")
names(avg_x_param_posts)<-c("Gene", "s_unscaled_hpd_low", "s_unscaled_hpd_high", "s_unscaled_hpd_map", "s_log10_ci_low", "s_log10_ci_high", "s_log10_map","hs_unscaled_hpd_low", "hs_unscaled_hpd_high", "hs_unscaled_hpd_map", "hs_log10_ci_low", "hs_log10_ci_high", "hs_log10_map","avg_hs_unscaled_hpd_low", "avg_hs_unscaled_hpd_high", "avg_hs_unscaled_hpd_map", "avg_hs_log10_ci_low", "avg_hs_log10_ci_high", "avg_hs_log10_map")
avg_x_param_posts<-subset(avg_x_param_posts, avg_x_param_posts$Gene %in% x_param_posts$Gene)
avg_x_param_posts<-avg_x_param_posts[!(duplicated(avg_x_param_posts$Gene)),]

#Autosome Posteriors
exp_posteriors<-read.table("/Volumes/SeagateBU/Backups.backupdb/Zach’s iMac/2022-05-25-074104/Macintosh HD - Data/Users/zachfuller/autosome_hs.summary_stats.3_30.tsv",header=T)

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
abline(v=-2,lty=3,col="red")
rbPal <- colorRampPalette(c('blue','cyan'))
exp_posteriors$ci_width<-exp_posteriors$log10_ci_high - exp_posteriors$log10_ci_low
ci_levels<-seq(6,0,-0.1)
exp_posteriors$Col <- rbPal(length(ci_levels))[as.numeric(cut(exp_posteriors$ci_width,breaks=ci_levels))]
segments(exp_posteriors$log10_ci_low,exp_posteriors$ecdf,exp_posteriors$log10_ci_high,lwd=.125,col=exp_posteriors$Col)
legend.col(col = rbPal(length(ci_levels)), lev = ci_levels)

#Distribution of hs for X chromosome
x_P<-ecdf(avg_x_param_posts$avg_hs_log10_map)
avg_x_param_posts$ecdf<-x_P(avg_x_param_posts$avg_hs_log10_map)
min(avg_x_param_posts$ecdf)
avg_x_param_posts$ci_width<-avg_x_param_posts$avg_hs_log10_ci_high - avg_x_param_posts$avg_hs_log10_ci_low
rbPal <- colorRampPalette(c('blue','cyan'))
avg_x_param_posts$Col <- rbPal(length(ci_levels))[as.numeric(cut(avg_x_param_posts$ci_width,breaks=ci_levels))]
#Make plot for X
plot(avg_x_param_posts$avg_hs_log10_map,avg_x_param_posts$ecdf,xlim=c(-7,0),cex=0.5,xlab="log10(hs)",ylab="Cumulative Probability",main="X Chromosome",ylim=c(0,1))
abline(v=-2,lty=3,col="red")
segments(avg_x_param_posts$avg_hs_log10_ci_low,avg_x_param_posts$ecdf,avg_x_param_posts$avg_hs_log10_ci_high,lwd=.35,col=avg_x_param_posts$Col)
legend.col(col = rbPal(length(ci_levels)), lev = ci_levels)

#Distribution of hs for X chromosome
x_P<-ecdf(avg_x_param_posts$hs_log10_map)
avg_x_param_posts$ecdf<-x_P(avg_x_param_posts$hs_log10_map)
min(avg_x_param_posts$ecdf)
avg_x_param_posts$ci_width<-avg_x_param_posts$hs_log10_ci_high - avg_x_param_posts$hs_log10_ci_low
rbPal <- colorRampPalette(c('blue','cyan'))
avg_x_param_posts$Col <- rbPal(length(ci_levels))[as.numeric(cut(avg_x_param_posts$ci_width,breaks=ci_levels))]
#Make plot for X
plot(avg_x_param_posts$hs_log10_map,avg_x_param_posts$ecdf,xlim=c(-7,0),cex=0.5,xlab="log10(hs)",ylab="Cumulative Probability",main="X Chromosome",ylim=c(0,1))
abline(v=-2,lty=3,col="red")
segments(avg_x_param_posts$hs_log10_ci_low,avg_x_param_posts$ecdf,avg_x_param_posts$hs_log10_ci_high,lwd=.35,col=avg_x_param_posts$Col)
legend.col(col = rbPal(length(ci_levels)), lev = ci_levels)
