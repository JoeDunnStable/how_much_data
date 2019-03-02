require("ggplot2")

setwd("~/Documents/XCode/how_much_data/output")
central_limit_data<-data.frame(alpha=c(1,2,4),kappa_mad=c(1,0,0), n=rep(NA,3))
pt_data<-read.table("pinelis_taleb_dump.out",header=T,strip.white=T,sep="," )
ns <- unique(pt_data$n)
n_levels<-as.character(ns)
pt_data$n = factor(pt_data$n, levels=n_levels)

pdf("pinelis_taleb_graphs.pdf")


qplot(x=param_value, y=kappa_mad, data=pt_data[pt_data$param_name=="alpha",], 
      color=n, facets=distribution ~ .,
      main=expression(paste(kappa, " vs ", alpha, " and n")), 
      xlab=expression(alpha),
      ylab=expression(kappa),
      geom="point", ylim=c(0,1)) +
     geom_line(aes(x=alpha, y=kappa_mad, color=n), data=central_limit_data,
               color="black")

qplot(x=param_value, y=kappa_mad, data=pt_data[pt_data$param_name=="sigma",], 
      color=n, facets=distribution ~ .,
      main=expression(paste(kappa, " vs ", sigma, " and n")), 
      xlab=expression(sigma),
      ylab=expression(kappa),
      geom="point", ylim=c(0,1)) +
  geom_hline(yintercept=0, color="black")

qplot(x=param_value, y=kappa_mad, data=pt_data[pt_data$param_name=="lambda",], 
      color=n, facets=distribution ~ .,
      main=expression(paste(kappa, " vs ", lambda, " and n")), 
      xlab=expression(lambda),
      ylab=expression(kappa),
      geom="point", ylim=c(0,1)) +
  geom_hline(yintercept=0, color="black")

dev.off()