#Batched model for all salmon stocks
library(here);library(dplyr)
stock_dat<- read.csv(here('data','salmon_productivity_compilation_may2023.csv'))
stock_info<- read.csv(here('data','all_stocks_info_may2023.csv'))
library(rstan);library(loo);library(dlm)

source(here('code','functions.R'))

##steps:
#fits saila-lorde (modified Ricker model) to 244 P. salmon stocks
#plot S-R fits from model, overview of 'c' depensation parameter (c > 1 = depensation)

stock_d=distinct(stock_dat,stock.id,.keep_all = T)
stock_name=paste(stock_d$stock,stock_d$species,sep='-')
depensation_summary=data.frame(stock=stock_name,sp=stock_d$species,c_median=NA,c.l90=NA,c.u90=NA,c.l95=NA,c.u95=NA,Sk=NA,Sk.l90=NA,Sk.u90=NA)
for(i in 1:length(stock_name)){
  stock1<- subset(stock_dat,stock.id==unique(stock.id)[i])
  if(any(stock1$spawners==0)||any(stock1$recruits==0)){
    stock1$spawners=stock1$spawners+1
    stock1$recruits=stock1$recruits+1
  }
  
  sl_test=rstan::stan(file = here('code','Perala depensation stan models','sl_uniform_sigma2_model.stan'), 
                      data=list(R=stock1$recruits,
                                S=stock1$spawners,
                                T=nrow(stock1),
                                kLow=1,
                                kUp=max(stock1$recruits)*3,
                                SkLow=1,
                                SkUp=max(stock1$spawners)*3,
                                sigmaLow=0,
                                sigmaUp=3,
                                cLow=0,
                                cUp=5),
                      pars=c('k','Sk','sigma2','c','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)
  
  params=rstan::extract(sl_test)
  
  depensation_summary[i,3]=median(params$c)
  depensation_summary[i,4]=quantile(params$c,0.05)
  depensation_summary[i,5]=quantile(params$c,0.95)
  depensation_summary[i,6]=quantile(params$c,0.025)
  depensation_summary[i,7]=quantile(params$c,0.975)
  depensation_summary[i,8]=median(params$Sk)
  depensation_summary[i,9]=quantile(params$Sk,0.05)
  depensation_summary[i,10]=quantile(params$Sk,0.95)
  
  sr_plot(params=params,x=stock1,pdf=1,i=i)
  
}

write.csv(depensation_summary,here('outputs','tables','depensation_parameter_overview.csv'))


pdf(here('outputs','figures','c_parameter_summary_plot.pdf'),width=14,height=8.5)
plot(c_median~seq(1:nrow(depensation_summary)),data=depensation_summary,ylim=c(min(na.omit(depensation_summary$c.l90)),max(na.omit(depensation_summary$c.u90))),type='n',bty='l',xlab='',xaxt='n',ylab='c parameter estimate')
abline(h=1,lty=3)
palette(c('darkgray','darkgreen','navy','salmon','darkred'))
for(i in 1:nrow(depensation_summary)){
  lines(c(depensation_summary$c.l90[i],depensation_summary$c.u90[i])~c(rep(i,2)))
}
points(c_median~seq(1:nrow(depensation_summary)),bg=as.factor(depensation_summary$sp),data=depensation_summary,cex=1,pch=21)
dev.off()

pdf(here('outputs','figures','c_vs_K_plot.pdf'),width=11,height=8.5)
plot(c_median~log10(Sk),data=depensation_summary,bty='l',ylab='Median c parameter estimate',xlab='Median log10(K)',type='n')
abline(h=1,lty=3)
points(c_median~log10(Sk),bg=as.factor(depensation_summary$sp),data=depensation_summary,cex=1,pch=21)
dev.off()


##extra
stock1<- subset(stock_dat,stock.id==unique(stock.id)[1])

#
bh_test=rstan::stan(file = here('code','Perala depensation stan models','bh_uniform_sigma2_model.stan'), 
                    data=list(R=stock1$recruits,
                              S=stock1$spawners,
                              T=nrow(stock1),
                              q=0.5,
                              RinfLow=1,
                              RinfUp=max(stock1$recruits)*3,
                              SqLow=1,
                              SqUp=max(stock1$spawners)*3,
                              sigmaLow=0,
                              sigmaUp=3),
                    pars=c('Rinf','Sq','sigma2','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sbh_test=rstan::stan(file = here('code','Perala depensation stan models','sbh_uniform_sigma2_model.stan'), 
                     data=list(R=stock1$recruits,
                               S=stock1$spawners,
                               T=nrow(stock1),
                               q=0.5,
                               RinfLow=1,
                               RinfUp=max(stock1$recruits)*3,
                               SqLow=1,
                               SqUp=max(stock1$spawners)*3,
                               sigmaLow=0,
                               sigmaUp=3,
                               cLow=0,
                               cUp=5),
                     pars=c('Rinf','Sq','sigma2','c','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

# stan code
bh_test2=rstan::stan(file = here('code','stan code','beverton_holt.stan'), 
                     data=list(R=stock1$recruits,
                               S=stock1$spawners,
                               N=nrow(stock1),
                               TT=as.numeric(stock1$broodyear)),
                     pars=c('a','b','sigma_e','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

shep_test2=rstan::stan(file = here('code','stan code','shepherd_sr.stan'), 
                       data=list(R=stock1$recruits,
                                 S=stock1$spawners,
                                 N=nrow(stock1),
                                 TT=as.numeric(stock1$broodyear)),
                       pars=c('a','b','c','sigma_e','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


ricker_test=rstan::stan(file = here('code','Perala depensation stan models','r_uniform_sigma2_model.stan'), 
                        data=list(R=stock1$recruits,
                                  S=stock1$spawners,
                                  T=nrow(stock1),
                                  kLow=1,
                                  kUp=max(stock1$recruits)*3,
                                  SkLow=1,
                                  SkUp=max(stock1$spawners)*3,
                                  sigmaLow=0,
                                  sigmaUp=3),
                        pars=c('k','Sk','sigma2','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sl_test=rstan::stan(file = here('code','Perala depensation stan models','sl_uniform_sigma2_model.stan'), 
                    data=list(R=stock1$recruits,
                              S=stock1$spawners,
                              T=nrow(stock1),
                              kLow=1,
                              kUp=max(stock1$recruits)*3,
                              SkLow=1,
                              SkUp=max(stock1$spawners)*3,
                              sigmaLow=0,
                              sigmaUp=3,
                              cLow=0,
                              cUp=5),
                    pars=c('k','Sk','sigma2','c','mu'), control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


###Testing HMM framework
ale_hmm=rstan::stan_model(file=here('code','HMM','allee_hmm.stan'))


df<- subset(stock_dat,stock.id==unique(stock.id)[1])
df$logRS=log(df$recruits/df$spawners)
df$S=df$spawners
df$R=df$recruits

ahm_test=rstan::stan(here('code','HMM','allee_hmm.stan'), 
                    data=list(N=nrow(df),
                              R_S=df$logRS,
                              S=df$spawners,
                              K=2,
                              pSmax_mean=max(df$spawners)*0.5,
                              pSmax_sig=max(df$spawners)*0.5
                    ),
                    control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

x_new=seq(min(df$spawners),max(df$spawners),length.out=200)
pred_df=data.frame(x_new)

post=rstan::extract(ahm_test)
pred_df[,2]=exp(median(post$log_a[,1])-median(post$b)*x_new)*x_new
pred_df[,3]=exp(median(post$log_a[,2])-median(post$b)*x_new)*x_new
df$gamma=apply(post$gamma[,,2],2,median)
gamma_df=data.frame(by=df$broodyear,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))

plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
  ylim(0,1)+
  geom_hline(yintercept=0.5,linetype='dashed')+
  geom_line(aes(x=by,y=gamma),linewidth=1.3)+
  geom_point(aes(colour = gamma),size=4) +
  scale_colour_viridis_c(name='p')+
  geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
  xlab("Year") + 
  ylab("Prob. of high prod. regime")+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

plot1=ggplot2::ggplot(df, aes(S, R)) +
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
  geom_point(aes(colour = gamma_df$gamma),size=2.5) +
  scale_colour_viridis_c(name='p')+
  xlab("Spawners") + 
  ylab("Recruits")+
  theme_classic(14)+
  xlim(0, max(df$S))+
  ylim(0, max(df$R))+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

legend = cowplot::get_legend(plot1)

plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                              plot2 + theme(legend.position="none"),
                              ncol=2,nrow=1,labels=c("A","B"))
plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
