library(readr)
library(tidyverse)
library(ggpubr)
library(ggtext)
####Import dataset: environmental####
setwd("C:/Users/manna/OneDrive - OGS/dpp/new R")
metadata_DPP <- read_csv("C:/Users/manna/OneDrive - OGS/dpp/New Folder/metadata DPP c1_1592022.csv")
View(metadata_DPP)
#Define factor and levels for plotting and handling
colnames(metadata_DPP)
metadata_DPP<-metadata_DPP%>%
  mutate(depth=factor(depth,levels=c("S","B"),
                      ordered = TRUE),
         month=factor(month,levels=c("Jan","Feb","Mar","Apr","May","Jun",
                                     "Jul","Aug","Sep","Oct","Nov","Dec"),
                      ordered=TRUE),
         Season=factor(Season,levels=c("Spr","Sum","Aut","Win"),
                       ordered = TRUE),
         year=factor(year,levels=c("2018","2019","2020","2021"),
                     ordered=TRUE),
         mm_yy=factor(
           interaction(month,year,sep="_"),
           levels=c("Oct_2018","Nov_2018","Dec_2018",
                    "Jan_2019","Feb_2019","Mar_2019","Apr_2019","May_2019","Jun_2019",
                    "Jul_2019","Aug_2019","Sep_2019","Oct_2019","Nov_2019","Dec_2019",
                    "Jan_2020","Feb_2020","Mar_2020","Apr_2020","May_2020","Jun_2020",
                    "Jul_2020","Aug_2020","Sep_2020","Oct_2020","Nov_2020","Dec_2020",
                    "Jan_2021","Feb_2021","Mar_2021","Apr_2021")))

#unify unit of measurement. DPP and HCP both ugC L d
metadata_DPP$dpp_part_ug<-metadata_DPP$dpp_part/1000
metadata_DPP$dpp_part_ug_sd<-metadata_DPP$dpp_part_sd/1000

metadata_DPP$dpp_diss_ug<-metadata_DPP$dpp_diss/1000
metadata_DPP$dpp_diss_ug_sd<-metadata_DPP$dpp_diss_sd/1000

metadata_DPP$dpp_tot_ug<-metadata_DPP$dpp_tot/1000
metadata_DPP$dpp_tot_ug_sd<-metadata_DPP$dpp_tot_sd/1000

metadata_DPP$hcp_ugd<-metadata_DPP$hcp*24
metadata_DPP$hcp_ugd_sd<-metadata_DPP$std_hcp*24

metadata_DPP<-metadata_DPP%>%mutate(cn_diss=doc/don,
                                    cn_part=poc/pn,
                                    dpp_hcp_ratio=dpp_part_ug/hcp_ugd,
                                    amoA=ifelse(amoA_status=="POS",amoA_tot,NA))
metadata_DPP$amoA_SD<-as.numeric(metadata_DPP$amoA_SD)
####Drawing Figure 1####
#Composite image of dDIC fixation and HCP rates (all data+seasonal)
dpp_seasonal<-metadata_DPP%>%
  ggplot(aes(x=month,y=dpp_part_ug))+
  geom_point(aes(color=depth),size=3)+
  scale_color_manual(name = "Depth",
                     values=c("#FCCf14","#008082"),
                     labels=c("Surface", "Bottom"))+
  scale_fill_manual(name = "Depth",
                    values=c("#FCCf14","#008082"),
                    labels=c("Surface", "Bottom"))+
  geom_errorbar(aes(ymin=(dpp_part_ug-dpp_part_ug_sd),
                    ymax=(dpp_part_ug+dpp_part_ug_sd)),
                color="black",
                size=0.8,
                width=0.5)+
  geom_smooth(method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              lwd=1.7,
              aes(x=month,
                  group=depth,
                  color=depth,
                  fill=depth))+
  facet_wrap(~depth,nrow = 2,
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  #ylab("dDIC fixation (μg C L<sup>-1</sup> d<sup>-1</sup>)")+
  #xlab("\nMonth")+
  theme(axis.title = element_markdown(size=13,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1))
dpp_seasonal

dpp_all<-metadata_DPP%>%
  ggplot(aes(x=mm_yy,y=dpp_part_ug,group=depth))+
  geom_point(aes(color=depth),size=4)+
  geom_line(aes(color=depth),size=0.8)+
  ylab("dDIC fixation (μg C L<sup>-1</sup> d<sup>-1</sup>)\n")+
  xlab("Month_Year")+
  scale_color_manual(name = "Depth",
                     values=c("#FCCf14","#008082"),
                     labels=c("Surface", "Bottom"))+
  geom_errorbar(aes(ymin=(dpp_part_ug-dpp_part_ug_sd),
                    ymax=(dpp_part_ug+dpp_part_ug_sd)),
                    color="black",
                    size=0.8,
                    width=0.5)+
  facet_wrap(~depth,nrow=2,
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  theme(axis.title = element_markdown(size=13,face="bold"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10),
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.position = "none")

dpp_all

hcp_seasonal<-metadata_DPP%>%
  ggplot(aes(x=month,y=hcp_ugd))+
  geom_point(aes(color=depth),size=3)+
  scale_color_manual(name = "Depth",
                     values=c("#FCCf14","#008082"),
                     labels=c("Surface", "Bottom"))+
  scale_fill_manual(name = "Depth",
                    values=c("#FCCf14","#008082"),
                    labels=c("Surface", "Bottom"))+
  geom_errorbar(aes(ymin=(hcp_ugd-hcp_ugd_sd),
                    ymax=(hcp_ugd+hcp_ugd_sd)),
                color="black",
                size=0.8,
                width=0.5)+
  geom_smooth(method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              lwd=1.7,
              aes(x=month,
                  group=depth,
                  color=depth,
                  fill=depth))+
  facet_wrap(~depth,nrow = 2,
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  #ylab("dDIC fixation (μg C L<sup>-1</sup> d<sup>-1</sup>)")+
  xlab("\nMonth")+
  theme(axis.title = element_markdown(size=13,face="bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1))
hcp_seasonal


hcp_all<-metadata_DPP%>%
  ggplot(aes(x=mm_yy,y=hcp_ugd,group=depth))+
  geom_point(aes(color=depth),size=4)+
  geom_line(aes(color=depth),size=0.8)+
  ylab("HCP (μg C L<sup>-1</sup> d<sup>-1</sup>)\n")+
  xlab("Month_Year")+
  scale_color_manual(name = "Depth",
                     values=c("#FCCf14","#008082"),
                     labels=c("Surface", "Bottom"))+
  geom_errorbar(aes(ymin=(hcp_ugd-hcp_ugd_sd),
                    ymax=(hcp_ugd+hcp_ugd_sd)),
                color="black",
                size=0.8,
                width=0.5)+
  facet_wrap(~depth,nrow=2,
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  theme(axis.title = element_markdown(size=13,face="bold"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.position = "none")
hcp_all


both_graphs<-ggarrange(dpp_all,dpp_seasonal,
                       hcp_all,hcp_seasonal,
                       align = "hv",
                       ncol = 2,nrow = 2,
                       common.legend = TRUE,
                       legend = "bottom",
                       widths = c(1,0.4),
                       labels = c("a","b","c","d"))
both_graphs

####Drawing Figure 3####
#amoA tilemap
amo_tile<-metadata_DPP%>%
  ggplot(aes(x=(month),
             y=(year),
             fill=amoA))+
  geom_tile(aes(color=""))+
  scale_fill_viridis_c(na.value="grey90")+
  scale_color_manual(values = NA)+
  labs(fill="amoA copies L<sup>-1</sup>")+
  guides(colour=guide_legend("NQ", override.aes=list(fill="grey90")))+
  facet_wrap(~depth,
             nrow = 2,
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  xlab("Month")+
  ylab("Year")+
  theme(axis.title = element_markdown(size=13,face="bold"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.title = element_markdown(face="bold"))
amo_tile 

####Drawing Figure S3####
#dDICfix:HCP ratio
ggplot(metadata_DPP,aes(x=mm_yy,y=dpp_hcp_ratio,group=depth))+
  geom_line(aes(color=depth),lwd=3)+
  scale_color_manual(name = "Depth",values=c("#FCCf14","#008082"),labels=c("Surface", "Bottom"))+
  xlab("Month_Year")+
  ylab("dDICfixation:HCP")+
  theme(axis.text.x = element_text(angle=45,hjust =1),
        axis.title =element_text(hjust=0.5,size=14),
        axis.ticks = element_line(linewidth = 1,lineend = "round",),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        axis.line = element_line(colour ="#607d8b"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        legend.position = "bottom")


####Partial least squares regression (PLS)####
library(pls)
library(caret)
#Partial least squares regression
colnames(metadata_DPP)
pls_both<-metadata_DPP%>%
  filter(dpp_part_ug<1)%>%
  column_to_rownames("sample_id")%>%
  dplyr::select(c(46,5,8:9,12,14:16,17,52))

plsFit_both<-plsr(dpp_part_ug~ ., data=pls_both, 
                  validation="LOO",
                  scale=TRUE,
                  model=TRUE)
validationplot(plsFit_both, val.type="RMSEP")#2 component 
#rerun the models with the correct number of components
plsFit_both<-plsr(dpp_part_ug~ ., data=pls_both, 
                  validation="LOO",
                  scale=TRUE,
                  model=TRUE,
                  ncomp=2)
#retrieving fit goodness(R^2)
#both
predict(plsFit_both,ncomp=2,newdata = pls_both)
RMSEP(plsFit_both, newdata = pls_both)
pls.pred2_both = predict(plsFit_both, pls_both, ncomp=2)
plot(pls_both$dpp_part, pls.pred2_both,main="Test Dataset")
abline(0, 1, col="red")
pls.eval_both<-data.frame(obs=pls_both$dpp_part, pred=pls.pred2_both)
colnames(pls.eval_both)<-c("obs","pred")
defaultSummary(pls.eval_both,model = plsFit_both)
#RMSE=0.11 ;R2=0.72;MAE=0.08
#retrieve coefficients
coefficients_both<-(coef(plsFit_both,ncomp=2))
sum.coef_both<-sum(sapply(coefficients_both, abs))
coefficients_both<-coefficients_both * 100 / sum.coef_both
coefficients_both<-sort(coefficients_both[, 1 , 1])
coefficients_both<-data.frame(coefficients_both)
coefficients_both$variables<-factor(rownames(coefficients_both),
                                    levels=c("temp", "nh4", "no2","don",
                                             "doc","poc","pn","chla",
                                             "hcp_ugd"),ordered = TRUE)
####Drawing Figure 2####
#Normalized regression coefficients for pls
coefficients_both%>%
  ggplot(aes(x=variables,y=coefficients_both))+
  geom_bar(stat="identity",
           position="dodge",
           fill="#607d8b")+
  #facet_grid(~category,scales = "free_x",space = "free_x")+
  ylab("Normalized regression coefficients")+
  scale_x_discrete(labels=c("temp" = "Temperature",
                            "nh4"="NH<sup>+</sup><sub>4</sub>",
                            "no2"="NO<sup>-</sup><sub>2</sub>",
                            "don"="DON",
                            "doc"="DOC",
                            "poc"="POC",
                            "pn"="PN",
                            "chla"="Chl *a*",
                            "hcp_ugd"="HCP"))+
  theme(axis.text.x = element_markdown(angle=90,
                                       vjust=0.5,
                                       hjust=1,
                                       face="bold"),
        axis.title.x =element_blank(),
        axis.title = element_text(size=13,face="bold"),
        #axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,size=0.5,colour = "#607d8b"),
        #panel.grid = element_line(size=0.1,colour = "#95a5a6"),
        legend.text = element_text(size=12),
        legend.title =element_text(size=12),
        legend.position = "right")

