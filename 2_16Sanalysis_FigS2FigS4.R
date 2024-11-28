###16S data analysis###
library(phyloseq)
library(microbiome)
C1_all<-readRDS("C:/Users/manna/OneDrive - OGS/dpp/New folder/C1.all.RDS")
#extract feature and tax tables
count_tab<-data.frame(otu_table(C1_all))
tax_tab<-data.frame(tax_table(C1_all))
#store reference sequences for future needs
asv_refseq<-data.frame("seq"=colnames(count_tab),"ASVs"=paste0("ASV",1:ncol(count_tab)))
#16S sequences are used as feature names. This is very difficult to work with, replace with ASV1 ecc.
names(count_tab) <- asv_refseq$ASVs[match(names(count_tab), asv_refseq$seq)]
#Tax table has to be transposed before this step
trasp_tax_tab<-t(tax_tab)
colnames(trasp_tax_tab)<-asv_refseq$ASVs[match(colnames(trasp_tax_tab),asv_refseq$seq)]
tax_tab<-as.data.frame(t(trasp_tax_tab))
#subset feature table only keeping samples from Oct18 to Apr21, when we have measurements of dDIC fix
rownames(count_tab)<-str_replace_all(rownames(count_tab),pattern = "Ago","Aug")
rownames(count_tab)<-str_replace_all(rownames(count_tab),pattern = "Sept","Sep")
count_tab_subs<-count_tab%>%
  filter(rownames(count_tab)%in%metadata_DPP$sample_id)
####Create the phyloseq object####
otu_tab<-otu_table(count_tab_subs,taxa_are_rows = FALSE)
tax_tab<-tax_table(as.matrix(tax_tab))
metadata<-sample_data(metadata_DPP%>%
                        column_to_rownames(var="sample_id"))
dpp_pseq<-phyloseq(otu_tab,tax_tab,metadata)
sample_names(dpp_pseq)
#Revise taxonomy table
#If a cell contains NA, this loop replace the NA with the prefix "unassigned_" followed by the previously known taxonomic level 
#In this way, we have all formatted in the right way and avoid problems in merging by taxonomic levels
tax_tab_torevise<-tax_table(dpp_pseq)

for (col in 1:ncol(tax_tab_torevise)) {
  for (row in 1:nrow(tax_tab_torevise)) {
    if (is.na(tax_tab_torevise[row,col])) {
      if (!grepl("feature_id", tax_tab_torevise[row,col-1]) & !grepl("unassigned", tax_tab_torevise[row,col-1])) {
        tax_tab_torevise[row,col] <- paste0("unassigned_", tax_tab_torevise[row,col-1])
      } else {
        tax_tab_torevise[row,col] <- tax_tab_torevise[row,col-1]
      }
    }
  }
}
#replace taxonomy table
tax_table(dpp_pseq)<-tax_tab_torevise
####Subset to >1% relative abundance in at least one sample
#Transform in relative abundance (0-1 range)
dpp_pseq_rel<-microbiome::transform(dpp_pseq,"compositional")
saveRDS(dpp_pseq_rel,"dpp_pseq_rel.rds")
#Filtering out rare taxa. Keeping only those >1% in at least one sample
taxa_greater1_rel<-filter_taxa(dpp_pseq_rel,function(x) max(x)>0.01,TRUE)
taxa_greater1_rel
#Prepare data for plotting
dpp_pseq_rel1_melt<-psmelt(taxa_greater1_rel)
dpp_pseq_rel1_melt$Abundance<-dpp_pseq_rel1_melt$Abundance*100

####Drawing Figure S2####
#Nitrates+amoA+nitrosopumilus
nitrites<-metadata_DPP%>%
  ggplot(aes(x=mm_yy,y=no2,group=depth))+
  geom_point(aes(color=depth),size=4)+
  geom_line(aes(color=depth,group=depth),size=1.5)+
  scale_color_manual(name = "Depth",
                     values=c("#FCCf14","#008082"),
                     labels=c("Surface", "Bottom"))+
  scale_fill_manual(name = "Depth",
                    values=c("#FCCf14","#008082"),
                    labels=c("Surface", "Bottom"))+
  facet_wrap(~depth,nrow = 2,strip.position = "right",
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  ylab("NO<sub>2</sub><sup>-</sub> (μmol L<sup>-1</sup>)")+
  xlab("\nMonth_Year")+
  theme(axis.title = element_markdown(size=11,face="bold"),
        #axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.position = "none")

nitrosopumulus<-dpp_pseq_rel1_melt%>%
  filter(Genus=="Candidatus Nitrosopumilus")%>%
  ggplot(aes(x=mm_yy,y=Abundance))+
  geom_bar(stat="identity")+
  facet_wrap(~depth,nrow = 2,strip.position = "right",
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  ylab("*Ca.* Nitrosopumilus Relative Abundance (%)")+
  xlab("\nMonth_Year")+
  theme(axis.title = element_markdown(size=11,face="bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1))
nitrosopumulus

amoA<-metadata_DPP%>%
  ggplot(aes(x=mm_yy,y=amoA))+
  geom_bar(stat="identity",aes(fill=depth))+
  geom_errorbar(aes(ymin=(amoA-amoA_SD),
                    ymax=(amoA+amoA_SD)),
                color="black",
                size=0.8,
                width=0.8)+
  scale_fill_manual(name = "Depth",
                    values=c("#FCCf14","#008082"),
                    labels=c("Surface", "Bottom"))+
  facet_wrap(~depth,nrow = 2,strip.position = "right",
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  ylab("amoA copies L<sup>-1</sup>")+
  xlab("\nMonth_Year")+
  theme(axis.title = element_markdown(size=11,face="bold"),
        #axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.position = "none")

amoA

ggarrange(nitrites,amoA,nitrosopumulus,nrow=3,
          labels = c("a","b","c"),
          label.x = c(0.05,0.05,0.05),
          align = "hv")

####Drawing Figure S4####
#relative abundance of slected taxa
dpp_pseq_rel1_melt%>%
  filter(Genus%in%c("SUP05 cluster",
                    "unassigned_SAR324 clade(Marine group B)",
                    "OM60(NOR5) clade","HIMB11",
                    "NS4 marine group","NS5 marine group"))%>%
  ggplot(aes(x=mm_yy,y=Abundance))+
  geom_bar(stat="identity",position = "stack",aes(fill=Genus))+
  facet_wrap(~depth,nrow = 2,strip.position = "right",
             labeller= as_labeller(c(S="Surface",B="Bottom")))+
  ylab("Relative Abundance (%)")+
  xlab("\nMonth_Year")+
  theme(axis.title = element_markdown(size=11,face="bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text=element_text(size=10),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,size=1),
        legend.position = "bottom")
