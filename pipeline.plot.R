library(ggplot2)
library(plyr)
library(forcats)
library(dplyr)
library(Rmisc)
library(VennDiagram)
library(stringr)
library(openxlsx)
library(ggfortify)
library(radmixture)
library(yaml)
library(kableExtra)
library(scales) # required for comma
library(dplyr)
library(nanopoRe) # devtools::install_github("sagrudd/nanopoRe") 
library(VariantAnnotation)
library(vcfR)
library(CePa)

###################### Sample Information ######################

library(openxlsx)
data_samples_info <- read.xlsx(paste('example','samples.xlsx',sep = '/'))

data_samples_info$Group <- 'HAN'
data_samples_info[data_samples_info$Races=='TIB',]$Group <- 'TIB'
data_samples_info[data_samples_info$Races=='Biaka',]$Group <- 'AFR'
data_samples_info[data_samples_info$Races=='HUM',]$Group <- 'HUM'
data_samples_info$Group <- factor(data_samples_info$Group, levels = c('HUM', 'AFR', 'HAN', 'TIB'))

data_samples_info$Group_Detail <- as.character(data_samples_info$Group)
data_samples_info[data_samples_info$Races=='HAN'&data_samples_info$Location%in%c('Guangxi','Guangdong'),]$Group_Detail <- 'HANS'
data_samples_info[data_samples_info$Races=='HAN'&!data_samples_info$Location%in%c('Guangxi','Guangdong'),]$Group_Detail <- 'HANN'
data_samples_info[data_samples_info$Altitude>=4000&data_samples_info$Races=='TIB',]$Group_Detail <- 'TIBG'
data_samples_info[data_samples_info$Altitude<4000&data_samples_info$Races=='TIB',]$Group_Detail <- 'TIBL'
data_samples_info[data_samples_info$Races=='Biaka',]$Group_Detail <- 'AFR'
data_samples_info[data_samples_info$Races=='HUM',]$Group_Detail <- 'HUM'
data_samples_info$Group_Detail <- factor(data_samples_info$Group_Detail, levels = c('HUM', 'AFR', 'HANN', 'HANS', 'TIBL', 'TIBG'))

data_samples_info$Source <- 'NA'
data_samples_info[str_detect(data_samples_info$SampleID, 'CQ'),]$Source <- 'T15H10'
data_samples_info[str_detect(data_samples_info$SampleID, 'SAM'),]$Source <- 'T38H39'
data_samples_info[str_detect(data_samples_info$SampleID, 'WGC'),]$Source <- 'T48H50'
data_samples_info[str_detect(data_samples_info$SampleID, 'HGDP'),]$Source <- 'HGDP'
data_samples_info[data_samples_info$Races=='HUM',]$Source <- 'HUM'
data_samples_info$Source <- factor(data_samples_info$Source)

data_samples_info[is.na(data_samples_info$Altitude),]$Altitude <- 0

###################### heterozygosity and homozygosity ######################

data_sample_het <- read.table(paste('example', 'merge.genotypes.corrected.delCHR.svtk.gte1ngs.LDpruned.het',sep = '/'), header = TRUE)
data_sample_het$PHt.LDpruned <- (data_sample_het$N_SITES-data_sample_het$O.HOM.)
data_sample_het$PHot.LDpruned <- data_sample_het$O.HOM.
data_sample_het <- data_sample_het[data_sample_het$INDV%in%data_samples_info$SampleID,]
data_samples_info$PHt.LDpruned <- NA
data_samples_info$PHot.LDpruned <- NA
data_samples_info[order(data_samples_info$SampleID),]$PHt.LDpruned <- data_sample_het[order(data_sample_het$INDV),]$PHt.LDpruned
data_samples_info[order(data_samples_info$SampleID),]$PHot.LDpruned <- data_sample_het[order(data_sample_het$INDV),]$PHot.LDpruned

data_sample_het_snp <- read.table(paste('example', 'merge.T48H50andT38H39.reserved.LDpruned.het',sep = '/'), header = TRUE)
data_sample_het_snp <- data_sample_het_snp[data_sample_het_snp$INDV%in%data_samples_info$SampleID,]
data_sample_het_snp$PHt.LDpruned.SNP <- (data_sample_het_snp$N_SITES-data_sample_het_snp$O.HOM.)/data_sample_het_snp$N_SITES
data_samples_info <- merge(data_samples_info,data_sample_het_snp[,c('INDV','PHt.LDpruned.SNP')],by.x='SampleID',by.y='INDV',all.x=TRUE)

###################### sample stats ######################
data_sv_dist <- read.table(paste('example', 'merge.genotypes.corrected.delCHR.svtk.dist', sep = '/'))
colnames(data_sv_dist) <- c('Dist','Chrom','All','Shared','Major','Polymorphic','Singleton')
data_sv_dist <- data_sv_dist[data_sv_dist$All>0, ]

###################### sv details ######################
source('scripts/funcs_load_data.R')
path_TIBvsHAN <- paste('example', 'merge.genotypes.corrected.delCHR.svtk.annotsv.tsv', sep = '/')
data_all <- func.load.svtk.annotsv(path = path_TIBvsHAN, samplelist = data_samples_info)

data_sv_details_all <- data_all$data_sv
data_sample_details <- data_all$data_samples
data_genotypes <- data_all$data_genotypes
rownames(data_genotypes) <- data_sv_details_all$SVID
data_genotypes <- data_genotypes[order(data_sv_details_all$SVID),]
data_sv_details_all <- data_sv_details_all[order(data_sv_details_all$SVID),]
data_sample_details_ONT <- subset(data_sample_details, data_sample_details$id%in%data_samples_info[data_samples_info$Platform%in%c('ONT'),]$SampleID)

###################### sv vcf details ######################
source('funcs_load_data.R')
path_TIBvsHAN_vcf <- paste(path_data, 'merge.genotypes.corrected.delCHR.svtk.vcf.gz', sep = '/')
data_all_vcf <- func.load.svtk.vcf.local(path_TIBvsHAN_vcf, reference_fasta=paste(path_data, '../genome/hg19.fa', sep = '/'))
data_sv_details_raw <- data_all_vcf$data_sv_vcf
data_sv_details_all_karyo <- data_all_vcf$data_sv_karyo


###################### confguration for plots ######################
prefix_filename = 'plots'
config_color_repeat <- c(
  'Low_complexity'='light blue',
  'Simple_repeat'='blue',
  'Satellite'='dark blue',
  'SINE'='green',
  'LINE'='yellow green',
  'LTR'='dark green',
  'DNA'='dark red',
  'RNA'='purple',
  'RC'='yellow',
  'None'='gray',
  'Other'='darkgray')

config_color_svtype <- c('DEL'='dark red',
                         'INS'='purple',
                         'DUP'='blue',
                         'INV'='dark green',
                         'TRA'='light blue')

config_color_group_supp <- c('Singleton'='light blue',
                             'Polymorphic'='blue',
                             'Major'='purple',
                             'Shared'='dark red')

config_color_group_pop <- c('All'='black',
                            'Shared'='gray',
                            'AFR'='dark red',
                            'HAN'='dark green',
                            'TIB'='dark blue',
                            'HANN'='green',
                            'HANS'='dark green',
                            'TIBG'='blue',
                            'TIBL'='dodger blue')

config_color_group_supp_han <- c('No Support'='gray',
                                 'Singleton'='light blue',
                                 'Polymorphic'='blue',
                                 'Major'='purple',
                                 'Shared'='dark red')

# config_color_group_len <- data.frame(GROUP_LEN=c('0-1k',    '1k-20k', 'above 20k'),
                                     # Color=c('light blue', 'purple','dark red'))
config_color_group_len <- c('0-100bp'='dark red',
                            '100bp-250bp'='purple',
                            '250bp-500bp'='blue',
                            '500bp-1kb'='light blue',
                            '50bp-1kb'='dark red',
                            '1kb-10kb'='purple',
                            '10kb-100kb'='blue',
                            '100kb-1Mb'='dark green',
                            'above 1Mb'='light blue')

config_color_database <- c('1000Genome'='dark red',
                           'DDD'='purple',
                           'DGV'='blue',
                           'GnomAD'='dark green',
                           'IMH'='light blue',
                           'None'='gray')

config_color_gene <- c('COPY_GAIN'='blue',
                       'DUP_PARTIAL'='purple',
                       'LOF'='yellow green',
                       'INV_SPAN'='dark red',
                       'UTR'='dark green',
                       'INTRONIC'='gray',
                       'INTERGENIC'='black')

val_res = 300
data_plot = data_sv_details_all

###################### SV locations (circos)  ######################
library(circlize)
karyo_plot <- as.data.frame(data_sv_details_all_karyo)
karyo_plot$SUPP <- as.integer(karyo_plot$SUPP)

karyo_plot_BND <- karyo_plot[karyo_plot$SVTYPE=='TRA',]
karyo_plot_BND$seqnames <- paste('chr',karyo_plot_BND$seqnames, sep = '')
karyo_plot_BND_link <- karyo_plot_BND[,c('CHR2','END','END')]
colnames(karyo_plot_BND_link) <- c('seqnames', 'start', 'end')
karyo_plot_BND_link$seqnames <- paste('chr', as.character(karyo_plot_BND_link$seqnames), sep = "")
karyo_plot_BND_link$start <- as.integer(karyo_plot_BND_link$start)
karyo_plot_BND_link$end <- as.integer(karyo_plot_BND_link$end)

array_seqnames <- paste('chr', karyo_plot$seqnames,sep = '')

pdf(file = paste(prefix_filename,'sv.circos','pdf',sep = '.'), width = 8, height = 8)
  circos.initializeWithIdeogram(species = 'hg19')
  
  group_shared <- karyo_plot$SUPP==25
  group_major <- karyo_plot$SUPP%in%seq(13,24)
  group_polymorphic <- karyo_plot$SUPP%in%seq(2,12)
  group_singleton <- karyo_plot$SUPP==1
  
  circos.trackHist(factors=array_seqnames, 
                   track.height = 0.1,
                   x=karyo_plot$start,col = "#999999", border = '#999999', 
                   bg.border = NA, bin.size = 500000)
  circos.trackHist(factors=array_seqnames[group_shared], 
                   track.height = 0.1, 
                   x=karyo_plot[group_shared,]$start,col = "red", border = 'red', 
                   bg.border = NA,bin.size = 500000)
  circos.trackHist(factors=array_seqnames[group_major], track.height = 0.1, x=karyo_plot[group_major,]$start,col = "purple", border = "purple",bg.border = NA,bin.size = 500000)
  circos.trackHist(factors=array_seqnames[group_polymorphic], track.height = 0.1, x=karyo_plot[group_polymorphic,]$start,col = "blue",border = "blue",bg.border = NA,bin.size = 500000)
  circos.trackHist(factors=array_seqnames[group_singleton], track.height = 0.1, x=karyo_plot[group_singleton,]$start,col = "light blue",border = "light blue",bg.border = NA,bin.size = 500000)
  
  group_shared <- karyo_plot_BND$SUPP==25
  group_major <- karyo_plot_BND$SUPP%in%seq(13,24)
  group_polymorphic <- karyo_plot_BND$SUPP%in%seq(2,12)
  group_singleton <- karyo_plot_BND$SUPP==1
dev.off()

###################### samples NGS location ######################
library(maps)
library(mapdata)
library(maptools);
china_map=readShapePoly('china-province-border-data/bou2_4p.shp');
china_map@data$cName <- iconv(china_map@data$NAME, from = "GBK")

x <- china_map@data                          #读取行政信息
xs <- data.frame(x,id=seq(0:924)-1)          #含岛屿共925个形状
china_map1 <- fortify(china_map)             #转化为数据框
china_map_data <- join(china_map1, xs, type = "full", )       #合并两个数据框

tmp <- summary(factor(data_samples_info[data_samples_info$Platform=='NGS',]$cLocation))

NAME <- names(tmp)
pop <- tmp
pop <- data.frame(NAME, pop)
colnames(pop) <- c('cName', 'pop')
china_map_pop <- join(china_map_data, pop, type = "full")

pdf(file = paste(prefix_filename,'samples.map','pdf',sep = '/'), width = 8, height = 6)
postscript(file = paste(prefix_filename,'samples.map','eps',sep = '.'), width = 8, height = 6)
ggplot(china_map_pop, aes(x = long, y = lat, group = group, fill = pop)) +
  geom_polygon() +
  geom_path(color = "grey40") +
  coord_map() +
  xlab('') + 
  ylab("") + 
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 8,face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
dev.off()

###################### samples individual ######################
##need depackage scater
df_sorted <- arrange(data_sample_details_ONT, id, Type) 
df_cumsum <- ddply(df_sorted, "id",
                   transform, 
                   ypos=cumsum(Discovery) - 0.5*Discovery)

pdf(file = paste(prefix_filename,'sv.samples.bar','pdf',sep = '/'), width = 4, height = 3)
ggplot(data=df_cumsum, aes(x=id, y=Discovery, fill = Type)) +
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.direction = 'horizontal',
        legend.position=c(0.5,0.98),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold", angle = 60, vjust = 0.6),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=c('light blue', 'blue', 'purple','dark red')) + 
  xlab("") + 
  ylab("Discovery")
dev.off()

###################### telomere enrichment ######################

pdf(file = paste(prefix_filename,'Svs_dist.gt.update','pdf',sep = '/'), width = 4, height = 3)
ggplot(data = data_sv_dist[data_sv_dist$All!=0,],
       aes(x = Dist, y = All)) +
  geom_point(size = 0.05, color="blue") +
  annotate("rect", fill = "dark gray", alpha = 0.5, 
           xmin = 0, xmax = 5,
           ymin = -Inf, ymax = Inf)+
  xlim(0, 150) +
  # ylim(0, 250) +
  ylab('SVs Per 500kbp') +
  xlab('Telomere Distance') +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 6,face = "bold"),
    axis.text.y = element_text(size = 6,face = "bold"),
    axis.title.y = element_text(size = 6,face = "bold"),
    axis.title.x = element_text(size = 6,face = "bold"))
dev.off()

###################### SV & GROUP_SUPP ######################
data_sv_group <- plyr::count(data_plot,'GROUP_SUPP')
data_sv_group$GROUP_SUPP <- factor(data_sv_group$GROUP_SUPP)
pdf(file = paste(prefix_filename,'sv.group_supp.pie','pdf',sep = '/'), width = 3, height = 3)
ggplot(data_sv_group, aes(x="",y=freq,fill=GROUP_SUPP)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/4+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_sv_group$freq/sum(data_sv_group$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
scale_fill_manual(values=config_color_group_supp) +
  theme_minimal() + 
  theme(legend.position = c(0.5,0),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.title=element_blank(),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

###################### SVLEN ######################
pdf(file = paste(prefix_filename,'svlen.freqploy','pdf',sep = '、'), width = 4, height = 3)
ggplot(data=data_plot[data_plot$SVTYPE!='TRA',], aes(x=SVLEN, color = SVTYPE)) +
  geom_freqpoly(binwidth=1/5, size=0.8) +  
  scale_x_continuous(trans = 'log10',
                     breaks=c(100,1000,10000,100000,1000000,10000000),
                     labels =c('100bp','1kb','10kb','100kb','1Mb', '10Mb'),
                     limits = c(50, 10000000)
  ) + 
  scale_y_continuous(trans = 'log10') + 
  theme_minimal() +
  theme(legend.title=element_blank(),
        legend.position=c(0.9,0.9),
        legend.direction = 'vertical',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_svtype) + 
  xlab('') + 
  ylab("SV Count")
dev.off()

###################### SVLEN & GROUP_SUPP ######################
data_sv_group <- plyr::count(data_plot,c('GROUP_LEN','GROUP_SUPP'))
data_sv_group$label <- 0

for (type_sv in levels(data_sv_group$GROUP_LEN)){
  data_sv_group[data_sv_group$GROUP_LEN==type_sv,]$label <- sum(data_sv_group[data_sv_group$GROUP_LEN==type_sv,]$freq)
}
pdf(file = paste(prefix_filename,'svlen.group_supp.bar','pdf',sep = '.'), width = 3, height = 3)
ggplot(data_sv_group, aes(x = GROUP_LEN, y = freq, fill = GROUP_SUPP)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=1)
dev.off()

###################### SVTYPE & GROUP_SUPP ######################
data_sv_group <- plyr::count(data_plot,c('SVTYPE','GROUP_SUPP'))
data_sv_group$SVTYPE <- factor(data_sv_group$SVTYPE)
data_sv_group$label <- 0

for (type_sv in levels(data_sv_group$SVTYPE)){
  data_sv_group[data_sv_group$SVTYPE==type_sv,]$label <- sum(data_sv_group[data_sv_group$SVTYPE==type_sv,]$freq)
}
data_sv_group$SVTYPE <- factor(data_sv_group$SVTYPE, levels = c( 'DEL', 'INS', 'DUP', 'INV', 'TRA'))

pdf(file = paste(prefix_filename,'svtype.group_supp.bar','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_sv_group, aes(x = SVTYPE, y = freq, fill = GROUP_SUPP)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2)
dev.off()

###################### Brepoints CI ###################################
data_tmp <- data_plot 

data_tmp$GROUP_CI <- NA
data_tmp[data_tmp$MAXCI_POS<=1000&data_tmp$MAXCI_END<=1000,]$GROUP_CI <- '500bp-1kb'
data_tmp[data_tmp$MAXCI_POS<=500&data_tmp$MAXCI_END<=500,]$GROUP_CI <- '250bp-500bp'
data_tmp[data_tmp$MAXCI_POS<=250&data_tmp$MAXCI_END<=250,]$GROUP_CI <- '100bp-250bp'
data_tmp[data_tmp$MAXCI_POS<=100&data_tmp$MAXCI_END<=100,]$GROUP_CI <- '0-100bp'
data_tmp$GROUP_CI <- factor(data_tmp$GROUP_CI, levels = c('0-100bp', '100bp-250bp', '250bp-500bp', '500bp-1kb'))


df.new<-ddply(data_tmp,.(GROUP_SUPP),plyr::summarise,
              prop=prop.table(table(GROUP_CI)),
              SUPP=names(table(GROUP_CI)))
df.new$SUPP <- factor(df.new$SUPP, levels = c('0-100bp', '100bp-250bp', '250bp-500bp', '500bp-1kb'))
pdf(file = paste(prefix_filename,'breakpoints.ci.group_supp.prop','pdf',sep = '.'), width = 4, height =3)
ggplot(df.new, aes(SUPP, prop, fill=GROUP_SUPP)) + 
  geom_bar(stat="identity",position = 'dodge') + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.5,0.9),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold")
        ) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("Max Interval of Breakpoints") + 
  ylab("Prop")
dev.off()

###################### Brepoints & Repeat ###################################
data_tmp <- data_plot[,c('SVID','SVTYPE','GROUP_LEN','SVLEN','Repeats_type_right', 'Repeats_type_left', 'GCcontent_left','GCcontent_right', 'GROUP_SUPP')]
data_tmp <- rbind(data_tmp, data_tmp)
data_tmp$RepClass <- factor(c(as.character(data_plot$Repeats_type_left),as.character(data_plot$Repeats_type_right)), levels = c('Low_complexity', 'Simple_repeat', 'Satellite','SINE','LINE','LTR','DNA','RNA','RC','Other','None'))
data_tmp$GCcontent <- c(data_plot$GCcontent_left,data_plot$GCcontent_right)

data_tmp <- data_tmp[data_tmp$RepClass!='Other',]
data_tmp$RepClass <- factor(as.character(data_tmp$RepClass),levels = c('Low_complexity', 'Simple_repeat', 'Satellite','SINE','LINE','LTR','DNA','RNA','RC','None'))

data_sv_group <- plyr::count(data_tmp,c('RepClass','GROUP_SUPP'))
data_sv_group$label <- 0

for (type_sv in levels(data_sv_group$RepClass)){
  data_sv_group[data_sv_group$RepClass==type_sv,]$label <- sum(data_sv_group[data_sv_group$RepClass==type_sv,]$freq)
}

# tiff(filename = paste(prefix_filename,'sv_summary.tiff',sep = '.'), width = 1200, height = 1000, res = 300)
# png(filename = paste(prefix_filename,'breakpoints.repeat.group_supp.bar.png',sep = '.'), width = 1500, height = 2000, res = val_res)
pdf(file = paste(prefix_filename,'breakpoints.repeat.group_supp.bar','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_group, aes(x = RepClass, y = freq, fill = GROUP_SUPP)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2.5)
dev.off()

data_sv_repeat <- subset(data_tmp, SVTYPE%in%c('DEL','INS')&GROUP_LEN%in%c('50bp-1kb','1kb-10kb'))
# data_sv_repeat <- subset(data_tmp, SVTYPE%in%c('DEL','INS')&SVLEN<=1000000&GROUP_SUPP%in%c('Major','Shared'))
data_sv_repeat_ins <- subset(data_sv_repeat, SVTYPE=='INS')
data_sv_repeat_del <- subset(data_sv_repeat, SVTYPE=='DEL')


data_sv_group <- plyr::count(data_sv_repeat,'RepClass')
# png(filename = paste(prefix_filename,'breakpoints.repeat.pie.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'breakpoints.repeat.pie','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_sv_group, aes(x="",y=freq,fill=RepClass)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/10+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_sv_group$freq/sum(data_sv_group$freq)),
            size=2,
            color='white', 
            fontface="bold")+
  scale_fill_manual(values=config_color_repeat) +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

data_sv_group <- plyr::count(data_sv_repeat_ins,'RepClass')
# png(filename = paste(prefix_filename,'breakpoints.repeat.ins.pie.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'breakpoints.repeat.ins.pie','pdf',sep = '.'), width = 4, height = 3)
  ggplot(data_sv_group, aes(x="",y=freq,fill=RepClass)) + 
    geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
    scale_fill_manual(values=config_color_repeat) +
    theme_minimal() + 
    theme(legend.title=element_blank(),
          legend.position="bottom",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.key.size=unit(0.5, 'cm'),
          legend.text = element_text(size = 6,face = "bold"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
dev.off()

data_sv_group <- plyr::count(data_sv_repeat_del,'RepClass')
pdf(file = paste(prefix_filename,'breakpoints.repeat.del.pie','pdf',sep = '.'), width = 4, height = 3)
# png(filename = paste(prefix_filename,'breakpoints.repeat.del.pie.png',sep = '.'), width = 1200, height = 1000, res = val_res)
  ggplot(data_sv_group, aes(x="",y=freq,fill=RepClass)) + 
    geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
    scale_fill_manual(values=config_color_repeat) +
    theme_minimal() + 
    theme(legend.title=element_blank(),
          legend.position="bottom",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.key.size=unit(0.5, 'cm'),
          legend.text = element_text(size = 6,face = "bold"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
dev.off()

###################### Brepoints & Repeat & formation ###################################
data_repeat_set <- data_plot[,c('SVTYPE','GROUP_LEN','SVLEN','Repeats_type_right', 'Repeats_type_left', 'GCcontent_left','GCcontent_right', 'GROUP_SUPP', 'RepClass')]
data_repeat_set <- data_repeat_set[data_repeat_set$SVLEN<=1000000,]
data_repeat_set_del <- data_repeat_set[data_repeat_set$SVTYPE=='DEL',]
data_repeat_set_ins <- data_repeat_set[data_repeat_set$SVTYPE=='INS',]

num_sv <- nrow(data_repeat_set_del)
condition_vntr <- data_repeat_set_del$RepClass%in%c('Low_complexity','Simple_repeat','Satellite')
condition_nahr <- (!condition_vntr)&data_repeat_set_del$Repeats_type_right==data_repeat_set_del$Repeats_type_left&data_repeat_set_del$Repeats_type_left!='None'
condition_te <- !condition_vntr & ! condition_nahr & 
  (data_repeat_set_del$Repeats_type_right%in%
     c('SINE','LINE','LTR','DNA')|
     data_repeat_set_del$Repeats_type_left%in%
     c('SINE','LINE','LTR','DNA'))

# data_rep_tmp <- data.frame(mech=c('VNTR','NAHR','TE','Others'),
#                            freq=c(
#                              sum(condition_vntr),
#                              sum(condition_nahr),
#                              sum(condition_te),
#                              num_sv -  sum(condition_vntr) - sum(condition_nahr) - sum(condition_te)))


ids_pie = c('SV','DEL', 'INS')
labels_pie = c('SV formation','DEL', 'INS')
parents_pie = c('','SV','SV')
values_pie = c(38216,nrow(data_repeat_set_del), nrow(data_repeat_set_ins))

ids_pie <- c(ids_pie, 'VNTR-DEL')
labels_pie <- c(labels_pie, 'VNTR')
parents_pie <- c(parents_pie, 'DEL')
values_pie <- c(values_pie, sum(condition_vntr))

for (label_condition in c('Low_complexity','Simple_repeat','Satellite')) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'VNTR-DEL', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'VNTR-DEL')
  values_pie <- c(values_pie, sum(data_repeat_set_del[condition_vntr,]$RepClass==label_condition) )
}

ids_pie <- c(ids_pie, 'NAHR')
labels_pie <- c(labels_pie, 'NAHR')
parents_pie <- c(parents_pie, 'DEL')
values_pie <- c(values_pie, sum(condition_nahr))

for (label_condition in c('SINE','LINE','LTR', 'Other')) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'NAHR', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'NAHR')
  if (label_condition == 'Other'){
    values_pie <- c(values_pie, sum(!data_repeat_set_del[condition_nahr,]$Repeats_type_left%in%c('SINE','LINE','LTR')) )
  } else {
    values_pie <- c(values_pie, sum(data_repeat_set_del[condition_nahr,]$Repeats_type_left==label_condition) )
  }
}


ids_pie <- c(ids_pie, 'TE-DEL')
labels_pie <- c(labels_pie, 'TE')
parents_pie <- c(parents_pie, 'DEL')
values_pie <- c(values_pie, sum(condition_te))

for (label_condition in c('SINE','LINE', 'Other')) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'TE-DEL', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'TE-DEL')
  
  if (label_condition == 'Other'){
    values_pie <- c(values_pie, 
                    sum(condition_te)-
                    sum((data_repeat_set_del[condition_te,]$Repeats_type_left=='SINE'|
                           data_repeat_set_del[condition_te,]$Repeats_type_right=='SINE') &
                          data_repeat_set_del[condition_te,]$Repeats_type_left!='LINE' &
                          data_repeat_set_del[condition_te,]$Repeats_type_right!='LINE')-
                    sum((data_repeat_set_del[condition_te,]$Repeats_type_left=='LINE'|
                            data_repeat_set_del[condition_te,]$Repeats_type_right=='LINE') &
                          data_repeat_set_del[condition_te,]$Repeats_type_left!='SINE' &
                          data_repeat_set_del[condition_te,]$Repeats_type_right!='SINE'))
                    
  } else {
    values_pie <- c(values_pie, 
                    sum((data_repeat_set_del[condition_te,]$Repeats_type_left==label_condition|
                          data_repeat_set_del[condition_te,]$Repeats_type_right==label_condition) &
                          !data_repeat_set_del[condition_te,]$Repeats_type_left%in%setdiff(c('SINE','LINE'), label_condition) &
                          !data_repeat_set_del[condition_te,]$Repeats_type_right%in%setdiff(c('SINE','LINE'), label_condition)))
                   
  }
}

# ids_pie <- c(ids_pie, 'Other-DEL')
# labels_pie <- c(labels_pie, 'Others')
# parents_pie <- c(parents_pie, 'DEL')
# values_pie <- c(values_pie,  num_sv -  sum(condition_vntr) - sum(condition_nahr) - sum(condition_te))
# 

num_sv <- nrow(data_repeat_set_ins)
condition_vntr <- data_repeat_set_ins$RepClass%in%c('Low_complexity','Simple_repeat','Satellite')
# condition_nahr <- !condition_vntr&data_repeat_set_ins$Repeats_type_right==data_repeat_set_ins$Repeats_type_left
condition_te <- !condition_vntr & (data_repeat_set_ins$Repeats_type_right%in%c('SINE','LINE','LTR','DNA')|data_repeat_set_ins$Repeats_type_left%in%c('SINE','LINE','LTR','DNA'))

ids_pie <- c(ids_pie, 'VNTR-INS')
labels_pie <- c(labels_pie, 'VNTR')
parents_pie <- c(parents_pie, 'INS')
values_pie <- c(values_pie, sum(condition_vntr))

for (label_condition in c('Low_complexity','Simple_repeat','Satellite')) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'VNTR-INS', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'VNTR-INS')
  values_pie <- c(values_pie, sum(data_repeat_set_ins[condition_vntr,]$RepClass==label_condition) )
}

ids_pie <- c(ids_pie, 'TE-INS')
labels_pie <- c(labels_pie, 'TE')
parents_pie <- c(parents_pie, 'INS')
values_pie <- c(values_pie, sum(condition_te))

for (label_condition in c('SINE','LINE', 'Other')) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'TE-INS', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'TE-INS')
  
  if (label_condition == 'Other'){
    values_pie <- c(values_pie, 
                    sum(condition_te)-
                      sum((data_repeat_set_ins[condition_te,]$Repeats_type_left=='SINE'|
                             data_repeat_set_ins[condition_te,]$Repeats_type_right=='SINE') &
                            data_repeat_set_ins[condition_te,]$Repeats_type_left!='LINE' &
                            data_repeat_set_ins[condition_te,]$Repeats_type_right!='LINE')-
                      sum((data_repeat_set_ins[condition_te,]$Repeats_type_left=='LINE'|
                             data_repeat_set_ins[condition_te,]$Repeats_type_right=='LINE') &
                            data_repeat_set_ins[condition_te,]$Repeats_type_left!='SINE' &
                            data_repeat_set_ins[condition_te,]$Repeats_type_right!='SINE'))
    
  } else {
    values_pie <- c(values_pie, 
                    sum((data_repeat_set_ins[condition_te,]$Repeats_type_left==label_condition|
                           data_repeat_set_ins[condition_te,]$Repeats_type_right==label_condition) &
                          !data_repeat_set_ins[condition_te,]$Repeats_type_left%in%setdiff(c('SINE','LINE'), label_condition) &
                          !data_repeat_set_ins[condition_te,]$Repeats_type_right%in%setdiff(c('SINE','LINE'), label_condition)))
    
  }
}

if (!require("processx")) install.packages("processx")
fig <- plot_ly(
              ids=ids_pie,
              labels=labels_pie, 
              parents=parents_pie,
              values=values_pie,
              type='sunburst',
              branchvalues = 'total'
              ) %>% add_surface()
orca(fig, "../data/result_genotypted/plots/surface-plot.svg")

# png(filename = paste(prefix_filename,'formation.pie.del.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'formation.pie.del','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_rep_tmp, aes(x="",y=freq,fill=mech)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/4+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_rep_tmp$freq/sum(data_rep_tmp$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

data_rep_tmp <- plyr::count(data_repeat_set_del[condition_nahr,],
                            'Repeats_type_left')
data_tmp <- data_repeat_set_del[condition_nahr,]
data_tmp[data_tmp$Repeats_type_left%in%data_rep_tmp[data_rep_tmp$freq/sum(data_rep_tmp$freq)<0.05,]$Repeats_type_left,]$Repeats_type_left <- 'Other'
data_rep_tmp <- plyr::count(data_tmp,
                            'Repeats_type_left')


pdf(file = paste(prefix_filename,'formation.pie.del.nahr','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_rep_tmp, aes(x="",y=freq,fill=Repeats_type_left)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/10+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_rep_tmp$freq/sum(data_rep_tmp$freq)),
            size=2,
            color='white', 
            fontface="bold")+
  scale_fill_manual(values=config_color_repeat) +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()


data_rep_tmp <- plyr::count(data_repeat_set_del[condition_vntr,],
                            'RepClass')

pdf(file = paste(prefix_filename,'formation.pie.del.vntr','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_rep_tmp, aes(x="",y=freq,fill=RepClass)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/10+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_rep_tmp$freq/sum(data_rep_tmp$freq)),
            size=2,
            color='white', 
            fontface="bold")+
  scale_fill_manual(values=config_color_repeat) +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

num_sv <- nrow(data_repeat_set_ins)
condition_vntr <- data_repeat_set_ins$RepClass%in%c('Low_complexity','Simple_repeat','Satellite')
# condition_nahr <- !condition_vntr&data_repeat_set_ins$Repeats_type_right==data_repeat_set_ins$Repeats_type_left
condition_te <- !condition_vntr & (data_repeat_set_ins$Repeats_type_right%in%c('SINE','LINE','LTR','DNA')|data_repeat_set_ins$Repeats_type_left%in%c('SINE','LINE','LTR','DNA'))

data_rep_tmp <- data.frame(mech=c('VNTR','TE','Others'),
                           freq=c(
                             sum(condition_vntr),
                             # sum(condition_nahr),
                             sum(condition_te),
                             num_sv -  sum(condition_vntr) - sum(condition_te)))

# png(filename = paste(prefix_filename,'formation.pie.ins.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'formation.pie.ins','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_rep_tmp, aes(x="",y=freq,fill=mech)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/4+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_rep_tmp$freq/sum(data_rep_tmp$freq)),
            size=2,
            color='white', 
            fontface="bold") +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

###################### Brepoints & Repeat & SVLEN ###################################

pdf(file = paste(prefix_filename,'breakpoints.repeat.len.lte1m','pdf',sep = '.'), 
    width = 3, height = 2.5)
ggplot(data=data_sv_repeat, aes(x=abs(SVLEN), fill = RepClass)) +
  geom_histogram() + 
  scale_x_continuous(trans = 'log10',breaks = c(100,300,2000,6000,20000),
                     labels = c('100bp','300bp','2kb','6kb','20kb'),limits = c(50,20000)) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.65),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold", angle = 45),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold")) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('') + 
  ylab("SV Count")
dev.off()

# tiff(filename = paste(prefix_filename,'sv_repeat_len_INS_lte1k.tiff',sep = '.'), width = 1200, height = 1000, res = 300)

# png(filename = paste(prefix_filename,'breakpoints.repeat.len.lte1m.ins.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'breakpoints.repeat.len.lte1m.ins','pdf',sep = '.'), width = 8, height = 6)
ggplot(data=data_sv_repeat_ins, aes(x=abs(SVLEN), fill = RepClass)) +
  geom_histogram() + 
  scale_x_continuous(trans = 'log10',breaks = c(100,300,2000,6000,20000),
                     labels = c('100bp','300bp','2kb','6kb','20kb'),limits = c(50,20000)) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.75),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold", angle = 45),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('') + 
  ylab("SV Count")
dev.off()

# tiff(filename = paste(prefix_filename,'sv_repeat_len_DEL_lte1k.tiff',sep = '.'), width = 1200, height = 1000, res = 300)
# png(filename = paste(prefix_filename,'breakpoints.repeat.len.lte1m.del.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'breakpoints.repeat.len.lte1m.del','pdf',sep = '.'), width = 8, height = 6)

ggplot(data=data_sv_repeat_del, aes(x=abs(SVLEN), fill = RepClass)) +
  geom_histogram() +
  scale_x_continuous(trans = 'log10',breaks = c(100,300,2000,6000,20000),
                     labels = c('100bp','300bp','2kb','6kb', '20kb'),limits = c(50,20000)) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.75),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold", angle = 45),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('') + 
  ylab("SV Count")
dev.off()

###################### Brepoints & GCcontent vs Reference ###################################


# library(gridExtra)
# png(filename = paste(prefix_filename,'breakpoints.repeat.enrichment.png',sep = '.'), width = 1200, height = 2000, res = val_res)
# # grid.arrange(plot.enrichment.ins,plot.enrichment.del,nrow=2)
# dev.off()

###################### Brepoints & GCcontent & Repeat###################################
pdf(file = paste(prefix_filename,'breakpoints.gc.repeat.ins','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'breakpoints.gc.repeat.ins.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data=data_sv_repeat_ins, aes(x=GCcontent, fill = RepClass)) +
  geom_histogram() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold", angle = 45),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold")) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('GC Content (INS)') + 
  ylab("SV Count")
dev.off()

pdf(file = paste(prefix_filename,'breakpoints.gc.repeat.del','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'breakpoints.gc.repeat.del.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data=data_sv_repeat_del, aes(x=GCcontent, fill = RepClass)) +
  geom_histogram() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold", angle = 45),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold")) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('GC Content (DEL)') + 
  ylab("SV Count")
dev.off()

pdf(file = paste(prefix_filename,'breakpoints.gc.repeat','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'breakpoints.gc.repeat.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data=data_sv_repeat, aes(x=GCcontent, fill = RepClass)) +
  geom_histogram() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold", angle = 45),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold")) +
  scale_fill_manual(values=config_color_repeat) + 
  xlab('GC Content') + 
  ylab("SV Count")
dev.off()

###################### Gene Annotation####################################################

c('COPY_GAIN','DUP_PARTIAL','LOF','INV_SPAN','UTR','INTRONIC','INTERGENIC')
data_sv_group <- plyr::count(data_plot[data_plot$SVLEN<1000000,],'GROUP_SV')
data_sv_group$GROUP_SV <- factor(data_sv_group$GROUP_SV ,levels = c('DUP_PARTIAL','LOF','INV_SPAN','UTR','INTRONIC','INTERGENIC'))
pdf(file = paste(prefix_filename,'gene.pie','pdf',sep = '.'), width = 2.5, height = 2.5)
ggplot(data_sv_group, aes(x="",y=freq,fill=GROUP_SV)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  scale_fill_manual(values=config_color_gene) +
  geom_text(aes(y=freq/6+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_sv_group$freq/sum(data_sv_group$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position=c(0.6,0.025),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

# png(filename = paste(prefix_filename,'gene.cds.pie.png',sep = '.'), width = 1200, height = 1000, res = val_res)
data_sv_group <- data_sv_group[!data_sv_group$GROUP_SV%in%c('INTRONIC','INTERGENIC'),]
pdf(file = paste(prefix_filename,'gene.cds.pie','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_sv_group, aes(x="",y=freq,fill=GROUP_SV)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  scale_fill_manual(values=config_color_gene) +
  geom_text(aes(y=freq/4+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_sv_group$freq/sum(data_sv_group$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()


data_sv_group <- plyr::count(data_plot[data_plot$SVLEN<1000000,],c('GROUP_SV','GROUP_SUPP'))
data_sv_group$label <- 0
data_sv_group$GROUP_SV <- factor(data_sv_group$GROUP_SV)

for (type_sv in levels(data_sv_group$GROUP_SV)){
  data_sv_group[data_sv_group$GROUP_SV==type_sv,]$label <- sum(data_sv_group[data_sv_group$GROUP_SV==type_sv,]$freq)
}

pdf(file = paste(prefix_filename,'gene.group_supp.bar','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'gene.group_supp.bar.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data_sv_group, aes(x = GROUP_SV, y = freq, fill = GROUP_SUPP)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2.5)
dev.off()




###################### Database ######################
data_sv_group <- plyr::count(data_plot,'Database')
pdf(file = paste(prefix_filename,'database.pie','pdf',sep = '.'), width = 3, height = 3)
ggplot(data_sv_group, aes(x="",y=freq,fill=Database)) + 
  geom_bar(width=1,stat="identity") + coord_polar("y",start=0) + 
  geom_text(aes(y=freq/6+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_sv_group$freq/sum(data_sv_group$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
  scale_fill_manual(values=config_color_database) +
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position=c(0.5,0),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

###################### Database & GROUP_SUPP ######################
data_sv_anno_group <- plyr::count(data_plot, c('GROUP_SUPP', 'Database'))
data_sv_anno_group$label <- 0

for (type_sv in levels(data_sv_anno_group$GROUP_SUPP)){
  data_sv_anno_group[data_sv_anno_group$GROUP_SUPP==type_sv,]$label <- sum(data_sv_anno_group[data_sv_anno_group$GROUP_SUPP==type_sv,]$freq)
}

pdf(file = paste(prefix_filename,'database.group_supp','pdf',sep = '.'), width = 4, height = 3)
ggplot(data_sv_anno_group, aes(x = GROUP_SUPP, y = freq, fill = forcats::fct_rev(Database))) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-20,0,0,0))) +
  scale_fill_manual(values=config_color_database) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2, fontface="bold")
dev.off()

###################### Database & SVTYPE ######################
data_sv_anno_group <- plyr::count(data_plot , c('SVTYPE', 'Database'))
data_sv_anno_group$SVTYPE <- factor(data_sv_anno_group$SVTYPE)
data_sv_anno_group$label <- 0

for (type_sv in levels(data_sv_anno_group$SVTYPE)){
  data_sv_anno_group[data_sv_anno_group$SVTYPE==type_sv,]$label <- sum(data_sv_anno_group[data_sv_anno_group$SVTYPE==type_sv,]$freq)
}

# tiff(filename = paste(prefix_filename,'sv_anno_group_svtype.tiff',sep = '.'), width = 1200, height = 1000, res = 300)
# png(filename = paste(prefix_filename,'database.svtype.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'database.svtype','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_anno_group, aes(x = SVTYPE, y = freq, fill = forcats::fct_rev(Database))) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-20,0,0,0))) +
  scale_fill_manual(values=config_color_database) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2, fontface="bold")
dev.off()

###################### Database & AF & GROUP_SUPP ######################
# tmp <- subset(data_plot, AF_Database > 0)
# tiff(filename = paste(prefix_filename,'sv_anno_af.tiff',sep = '.'), width = 1200, height = 1000, res = 300)
# png(filename = paste(prefix_filename,'database.af.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'database.af','pdf',sep = '.'), width = 8, height = 12)
ggplot(data=data_plot, aes(x=AF_Database, fill = GROUP_SUPP)) +
  geom_histogram() +
  theme_minimal()+
  # scale_y_continuous(trans = 'log', breaks = c(1000,10000),labels = c('1000','10000')) +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold", angle = 45),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold")) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab('AF') + 
  ylab("Discovery")
dev.off()



###################### Genotyping ######################
data_tmp <- data_plot
data_tmp$GROUP_SUPP_NGS <- as.character(data_tmp$GROUP_SUPP_NGS)
data_tmp <- data_tmp[!is.na(data_tmp$MR_NGS),]
# data_tmp <- data_tmp[!is.na(data_tmp$MR_NGS)&data_tmp$MR_NGS<0.05,]

data_tmp[is.na(data_tmp$GROUP_SUPP_NGS),]$GROUP_SUPP_NGS <- 'No Support'
data_tmp$GROUP_SUPP_NGS <- factor(data_tmp$GROUP_SUPP_NGS, levels = c('No Support','Singleton', 'Polymorphic', 'Major', 'Shared'))
data_tmp_1 <- plyr::count(data_tmp[data_tmp$MR_NGS<0.05,],'GROUP_SUPP_NGS')
data_tmp_2 <- plyr::count(data_tmp[data_tmp$MR_NGS<0.05&data_tmp$GROUP_SUPP_NGS=='No Support',],
                          'GROUP_SUPP')
data_tmp_3 <- plyr::count(data_tmp[data_tmp$MR_NGS<0.05&data_tmp$GROUP_SUPP_NGS!='No Support',],
                          'GROUP_SUPP')

ids_pie = c('SV', 'MR_0', 'MR_1', 'No_Support', 'Support')
labels_pie = c('SV genotyping', 'MR < 0.05', 'MR >= 0.05', 'genotyped AF = 0', 'genotyped AF > 0')
parents_pie = c('', 'SV', 'SV', 'MR_0', 'MR_0')
values_pie = c(nrow(data_tmp),
               sum(data_tmp$MR_NGS<0.05),
               sum(data_tmp$MR_NGS>=0.05),
               data_tmp_1[data_tmp_1$GROUP_SUPP_NGS=='No Support',]$freq,
               sum(data_tmp_1$freq)-data_tmp_1[data_tmp_1$GROUP_SUPP_NGS=='No Support',]$freq)

for (label_condition in data_tmp_2$GROUP_SUPP) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'No_Support', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'No_Support')
  values_pie <- c(values_pie, data_tmp_2[data_tmp_2$GROUP_SUPP==label_condition,]$freq)
}

for (label_condition in data_tmp_3$GROUP_SUPP) {
  
  ids_pie <- c(ids_pie, paste(label_condition, 'Support', sep = '_'))
  labels_pie <- c(labels_pie, label_condition)
  parents_pie <- c(parents_pie, 'Support')
  values_pie <- c(values_pie, data_tmp_3[data_tmp_3$GROUP_SUPP==label_condition,]$freq)
}

# data_sv_group <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & 
#                           MR_NGS <= 0.05 & 
#                           GROUP_SUPP%in%c('Shared', 'Major'))
data_sv_group <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & MR_NGS <= 0.05)

data_sv_group$GROUP_REP <- 'None'
data_sv_group[data_sv_group$RepClass!='None',]$GROUP_REP <- 'Repclass'
data_sv_group[data_sv_group$RepClass=='None'&data_sv_group$SD!='',]$GROUP_REP <- 'SD'
data_sv_group[data_sv_group$GROUP_REP=='None'&
                (data_sv_group$Repeats_type_left!='None'|data_sv_group$Repeats_type_right!='None'),]$GROUP_REP <- 'flanked'

data_tmp_4 <- plyr::count(data_sv_group,c('GROUP_SUPP','GROUP_REP'))

for (label_group in data_tmp_2$GROUP_SUPP){
  for (label_condition in unique(data_tmp_4$GROUP_REP)) {
    
    ids_pie <- c(ids_pie, paste(label_condition, label_group, 'No_Support', sep = '_'))
    labels_pie <- c(labels_pie, label_condition)
    parents_pie <- c(parents_pie, paste(label_group,'No_Support', sep = '_'))
    values_pie <- c(values_pie, data_tmp_4[data_tmp_4$GROUP_SUPP==label_group&
                                             data_tmp_4$GROUP_REP==label_condition,]$freq)
  }
}

if (!require("processx")) install.packages("processx")
fig <- plot_ly(
  ids=ids_pie,
  labels=labels_pie, 
  parents=parents_pie,
  values=values_pie,
  type='sunburst',
  branchvalues = 'total'
)
orca(fig, "../data/result_genotypted/plots/genotyping.pie.svg")



###################### Genotyping & HardyWeinberg ######################
library(HardyWeinberg)

plot.HWE <- function(dat,pop=NULL,title=NULL,full.legend=F,lab.cex=1){
  require(HardyWeinberg,quietly=T)
  #Gather HW p-values & colors
  HWE.mat <- dat
  HW.p <- HWChisqStats(X=HWE.mat,x.linked=F,pvalues=T)
  HW.cols <- rep("#4DAC26",times=length(HW.p))
  HW.cols[which(HW.p<0.05)] <- "#81F850"
  HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"
  
  #Generate HW plot frame
  par(mar=c(1,1,1,1),bty="n")
  plot(x=1.15*c(-1/sqrt(3),1/sqrt(3)),y=c(-0.15,1.15),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
  segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
           x1=c(0,1/sqrt(3),-1/sqrt(3)),
           y0=c(0,1,0),y1=c(1,0,0))
  HWTernaryPlot(X=HWE.mat,n=max(HWE.mat,na.rm=T),newframe=F,
                vbounds=F,mafbounds=F,
                region=1,vertexlab=NA,
                alpha=0.05,
                curvecols=c("#4DAC26","#81F850",NA,NA),pch=NA)
  
  #Add axes
  text(x=c(-1/sqrt(3),1/sqrt(3)),y=0,labels=c("0/0","1/1"),
       pos=1,cex=lab.cex,xpd=T,font=2)
  text(x=0,y=1,labels="0/1",pos=3,cex=lab.cex,xpd=T,font=2)
  
  #Finish HW plot
  HWTernaryPlot(X=HWE.mat,n=max(HWE.mat,na.rm=T),newframe=F,
                vbounds=F,mafbounds=F,
                region=1,vertexlab=NA,
                alpha=0.03/nrow(HWE.mat),
                curvecols=c("#4DAC26","#AC26A1",NA,NA),
                pch=21,cex=0.3,signifcolour=F,markercol=NA,
                markerbgcol=adjustcolor(HW.cols,alpha=0.25))
  segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
           x1=c(0,1/sqrt(3),-1/sqrt(3)),
           y0=c(0,1,0),y1=c(1,0,0))
  
  #Add legend
  n.pass <- length(which(HW.p>=0.05))
  print(paste("PASS: ",n.pass/length(HW.p),sep=""))
  n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HWE.mat)))
  print(paste("NOMINAL FAILS: ",n.nom/length(HW.p),sep=""))
  n.bonf <- length(which(HW.p<0.05/nrow(HWE.mat)))
  print(paste("BONFERRONI FAILS: ",n.bonf/length(HW.p),sep=""))
  legend("right",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=1,
         legend=c(paste(round(100*(n.pass/nrow(HWE.mat)),0),"%",sep=""),
                  paste(round(100*(n.nom/nrow(HWE.mat)),0),"%",sep=""),
                  paste(round(100*(n.bonf/nrow(HWE.mat)),0),"%",sep="")),
         bty="n",bg=NA,cex=lab.cex)
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=lab.cex,
       labels=paste(title,"\n \n ",sep=""),font=2)
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=lab.cex,
       labels=paste(" \n",prettyNum(max(apply(HWE.mat,1,sum),na.rm=T),big.mark=",")," Samples\n ",sep=""))
  text(x=par("usr")[2],y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),pos=2,cex=lab.cex,
       labels=paste(" \n \n",prettyNum(nrow(HWE.mat),big.mark=",")," SV",sep=""))
}

data_genotypes_normal <- data_genotypes[
  (!is.na(data_sv_details_all$AF_ALL_NGS))&
    data_sv_details_all$MR_NGS<0.05&
    data_sv_details_all$AF_ALL_NGS>0.01&
    data_sv_details_all$AF_ALL_NGS<1,
  colnames(data_genotypes)%in%data_samples_info[data_samples_info$Platform=='NGS',]$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0

data_genotypes_hw <- data_genotypes_normal
data_genotypes_hw$AA <- rowSums(data_genotypes_normal==0)
data_genotypes_hw$AB <- rowSums(data_genotypes_normal==1)
data_genotypes_hw$BB <- rowSums(data_genotypes_normal==2)

data_genotypes_hw <- data_genotypes_hw[,c('AA','AB','BB')]
# png(filename = paste(prefix_filename,'genotyping.hardyweinberg.png',sep = '.'), width = 1000, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.hardyweinberg','pdf',sep = '.'), width = 8, height = 6)
plot.HWE(data_genotypes_hw, lab.cex = 1)
dev.off()


data_genotypes_normal <- data_genotypes[
  (!is.na(data_sv_details_all$AF_ALL_NGS))&
    data_sv_details_all$MR_NGS<0.05&
    data_sv_details_all$AF_ALL_NGS>0.01&
    data_sv_details_all$AF_ALL_NGS<1&
  data_sv_details_all$SVTYPE=='INS',
  colnames(data_genotypes)%in%data_samples_info[data_samples_info$Platform=='NGS',]$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
data_genotypes_hw <- data_genotypes_normal
data_genotypes_hw$AA <- rowSums(data_genotypes_normal==0)
data_genotypes_hw$AB <- rowSums(data_genotypes_normal==1)
data_genotypes_hw$BB <- rowSums(data_genotypes_normal==2)

data_genotypes_hw <- data_genotypes_hw[,c('AA','AB','BB')]
# png(filename = paste(prefix_filename,'genotyping.hardyweinberg.ins.png',sep = '.'), width = 1000, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.hardyweinberg.ins','pdf',sep = '.'), width = 8, height = 6)
plot.HWE(data_genotypes_hw, lab.cex = 1)
dev.off()

data_genotypes_normal <- data_genotypes[
  (!is.na(data_sv_details_all$AF_ALL_NGS))&
    data_sv_details_all$MR_NGS<0.05&
    data_sv_details_all$AF_ALL_NGS>0.01&
    data_sv_details_all$AF_ALL_NGS<1&
    data_sv_details_all$SVTYPE=='DEL',
  colnames(data_genotypes)%in%data_samples_info[data_samples_info$Platform=='NGS',]$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
data_genotypes_hw <- data_genotypes_normal
data_genotypes_hw$AA <- rowSums(data_genotypes_normal==0)
data_genotypes_hw$AB <- rowSums(data_genotypes_normal==1)
data_genotypes_hw$BB <- rowSums(data_genotypes_normal==2)

data_genotypes_hw <- data_genotypes_hw[,c('AA','AB','BB')]
# png(filename = paste(prefix_filename,'genotyping.hardyweinberg.del.png',sep = '.'), width = 1000, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.hardyweinberg.del','pdf',sep = '.'), width = 8, height = 6)
plot.HWE(data_genotypes_hw, lab.cex = 1)
dev.off()

data_genotypes_normal <- data_genotypes[
  (!is.na(data_sv_details_all$AF_ALL_NGS))&
    data_sv_details_all$MR_NGS<0.05&
    data_sv_details_all$AF_ALL_NGS>0.01&
    data_sv_details_all$AF_ALL_NGS<1&
    data_sv_details_all$SVTYPE=='INV',
  colnames(data_genotypes)%in%data_samples_info[data_samples_info$Platform=='NGS',]$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
data_genotypes_hw <- data_genotypes_normal
data_genotypes_hw$AA <- rowSums(data_genotypes_normal==0)
data_genotypes_hw$AB <- rowSums(data_genotypes_normal==1)
data_genotypes_hw$BB <- rowSums(data_genotypes_normal==2)

data_genotypes_hw <- data_genotypes_hw[,c('AA','AB','BB')]
# png(filename = paste(prefix_filename,'genotyping.hardyweinberg.inv.png',sep = '.'), width = 1000, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.hardyweinberg.inv','pdf',sep = '.'), width = 8, height = 6)
plot.HWE(data_genotypes_hw, lab.cex = 1)
dev.off()

data_genotypes_normal <- data_genotypes[
  (!is.na(data_sv_details_all$AF_ALL_NGS))&
    data_sv_details_all$MR_NGS<0.05&
    data_sv_details_all$AF_ALL_NGS>0.01&
    data_sv_details_all$AF_ALL_NGS<1&
    data_sv_details_all$SVTYPE=='DUP',
  colnames(data_genotypes)%in%data_samples_info[data_samples_info$Platform=='NGS',]$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
data_genotypes_hw <- data_genotypes_normal
data_genotypes_hw$AA <- rowSums(data_genotypes_normal==0)
data_genotypes_hw$AB <- rowSums(data_genotypes_normal==1)
data_genotypes_hw$BB <- rowSums(data_genotypes_normal==2)

data_genotypes_hw <- data_genotypes_hw[,c('AA','AB','BB')]
pdf(file = paste(prefix_filename,'genotyping.hardyweinberg.dup','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'genotyping.hardyweinberg.dup.png',sep = '.'), width = 1000, height = 1000, res = val_res)
plot.HWE(data_genotypes_hw, lab.cex = 1)
dev.off()

###################### TODO Genotyping & AF correlation ##########################
data_sv_database <- data_sv_details_all[data_sv_details_all$AF_Database>0.01&data_sv_details_all$AF_ALL_NGS>0.01&data_sv_details_all$MR_NGS<0.05,]
ggplot(data_sv_database, aes(x=AF_ALL_NGS*100, y=AF_Database*100)) + geom_point() +
  scale_x_continuous(trans = 'log10', 
                     breaks = c(1,10,100)) + 
  scale_y_continuous(trans = 'log10', 
                     breaks = c(1,10,100))

data_sv_database <- data_sv_details_all[data_sv_details_all$AF_ALL>0.01&data_sv_details_all$MR<0.05&data_sv_details_all$AF_ALL_NGS>0.01&data_sv_details_all$MR_NGS<0.05,]
ggplot(data_sv_database, aes(x=AF_ALL_NGS*100, y=AF_ALL*100)) + geom_point() +
  scale_x_continuous(trans = 'log10', 
                     breaks = c(1,10,100)) + 
  scale_y_continuous(trans = 'log10', 
                     breaks = c(1,10,100))

###################### Genotyping & GROUP_SUPP ######################
tmp <- subset(data_plot, SUPP_NGS>0)
df.new<-ddply(tmp,.(SVTYPE),plyr::summarise,
              prop=prop.table(table(GROUP_SUPP_NGS)),
              SUPP=names(table(GROUP_SUPP_NGS)))
df.new$SUPP <- factor(df.new$SUPP, levels = c('Singleton', 'Polymorphic', 'Major', 'Shared'))
pdf(file = paste(prefix_filename,'genotyping.group_supp.bar','pdf',sep = '.'), width = 3, height = 3)
ggplot(df.new, aes(SUPP, prop, fill=SVTYPE)) + 
  geom_bar(stat="identity",position = 'dodge') + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.9,0.8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_svtype) + 
  xlab("") + 
  ylab("Prop")
dev.off()

###################### Genotyping & SUPP ######################
pdf(file = paste(prefix_filename,'genotyping.supp.bar','pdf',sep = '.'), width = 6, height = 3)
ggplot(data=data_tmp, aes(x=SUPP, fill = GROUP_SUPP_NGS)) +
  geom_bar() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.5,0.9),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold", angle = 45),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold")) +
  scale_fill_manual(values=config_color_group_supp_han) + 
  xlab('Samples Support') + 
  ylab("Discovery")
dev.off()

###################### Genotyping & NO SUPPORT ######################
data_tmp <- data_plot
data_tmp$GROUP_SUPP_NGS <- as.character(data_tmp$GROUP_SUPP_NGS)
data_tmp[is.na(data_tmp$GROUP_SUPP_NGS),]$GROUP_SUPP_NGS <- 'No Support'
data_tmp$GROUP_SUPP_NGS <- factor(data_tmp$GROUP_SUPP_NGS, levels = c('No Support','Singleton', 'Polymorphic', 'Major', 'Shared'))

data_tmp_4 <- subset(data_tmp, GROUP_SUPP_NGS=='No Support'& SVTYPE!='TRA')
data_tmp_4$Group_missing <- 'Others'
data_tmp_4[data_tmp_4$MR_NGS == 0, ]$Group_missing <- 'MR=0'
data_tmp_4[data_tmp_4$MR_NGS == 1, ]$Group_missing <- 'MR=1'
data_tmp_4$Group_missing <- factor(data_tmp_4$Group_missing, levels = c('MR=0','MR=1','Others'))


data_tmp_2 <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & MR_NGS==0)
# png(filename = paste(prefix_filename,'genotyping.nosupport.mr0.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.nosupport.mr0','pdf',sep = '.'), width = 3, height = 3)
ggplot(data_tmp_2, aes(GROUP_SUPP, fill=GROUP_SUPP)) + 
  geom_bar(position = 'dodge') + 
  # coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.9,0.9),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("Discovery") + 
  geom_text(aes(label = ..count..), stat= "count", size=2, vjust=1.6, color='white', fontface="bold")
dev.off()

data_tmp_3 <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & MR_NGS==1)
# png(filename = paste(prefix_filename,'genotyping.nosupport.mr1.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.nosupport.mr1','pdf',sep = '.'), width = 3, height = 3)
ggplot(data_tmp_3, aes(GROUP_SUPP, fill=GROUP_SUPP)) + 
  geom_bar(position = 'dodge') + 
  # coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.85,0.8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-15,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp) + 
  xlab("") + 
  ylab("Discovery") + 
  geom_text(aes(label = ..count..), stat= "count", size=2, vjust=1.6, color='white', fontface="bold")
dev.off()

###################### Genotyping & NO SUPPORT & MR=0 Major & Shared ######################

data_sv_group <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & 
                          MR_NGS <= 0.05 & 
                          GROUP_SUPP%in%c('Shared', 'Major'))
data_sv_group <- subset(data_tmp, GROUP_SUPP_NGS=='No Support' & MR_NGS <= 0.05)
# data_sv_group <- data_tmp

data_sv_group$GROUP_REP <- 'None'
data_sv_group[data_sv_group$RepClass!='None',]$GROUP_REP <- 'Repclass'
data_sv_group[data_sv_group$RepClass=='None'&data_sv_group$SD!='',]$GROUP_REP <- 'SD'
data_sv_group[data_sv_group$GROUP_REP=='None'&
                (data_sv_group$Repeats_type_left!='None'|data_sv_group$Repeats_type_right!='None'),]$GROUP_REP <- 'flanked'

data_tmp_1 <- plyr::count(data_sv_group,'GROUP_REP')
pdf(file = paste(prefix_filename,'genotyping.nosupport.mr0.repclass.pie','pdf',sep = '.'), 
    width = 3, height = 3)
ggplot(data_tmp_1, aes(x="",y=freq,fill=GROUP_REP)) + 
  geom_bar(width=1,stat="identity") + 
  coord_polar("y",start=0) + 
  scale_fill_manual(values=c('None'='gray',
                             'Repclass'='dark red', 
                             'SD'='dark green',
                             'flanked'='dark blue')) + 
  geom_text(aes(y=freq/5+c(0,cumsum(freq)[-length(freq)])), 
            label=percent(data_tmp_1$freq/sum(data_tmp_1$freq)),
            size=2,
            color='white', 
            fontface="bold") + 
  theme_minimal() + 
  theme(legend.title=element_blank(),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

###################### Genotyping & SVTYPE & LD ######################
data_sv_group <- data_plot[data_plot$LD>0&data_plot$SVTYPE!='TRA'&!is.na(data_plot$LD),]
# png(filename = paste(prefix_filename,'genotyping.svtype.LD.boxplot.png',sep = '.'), width = 1200, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.svtype.LD.boxplot','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_group, aes(x=SVTYPE, y=LD,group=SVTYPE)) + 
  geom_boxplot(aes(fill=SVTYPE),outlier.size = 0.1) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        # panel.border = element_blank(),
        # panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_svtype) + 
  xlab('') + 
  ylab("max LD")
dev.off()

###################### Genotyping & Enrichment EQTL&GWAS ######################
data_tmp <- plyr::count(data_snps_eqtl,c('SVID', 'Gene'))
data_tmp <- data_tmp[data_tmp$Gene!='.',]
data_sv_group <- plyr::count(data_sv_details_all[data_sv_details_all$SVID%in%unique(data_tmp$SVID)&data_sv_details_all$AF_ALL_NGS>0.05&!is.na(data_sv_details_all$AF_ALL_NGS),],
                             c('GROUP_SUPP_NGS','SVTYPE'))
# png(filename = paste(prefix_filename,'genotyping.svtype.LD.eqtl.png',sep = '.'), width = 1200, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.svtype.LD.eqtl','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_group, aes(x=GROUP_SUPP_NGS, y=freq, fill=SVTYPE)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(legend.title=element_blank(),
        legend.position='bottom',
        # panel.border = element_blank(),
        # panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_svtype) + 
  xlab('') + 
  ylab("")
dev.off()

data_tmp <- plyr::count(data_snps_gwas,c('SVID', 'SNPID'))
data_tmp <- data_tmp[data_tmp$SNPID!='.',]
data_sv_group <- plyr::count(data_sv_details_all[data_sv_details_all$SVID%in%unique(data_tmp$SVID)&data_sv_details_all$AF_ALL_NGS>0.05&!is.na(data_sv_details_all$AF_ALL_NGS),],
                             c('GROUP_SUPP_NGS','SVTYPE'))
# png(filename = paste(prefix_filename,'genotyping.svtype.LD.gwas.png',sep = '.'), width = 1200, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.svtype.LD.gwas','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_group, aes(x=GROUP_SUPP_NGS, y=freq, fill=SVTYPE)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(legend.title=element_blank(),
        legend.position='bottom',
        # panel.border = element_blank(),
        # panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_svtype) + 
  xlab('') + 
  ylab("")
dev.off()

###################### target-90 EQTL ######################

# png(filename = paste(prefix_filename,'target90.eqtl.png',sep = '.'), width = 1200, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'target90.eqtl','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_pairs_target_90, aes(x=GROUP, y=slope, fill=GROUP)) + 
  geom_violin() + 
  theme_minimal() + 
  theme(legend.position = 'none') +
  xlab('') +
  ylab('effect size')
dev.off()

###################### target-90 enrichment ######################
data.gene.background <- read.xlsx('../data/result_genotypted/target-90/deng.background.xlsx')
data.gene.target90 <- read.table('../data/result_genotypted/target-90/Tables-qc.genelist.txt', header = FALSE)
colnames(data.gene.target90) <- 'Gene'
data.gene.background$Genes <- NA
data.deng.prior <- read.table('../data/result_genotypted/target-90/deng.prior.supplement.txt', header = FALSE)
colnames(data.deng.prior) <- 'Gene'

data.deng.prior.new <- data.deng.prior$Gene[!str_detect(data.gene.background[data.gene.background$Gene.Set=='prior literature',]$Gene.Symbol, as.character(data.deng.prior$Gene))]
data.gene.background$N1 <- data.gene.background$N1 + length(data.deng.prior.new)
data.gene.background$N0 <- data.gene.background$N0 - length(data.deng.prior.new)
data.gene.background$Number.of.Genes <-  data.gene.background$Number.of.Genes + length(data.deng.prior.new)
data.gene.background[data.gene.background$Gene.Set=='prior literature',]$Gene.Symbol <- paste(
  data.gene.background[data.gene.background$Gene.Set=='prior literature',]$Gene.Symbol,
  paste(data.deng.prior.new,collapse = ','),
  sep = ',')

for (index_row in seq(nrow(data.gene.background))) {
  
  data.gene.background[index_row,]$Genes <- paste(data.gene.target90$Gene[str_detect(data.gene.background[index_row,]$Gene.Symbol,
                                                       as.character(data.gene.target90$Gene))],collapse = ',')
  
  data.gene.background[index_row,]$A1 <- sum(str_detect(data.gene.background[index_row,]$Gene.Symbol,
                                                        as.character(data.gene.target90$Gene)))
  data.gene.background[index_row,]$A0 <- sum(!str_detect(data.gene.background[index_row,]$Gene.Symbol,
                                                        as.character(data.gene.target90$Gene)))
  
  tmp <- fisher.test(matrix(c(data.gene.background[index_row,]$A1,
                              data.gene.background[index_row,]$A0,
                              data.gene.background[index_row,]$N1,
                              data.gene.background[index_row,]$N0),nrow = 2))
  data.gene.background[index_row,]$Odds.Ratio <- tmp$estimate
  data.gene.background[index_row,]$p.value <- tmp$p.value
}

write.table(data.gene.background,
            file = '../data/result_genotypted/target-90/enrich.update.txt',
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

data.gene.background <- read.xlsx('../data/result_genotypted/target-90/enrich.xlsx')

for (index_row in seq(nrow(data.gene.background))) {

  tmp <- fisher.test(matrix(c(data.gene.background[index_row,]$A1,
                              data.gene.background[index_row,]$A0,
                              data.gene.background[index_row,]$N1,
                              data.gene.background[index_row,]$N0),nrow = 2))
  data.gene.background[index_row,]$Odds.Ratio <- tmp$estimate
  data.gene.background[index_row,]$p.value <- tmp$p.value
}

data.gene.background$Gene.Set <- factor(data.gene.background$Gene.Set, levels = rev(data.gene.background$Gene.Set[order(data.gene.background$Odds.Ratio,decreasing = TRUE)]))
pdf(file = paste(prefix_filename,'enrich.pathway','pdf',sep = '.'), width = 3.5, height = 2)
ggplot(data=data.gene.background[2:11,],aes(x=Gene.Set,
                                          y=A1,
                                          fill=p.value)) + 
  geom_bar(stat="identity") + coord_flip() + 
  theme(panel.background=element_rect(fill='transparent',colour = 'black'),
        legend.key.size=unit(0.2, 'cm'),
        legend.title=element_blank(),
        legend.text = element_text(size = 2.5,face = "bold"),
        axis.text.x = element_text(size = 2.5,face = "bold"),
        axis.text.y = element_text(color="black",size=2.5,face = "bold"),
        axis.title.y = element_text(size = 2.5,face = "bold")) + 
  scale_fill_gradient(low="red",high="blue", breaks=c(0.05, 0.1), limits=c(0, 0.15)) +
  xlab("") + 
  ylab("")
dev.off()

data.gene.background <- read.csv('../data/result_genotypted/target-90/metascape.target.90/metascape_result.tsv', sep = '\t')
data.gene.background$count <- as.integer(str_split_fixed(data.gene.background$InTerm_InList,pattern = '/', n=2)[,1])
data.gene.background$p <- 10 ** data.gene.background$LogP
data.gene.background$Description <- factor(data.gene.background$Description, levels = rev(data.gene.background$Description))
pdf(file = paste(prefix_filename,'enrich.go','pdf',sep = '.'), width = 3.5, height = 2)
ggplot(data=data.gene.background[1:10,],aes(x=Description,
                                            y=count,
                                            fill=p)) + 
  geom_bar(stat="identity") + coord_flip() + 
  theme(panel.background=element_rect(fill='transparent',colour = 'black'),
        legend.key.size=unit(0.2, 'cm'),
        legend.title=element_blank(),
        legend.text = element_text(size = 2.5,face = "bold"),
        axis.text.x = element_text(size = 2.5,face = "bold"),
        axis.text.y = element_text(color="black",size=2.5,face = "bold"),
        axis.title.y = element_text(size = 2.5,face = "bold")) + 
  scale_fill_gradient(low="red",high="blue") +
  xlab("") + 
  ylab("")
dev.off()

###################### Genotyping & Gene ####################################################
data_sv_group <- plyr::count(data_plot[data_plot$SVLEN<1000000&data_plot$SUPP_NGS>0&!is.na(data_plot$SUPP_NGS),],c('GROUP_SV','GROUP_SUPP_NGS'))
data_sv_group$label <- 0
data_sv_group$GROUP_SV <- factor(data_sv_group$GROUP_SV)

for (type_sv in levels(data_sv_group$GROUP_SV)){
  data_sv_group[data_sv_group$GROUP_SV==type_sv,]$label <- sum(data_sv_group[data_sv_group$GROUP_SV==type_sv,]$freq)
}

# tiff(filename = paste(prefix_filename,'sv_summary.tiff',sep = '.'), width = 1200, height = 1000, res = 300)
# png(filename = paste(prefix_filename,'genotyping.gene.group_supp.bar.png',sep = '.'), width = 1200, height = 1000, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.gene.group_supp.bar','pdf',sep = '.'), width = 8, height = 6)
ggplot(data_sv_group, aes(x = GROUP_SV, y = freq, fill = GROUP_SUPP_NGS)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), expand = expand_scale(mult = .1)) + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(-10,0,0,0))) +
  scale_fill_manual(values=config_color_group_supp_han) + 
  xlab("") + 
  ylab("") +
  geom_text(aes(label = label, y= ..prop..), stat= "count", hjust = -0.1, size=2.5)
dev.off()

###################### Genotyping & POP VAF ######################

data_tmp <- data_plot[!is.na(data_plot$AF_ALL_NGS)&data_plot$MR_NGS<0.05&data_plot$AF_ALL_NGS>0&data_plot$AF_ALL_NGS<0.1,]
pdf(file = paste(prefix_filename,'genotyping.pop.vaf.group_pop','pdf',sep = '.'), width = 3, height = 3)
# png(filename = paste(prefix_filename,'genotyping.pop.vaf.group_pop.png',sep = '.'), width = 1200, height = 1200, res = val_res)
ggplot(data=data_tmp, aes(AF_ALL_NGS, stat(count), fill = GROUP_POP_NGS)) +
  geom_density(position = "fill",bw=0.005) + 
  scale_fill_manual(values=config_color_group_pop)+ 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('VAF') + 
  ylab("Proportion")
dev.off()

data_tmp <- data_plot[!is.na(data_plot$AF_ALL_NGS)&data_plot$MR_NGS<0.05&data_plot$AF_ALL_NGS>0&data_plot$AF_ALL_NGS<0.1,]
# png(filename = paste(prefix_filename,'genotyping.pop.vaf.group_pop_detail.png',sep = '.'), width = 1200, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.pop.vaf.group_pop_detail','pdf',sep = '.'), width = 3, height = 3)
ggplot(data=data_tmp, aes(AF_ALL_NGS, stat(count), fill = GROUP_POP_DETAIL_NGS)) +
  geom_density(position = "fill",bw=0.005) + 
  scale_fill_manual(values=config_color_group_pop)+ 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('VAF') + 
  ylab("Proportion")
dev.off()

###################### Genotyping & POP PCA ######################

group_sample <- as.data.frame(data_samples_info[data_samples_info$Source=='T15H10',c('SampleID','Group','Group_Detail','Location')])
group_sample <- group_sample[order(group_sample$SampleID),]
rownames(group_sample) <- group_sample$SampleID
# group_sample <- as.data.frame(data_samples_info[1:25, c('SampleID','Group', 'Location')])
# rownames(group_sample) <- group_sample$SampleID

data_genotypes_normal <- data_genotypes[data_sv_details_all$MR<0.05,colnames(data_genotypes)%in%group_sample$SampleID]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
# tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/(1+(2*sqrt(data_sv_details_all$AF_ALL*(1-data_sv_details_all$AF_ALL))))
p <- (1 + rowSums(data_genotypes_normal))/(2+2*nrow(group_sample))
tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/2*sqrt(p*(1-p))
data_sv_details_all_matrix.pca <- prcomp(t(as.matrix(tmp)))
# png(filename = paste(prefix_filename,'pop.pca.ont.png',sep = '.'), width = 1200, height = 1000, res = 300)
pdf(file = paste(prefix_filename,'pop.pca.ont','pdf',sep = '.'), width = 8, height = 6)
autoplot(data_sv_details_all_matrix.pca, data = group_sample, colour='Group_Detail', label=TRUE, frame=TRUE)
dev.off()


group_sample <- as.data.frame(data_samples_info[!data_samples_info$Source%in%c('T15H10','HUM'),c('SampleID','Group','Group_Detail','Location')])
group_sample <- group_sample[order(group_sample$SampleID),]
rownames(group_sample) <- group_sample$SampleID

data_genotypes_normal <- data_genotypes[!is.na(data_sv_details_all$AF_ALL_NGS)&data_sv_details_all$MR_NGS<0.05&data_sv_details_all$AF_ALL_NGS>0.01,colnames(data_genotypes)%in%group_sample$SampleID]
data_genotypes_normal <- data_genotypes_normal[, order(colnames(data_genotypes_normal))]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
# tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/(1+(2*sqrt(data_sv_details_all$AF_ALL_NGS*(1-data_sv_details_all$AF_ALL_NGS))))
p <- (1 + rowSums(data_genotypes_normal))/(2+2*nrow(group_sample))
tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/2*sqrt(p*(1-p))
data_sv_details_all_matrix.pca <- prcomp(t(as.matrix(tmp)))
PCi<-data.frame(data_sv_details_all_matrix.pca$x,Group_Detail=group_sample$Group_Detail)

# png(filename = paste(prefix_filename,'genotyping.pop.pca.png',sep = '.'), width = 1200, height = 1000, res = 300)
pdf(file = paste(prefix_filename,'genotyping.pop.pca','pdf',sep = '.'), width = 3, height = 3)
p1 <- ggplot(PCi,aes(x=PC1,y=PC2,color=Group_Detail))+
  geom_point(size=0.8)+
  scale_color_manual(values = config_color_group_pop) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.1,0.9),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        # panel.background = element_blank(),
        # # panel.border =  element_rect(size=0.6, colour = "black"),
        # axis.line = element_line(size=0.6, colour = "black"),
        # axis.line.x.top = element_line(size=0.6, colour = "black"),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))
        )
p1
dev.off()

group_sample <- as.data.frame(data_samples_info[data_samples_info$Group_Detail!='AFR'&!data_samples_info$Source%in%c('T15H10','HUM'),c('SampleID','Group','Group_Detail','Location')])
group_sample <- group_sample[order(group_sample$SampleID),]
rownames(group_sample) <- group_sample$SampleID

data_genotypes_normal <- data_genotypes[!is.na(data_sv_details_all$AF_ALL_NGS)&data_sv_details_all$MR_NGS<0.05&data_sv_details_all$AF_ALL_NGS>0.01,colnames(data_genotypes)%in%group_sample$SampleID]
data_genotypes_normal <- data_genotypes_normal[, order(colnames(data_genotypes_normal))]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
# tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/(1+(2*sqrt(data_sv_details_all$AF_ALL_NGS*(1-data_sv_details_all$AF_ALL_NGS))))
p <- (1 + rowSums(data_genotypes_normal))/(2+2*nrow(group_sample))
tmp <- (data_genotypes_normal- rowMeans(data_genotypes_normal))/2*sqrt(p*(1-p))
data_sv_details_all_matrix.pca <- prcomp(t(as.matrix(tmp)))
PCi<-data.frame(data_sv_details_all_matrix.pca$x,Group_Detail=group_sample$Group_Detail)

pdf(file = paste(prefix_filename,'genotyping.pop.pca.noAFR','pdf',sep = '.'), width = 3, height = 3)
# png(filename = paste(prefix_filename,'genotyping.pop.pca.noAFR.png',sep = '.'), width = 1200, height = 1000, res = 300)
p2 <- ggplot(PCi,aes(x=PC1,y=PC2,color=Group_Detail))+
  geom_point(size=0.8)+ #Size and alpha just for fun
  scale_color_manual(values = config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position=c(0.1,0.9),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.3, 'cm'),
        # panel.background = element_blank(),
        # panel.border =  element_rect(size=0.6, colour = "black"),
        # axis.line = element_line(size=0.6, colour = "black"),
        # axis.line.x.top = element_line(size=0.6, colour = "black"),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold"),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0)))
p2
dev.off()

###################### SVTYPE & Watterson estimator ##########################

num_harmonic <- sum(1/c(1:22))
func_calc_harmonic <- function(array_chrom){
  
  res_harmonic <- 1:length(array_chrom)
  index_res <- 1
  for (num_chrom in array_chrom){
    res_harmonic[index_res] <- sum(1/c(1:num_chrom))
    index_res <- index_res + 1
  }
  return(res_harmonic)
}


estimator_watterson_hans <- nrow(data_sv_details_all[data_sv_details_all$SUPP_HANS_NGS>0,])/num_harmonic
estimator_watterson_hann <- nrow(data_sv_details_all[data_sv_details_all$SUPP_HANN_NGS>0,])/num_harmonic
estimator_watterson_tibg <- nrow(data_sv_details_all[data_sv_details_all$SUPP_TIBG_NGS>0,])/num_harmonic
estimator_watterson_tibl <- nrow(data_sv_details_all[data_sv_details_all$SUPP_TIBL_NGS>0,])/num_harmonic

rate_mutation <- c(estimator_watterson_hans, estimator_watterson_hann ,estimator_watterson_tibg, estimator_watterson_tibl)/(4*10000)

data_sv_group_hans <- plyr::count(data_sv_details_all[data_sv_details_all$SUPP_HANS_NGS>0&data_sv_details_all$SVTYPE!='TRA',],c('SVTYPE','Chrom'))
data_sv_group_hann <- plyr::count(data_sv_details_all[data_sv_details_all$SUPP_HANN_NGS>0&data_sv_details_all$SVTYPE!='TRA',],c('SVTYPE','Chrom'))
data_sv_group_tibg <- plyr::count(data_sv_details_all[data_sv_details_all$SUPP_TIBG_NGS>0&data_sv_details_all$SVTYPE!='TRA',],c('SVTYPE','Chrom'))
data_sv_group_tibl <- plyr::count(data_sv_details_all[data_sv_details_all$SUPP_TIBL_NGS>0&data_sv_details_all$SVTYPE!='TRA',],c('SVTYPE','Chrom'))

data_sv_group <- rbind(plyr::count(data_sv_group_hans, 'SVTYPE'), plyr::count(data_sv_group_hann, 'SVTYPE'),plyr::count(data_sv_group_tibg, 'SVTYPE'),plyr::count(data_sv_group_tibl, 'SVTYPE'))
data_sv_group$Group <- c(rep('HANS',4),rep('HANN',4),rep('TIBG',4),rep('TIBL',4))
data_sv_group$Chrom <- rbind(plyr::count(data_sv_group_hans[,-3],c('SVTYPE')),
                             plyr::count(data_sv_group_hann[,-3],c('SVTYPE')),
                             plyr::count(data_sv_group_tibg[,-3],c('SVTYPE')),
                             plyr::count(data_sv_group_tibl[,-3],c('SVTYPE')))$freq



data_sv_group$Rate <- data_sv_group$freq/(4*10000*func_calc_harmonic(data_sv_group$Chrom))


pdf(file = paste(prefix_filename,'svtype.mutation.ngs','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'svtype.mutation.ngs.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data=data_sv_group, aes(x=SVTYPE, y=Rate, group=Group)) +
  geom_line(aes(color=Group)) + 
  geom_point(aes(color=Group)) + 
  scale_color_manual(values=config_color_group_pop) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("Mutation Rate")
dev.off()


estimator_watterson <- nrow(data_sv_details_all[data_sv_details_all$SUPP_NGS>0&!is.na(data_sv_details_all$SUPP_NGS),])/num_harmonic
data_sv_group_mean <- plyr::count(data_sv_details_all[data_sv_details_all$SUPP_NGS>0&!is.na(data_sv_details_all$SUPP_NGS),], c('SVTYPE','Chrom'))

data_sv_group <- plyr::count(data_sv_group_mean, 'SVTYPE')
data_sv_group$Chrom <- plyr::count(data_sv_group_mean[,-3],c('SVTYPE'))$freq

data_sv_group$Rate <- data_sv_group$freq/(4*10000*func_calc_harmonic(data_sv_group$Chrom))
#
pdf(file = paste(prefix_filename,'svtype.mutation.mean','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'svtype.mutation.mean.png',sep = '.'), width = 1200, height = 1000, res = val_res)
ggplot(data=data_sv_group, aes(x=SVTYPE, y=Rate)) +
  geom_line() +
  geom_point() +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) +
  xlab('') +
  ylab("Mutation Rate")
dev.off()


###################### Genotyping & POP PHt ######################

group_sample <- as.data.frame(data_samples_info[!data_samples_info$Source%in%c('T15H10'),
                                                c('SampleID','Group','Group_Detail','Location')])
group_sample <- group_sample[order(group_sample$SampleID),]
group_sample$PHt <- 0
group_sample$PHot <- 0
rownames(group_sample) <- group_sample$SampleID

data_genotypes_normal <- data_genotypes[!is.na(data_sv_details_all$AF_ALL_NGS)&data_sv_details_all$MR_NGS<0.05&data_sv_details_all$AF_ALL_NGS>0.01,colnames(data_genotypes)%in%group_sample$SampleID]
data_genotypes_normal <- data_genotypes_normal[, order(colnames(data_genotypes_normal))]
data_genotypes_normal[data_genotypes_normal==-1] <- 0
data_genotypes_normal_type <- data_sv_details_all[!is.na(data_sv_details_all$AF_ALL_NGS)&data_sv_details_all$MR_NGS<0.05&data_sv_details_all$AF_ALL_NGS>0.01,]$SVTYPE


for (tmp_sample in group_sample$SampleID) {
  group_sample[group_sample$SampleID==tmp_sample,]$PHt <- sum(data_genotypes_normal[data_genotypes_normal_type%in%c('INS','DEL'),tmp_sample]==1)
  group_sample[group_sample$SampleID==tmp_sample,]$PHot <- sum(data_genotypes_normal[data_genotypes_normal_type%in%c('INS','DEL'),tmp_sample]==2)
}

group_sample_svtype <- rbind(group_sample,group_sample)
group_sample_svtype$SVTYPE <- c(rep('DEL', nrow(group_sample)),
                                rep('INS', nrow(group_sample)))

for (tmp_svtype in unique(group_sample_svtype$SVTYPE)) {
  for (tmp_sample in group_sample$SampleID) {
    group_sample_svtype[group_sample_svtype$SampleID==tmp_sample&group_sample_svtype$SVTYPE==tmp_svtype,]$PHt <- sum(data_genotypes_normal[data_genotypes_normal_type==tmp_svtype,tmp_sample]==1)
    group_sample_svtype[group_sample_svtype$SampleID==tmp_sample&group_sample_svtype$SVTYPE==tmp_svtype,]$PHot <- sum(data_genotypes_normal[data_genotypes_normal_type==tmp_svtype,tmp_sample]==2)
  }
}

# group_sample <-group_sample[!group_sample$Group_Detail%in%c('HUM','AFR'),]
group_sample <-group_sample[!group_sample$Group_Detail%in%c('HUM'),]
group_sample_svtype <-group_sample_svtype[!group_sample_svtype$Group_Detail%in%c('HUM'),]
# group_sample <- group_sample[!group_sample$SampleID%in%c('WGC025266D', 'WGC025273D'),]
group_sample$Location <- factor(group_sample$Location, levels = unique(c(
  group_sample[group_sample$Group_Detail=='TIBG',]$Location,
  group_sample[group_sample$Group_Detail=='TIBL',]$Location,
  group_sample[group_sample$Group_Detail=='HANN',]$Location,
  group_sample[group_sample$Group_Detail=='HANS',]$Location,
  group_sample[group_sample$Group_Detail=='AFR',]$Location
)))

pairwise.t.test(group_sample$PHot, g = group_sample$Group_Detail, p.adjust.method = 'bonferroni')
pairwise.t.test(group_sample$PHt, g = group_sample$Group_Detail, p.adjust.method = 'bonferroni')
pdf(file = paste(prefix_filename,'genotyping.pop.het.loc','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'genotyping.pop.het.loc.png',sep = '.'), width = 1200, height = 1000, res = 300)
p1 <- ggplot(group_sample, aes(x=Location, y=PHt, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.line.y = element_line(size=0.6, colour = "black"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#het")
p2 <- ggplot(group_sample, aes(x=Location, y=PHot, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#hom")
# ggarrange(p1, p2,ncol=2,nrow = 1)
multiplot(p1,p2,cols = 2)
dev.off()


pdf(file = paste(prefix_filename,'genotyping.pop.hetANDhom','pdf',sep = '.'), width = 6, height = 3)
# png(filename = paste(prefix_filename,'genotyping.pop.hetANDhom.png',sep = '.'), width = 2400, height = 1000, res = 300)
p1 <- ggplot(group_sample, aes(x=Group_Detail, y=PHot, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() +
  # geom_signif(comparisons = list(c('TIBG','TIBL')), map_signif_level=TRUE)
  #geom_signif(comparisons = list(levels(factor(group_sample$Group_Detail))), map_signif_level=TRUE)
  scale_fill_manual(values=config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 6,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.line.y = element_line(size=0.6, colour = "black"),
        axis.text.x = element_text(size = 6,face = "bold",vjust = 0.6),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#hom")
p2 <- ggplot(group_sample, aes(x=Group_Detail, y=PHt, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  scale_fill_manual(values=config_color_group_pop) + 
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,face = "bold",vjust = 0.6),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#het")
# ggarrange(p1, p2,ncol=2,nrow = 1)
multiplot(p1,p2,cols = 2)
dev.off()


pdf(file = paste(prefix_filename,'genotyping.pop.hetANDhom.ins','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'genotyping.pop.hetANDhom.ins.png',sep = '.'), width = 2400, height = 1000, res = 300)
p1 <- ggplot(group_sample_svtype[group_sample_svtype$SVTYPE=='INS',], aes(x=Group_Detail, y=PHot, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  scale_fill_manual(values=config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.line.y = element_line(size=0.6, colour = "black"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#hom")
p2 <- ggplot(group_sample_svtype[group_sample_svtype$SVTYPE=='INS',], aes(x=Group_Detail, y=PHt, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  scale_fill_manual(values=config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#het")
multiplot(p1,p2,cols = 2)
dev.off()


pdf(file = paste(prefix_filename,'genotyping.pop.hetANDhom.ins','pdf',sep = '.'), width = 8, height = 6)
# png(filename = paste(prefix_filename,'genotyping.pop.hetANDhom.del.png',sep = '.'), width = 2400, height = 1000, res = 300)
p1 <- ggplot(group_sample_svtype[group_sample_svtype$SVTYPE=='DEL',], aes(x=Group_Detail, y=PHot, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  scale_fill_manual(values=config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        legend.text = element_text(size = 10,face = "bold"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.line.y = element_line(size=0.6, colour = "black"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#hom")
p2 <- ggplot(group_sample_svtype[group_sample_svtype$SVTYPE=='DEL',], aes(x=Group_Detail, y=PHt, fill=Group_Detail)) + 
  geom_violin() + 
  coord_flip() + 
  scale_fill_manual(values=config_color_group_pop) +
  theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size=unit(0.5, 'cm'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.6, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 10,face = "bold",vjust = 0.6),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.x.bottom = element_text(margin = margin(5,0,0,0))) + 
  xlab('') + 
  ylab("#het")
multiplot(p1,p2,cols = 2)
dev.off()

###################### Genotyping & POP FST ######################
library(qqman)
library(Cairo)

  Fstfile <- read.table(paste(path_data, 'fst/merge.paragraph.genotypes.TIBvsHAN.20k_5k.windowed.weir.fst', sep = ''), header = T, stringsAsFactors = F)
  SNP <- paste(Fstfile[,1], Fstfile[,2], sep = ':')
  Fstfile <- cbind(SNP, Fstfile)
  colnames(Fstfile) <- c('SNP', 'CHR', 'POS','END','Bins', 'Fst','mean')
  Fstfile[Fstfile$CHR == 'chrX',]$CHR <- 'chr23'
  Fstfile$CHR <- as.numeric(str_split_fixed(Fstfile$CHR,'chr',2)[,2])
  
  # filePNG <-paste(prefix_filename,'genotyping.pop.fst.pdf',sep = '.')
  # filePNG <- paste(prefix_filename,'genotyping.pop.fst.tiff',sep = '.')
  CairoPNG(file=filePNG, width = 1500, height = 500)
  CairoPDF(file = filePNG,
           width = 8, height = 4, onefile = TRUE, family = "Helvetica")
  # CairoTIFF(file = filePNG,
           # width = 800, height = 400, onefile = TRUE, family = "Helvetica", dpi=300)
  colorset <- c('#FF0000', '#FFD700', '#2E8B57', '#7FFFAA', '#6495ED', '#0000FF', '#FF00FF')
  manhattan(Fstfile, chr='CHR', bp='POS', p='Fst', snp='SNP', col=colorset, 
            logp=FALSE, 
            suggestiveline=0.15, 
            genomewideline=FALSE, 
            ylab='Fst', 
            ylim=c(0,0.6), 
            font.lab=4,
            cex.lab=0.8,  
            cex=0.1, 
            chrlabs = c(1:22, "X"))
  dev.off()
  
  
  Fstfile <- read.table(paste(path_data, 'fst/merge.paragraph.genotypes.TIBvsAFR.20k_5k.windowed.weir.fst', sep = ''), header = T, stringsAsFactors = F)
  SNP <- paste(Fstfile[,1], Fstfile[,2], sep = ':')
  Fstfile <- cbind(SNP, Fstfile)
  colnames(Fstfile) <- c('SNP', 'CHR', 'POS','END','Bins', 'Fst','mean')
  Fstfile[Fstfile$CHR == 'chrX',]$CHR <- 'chr23'
  Fstfile$CHR <- as.numeric(str_split_fixed(Fstfile$CHR,'chr',2)[,2])
  
  filePNG <-paste(prefix_filename,'genotyping.pop.fst.afr.fst',sep = '.')
  # CairoPNG(file=filePNG, width = 1500, height = 500)
  CairoPDF(file = filePNG,
           width = 6, height = 3, onefile = TRUE, family = "Helvetica")
  colorset <- c('#FF0000', '#FFD700', '#2E8B57', '#7FFFAA', '#6495ED', '#0000FF', '#FF00FF')
  manhattan(Fstfile, chr='CHR', bp='POS', p='Fst', snp='SNP', col=colorset, 
            logp=FALSE, suggestiveline=FALSE, genomewideline=FALSE, ylab='Fst', ylim=c(0,1), 
            font.lab=4,cex.lab=1.2, main='TIBvsAFR', cex=0.8, chrlabs = c(1:22, "X"))
  dev.off()

###################### Genotyping & POP admixture ######################

# Assign the first argument to prefix
prefix=paste(path_data, '/admixture/merge.paragraph.genotypes.local.reformat.addPopInfo.addFST.LDpruned',sep = '')

# Get individual names in the correct order
labels <- data.frame('ind' <- colnames(data_sv_details_raw@gt)[27:217], 'pop' <- NA)
names(labels)<-c("ind","pop")
labels <- labels[labels$ind%in%data_samples_info$SampleID,]
index_sample <- 1
while (index_sample <= nrow(labels)) {
  labels[index_sample,]$pop <- as.character(data_samples_info[data_samples_info$SampleID==labels[index_sample,]$ind,]$Group_Detail)
  index_sample <- index_sample + 1
}

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit('AFR,HANS,HANN,TIBL,TIBG',",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
maxK=7
tbl<-lapply(2:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual
# tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1200,res=200)
# png(filename = paste(prefix_filename,'genotyping.pop.admixture.include_Biaka', 'Group_Detail', 'png',sep = '.'), width = 2000, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.pop.admixture.include_Biaka','Group_Detail','pdf',sep = '.'), 
    width = 6, height = 3)
  par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,9,1),
      mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
  bp<-barplot(t(as.matrix(tbl[[2]][order(labels$n),])), 
              col=rainbow(n=2),
              xaxt="n", 
              border=NA,
              ylab="K=2",
              yaxt="n",
              space=spaces)
  axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)
  lapply(3:(maxK-1), 
         function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), 
                             col=rainbow(n=x),xaxt="n", 
                             border=NA, 
                             ylab=paste0("K=",x),
                             yaxt="n",
                             space=spaces))
  axis(1,at=c(which(spaces==0.5),
              bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2, 
       tick = F,
       labels=unlist(strsplit('AFR,HANS,HANN,TIBL,TIBG',",")))
dev.off()


source('popcorn-master/R/admixture.R')
source('popcorn-master/R/io.R')

sort_order <- c(data_samples_info[data_samples_info$Group_Detail=='AFR'&data_samples_info$Platform=='NGS',]$SampleID,
                data_samples_info[data_samples_info$Group_Detail=='HANS'&data_samples_info$Platform=='NGS',]$SampleID,
                data_samples_info[data_samples_info$Group_Detail=='HANN'&data_samples_info$Platform=='NGS',]$SampleID,
                data_samples_info[data_samples_info$Group_Detail=='TIBL'&data_samples_info$Platform=='NGS',]$SampleID,
                data_samples_info[data_samples_info$Group_Detail=='TIBG'&data_samples_info$Platform=='NGS',]$SampleID)

plot_admixture_multi <- function(prefix, k, sort_order, flag_label=FALSE){
  pops <- read_fam(paste(prefix,'fam', sep = '.'))
  pops <- pops[pops$fid%in%data_samples_info$SampleID&
                 pops$iid%in%data_samples_info$SampleID,]
  Q <- read_Q_matrix(paste(prefix,k,'Q', sep = '.'))
  Q <- Q[Q$fid%in%data_samples_info$SampleID&
           Q$iid%in%data_samples_info$SampleID,]
  so <- sort_by_cluster(Q)
  QQ <- tidy(Q, pops)
  return(plot_admixture(QQ, label = flag_label, sort.order = sort_order))
}

p6 <- plot_admixture_multi(prefix, 6, sort_order) + theme(plot.margin = unit(c(0,0,0,0),'cm'))
p5 <- plot_admixture_multi(prefix, 5, sort_order) + theme(plot.margin = unit(c(0,0,0,0),'cm'))
p4 <- plot_admixture_multi(prefix, 4, sort_order) + theme(plot.margin = unit(c(0,0,0,0),'cm'))
p3 <- plot_admixture_multi(prefix, 3, sort_order) + theme(plot.margin = unit(c(0,0,0,0),'cm'))
p2 <- plot_admixture_multi(prefix, 2, sort_order) + theme(plot.margin = unit(c(0,0,0,0),'cm'))
pdf(file = paste(prefix_filename,'genotyping.pop.admixture.include_Biaka','Group_Detail','pdf',sep = '.'), width = 6, height = 3)
multiplot(p6, p5, p4, p3, p2, cols = 1)
dev.off()
# ggarrange(p6, p5, p4, ,p3, p2,ncol=1,nrow = 5)
###################### Genotyping & POP njtree ######################
# library(ape)
# # tree <- rtree(n=20)
# # plot(tree, edge.width = 2)
# 
# mat.genotypes <- t(data_genotypes[26:216])
# mat.nj <- nj(dist.gene(mat.genotypes))
# myBoots <- boot.phylo(mat.nj, mat.genotypes, function(xx){nj(dist.gene(xx))}, B=100, mc.cores = 1)


###################### Genotyping & POP treemix ######################

source('plotting_funcs.R')

# png(filename = paste(prefix_filename,'genotyping.pop', 'treemix', 'png',sep = '.'), width = 1000, height = 1200, res = val_res)
pdf(file = paste(prefix_filename,'genotyping.pop','treemix','pdf',sep = '.'), width = 8, height = 8)
prefix <- paste(path_data, '/treemix/merge.paragraph.genotypes.local.reformat.addPopInfo.addFST.noN.LDpruned', sep = '/') 
plot_tree(cex = 0.6, prefix, mbar=FALSE)
dev.off()
# par(mfrow=c(2,3))
# for(edge in 0:5){
#   plot_tree(cex = 0.8, paste0(prefix ,'.',edge))
#   title(paste(edge, "edges"))
# }
# for (edge in 0:5){
#   plot_resid(stem = paste0(prefix, '.', edge), pop_order = 'samples.clust')
# }
pdf(file = paste(prefix_filename,'genotyping.pop.treemix','detailLoc','pdf',sep = '.'), width = 8, height = 8)
# png(filename = paste(prefix_filename,'genotyping.pop.treemix', 'detailLoc', 'png',sep = '.'), width = 1000, height = 1200, res = val_res)
prefix <- paste(path_data, '/treemix/merge.paragraph.genotypes.local.reformat.addPopInfo.addFST.noN.LDpruned.detailLoc', sep = '/') 
plot_tree(cex = 0.5, prefix, mbar=FALSE)
dev.off()
