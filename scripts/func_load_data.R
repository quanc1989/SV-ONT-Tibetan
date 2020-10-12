func.load.svtk.annotsv.local.include_Biaka <- function(path, samplelist){
  
  num_samples <- nrow(samplelist)
  
  num_ont_samples <- nrow(samplelist[samplelist$Platform=='ONT',])
  num_ont_samples_TIB <- nrow(samplelist[samplelist$Platform=='ONT'&samplelist$Races=='TIB',])
  num_ont_samples_HAN <- nrow(samplelist[samplelist$Platform=='ONT'&samplelist$Races=='HAN',])
  
  
  num_ngs_samples_TIB <- nrow(samplelist[samplelist$Platform=='NGS'&samplelist$Races=='TIB',])
  num_ngs_samples_HAN <- nrow(samplelist[samplelist$Platform=='NGS'&samplelist$Races=='HAN',])
  num_ngs_samples_AFR <- nrow(samplelist[samplelist$Platform=='NGS'&samplelist$Races=='Biaka',])
  num_ngs_samples_HUM <- nrow(samplelist[samplelist$Platform=='NGS'&samplelist$Races=='HUM',])
  num_ngs_samples <- nrow(samplelist[samplelist$Platform=='NGS',]) - num_ngs_samples_HUM
  
  group_sv <- c('Singleton','Polymorphic','Major', 'Shared')
  data_samples_sv <- expand.grid(samplelist$SampleID, group_sv)
  colnames(data_samples_sv) <- c('id', 'Type')
  data_samples_sv$Discovery <- 0
  
  data.bed.raw <- read.csv(path, row.names = NULL, header = TRUE, sep = '\t')
  
  column_info <- c( 'SVLEN', 'END', 
                    'SUPP', 'SUPP_RATE','SUPP_NGS','SUPP_RATE_NGS', 
                   'AF_ALL', 'AF_ALL_NGS', 
                   'DP_ALL', 'DP_ALL_NGS', 
                   'MR', 'MR_NGS', 
                   'FST', 'FST_NGS', 'FST_TIBvsHAN', 'FST_TIBGvsHANB', 'FST_TIBGvsHANS', 'FST_TIBGvsHANN', 'FST_TIBLvsHANB',  'FST_TIBLvsHANS', 'FST_TIBLvsHANN',
                   'FST_TIBvsAFR', 'FST_TIBGvsAFR', 'FST_TIBLvsAFR', 'FST_HANvsAFR', 'FST_HANBvsAFR', 'FST_HANNvsAFR', 'FST_HANSvsAFR', 'FST_TIBGvsTIBL', 'FST_HANBvsHANH',
                   'FST_HANBvsHANS', 'FST_HANNvsHANS',
                   'MAXCI_POS', 'MAXCI_END')
  
  for (subgroup in c('HAN', 'TIB')){
    
    index_col <- length(column_info) + 1
    column_info[index_col] <-  paste('SUPP_', subgroup, sep = '')
    column_info[index_col + 1] <-  paste('SUPP_RATE_', subgroup, sep = '')
    column_info[index_col + 2] <-  paste('AF_', subgroup, sep = '')
    column_info[index_col + 3] <-  paste('DP_', subgroup, sep = '')
    column_info[index_col + 4] <-  paste('MR_', subgroup, sep = '')
  }
  
  for (subgroup in c('HAN', 'TIB', 'TIBG', 'TIBL', 'HANS', 'HANN', 'HANB', 'HUM', 'AFR', 'NEA', 'DEN')){
    
    index_col <- length(column_info) + 1
    column_info[index_col] <-  paste('SUPP_', subgroup, '_NGS', sep = '')
    column_info[index_col + 1] <-  paste('SUPP_RATE_', subgroup, '_NGS', sep = '')
    column_info[index_col + 2] <-  paste('AF_', subgroup, '_NGS', sep = '')
    column_info[index_col + 3] <-  paste('DP_', subgroup, '_NGS', sep = '')
    column_info[index_col + 4] <-  paste('MR_', subgroup, '_NGS', sep = '')
  }
  
  for (name_col in column_info){
    data.bed.raw[,name_col] <- as.numeric(str_split_fixed(str_split_fixed(str_extract_all(string = data.bed.raw$INFO, pattern = paste(name_col,'.+;', sep = '=')), '=', 2)[,2], ';',2)[,1])
  }
  
  column_info <- c('SVTYPE', 'CHR2','STRANDS','SEQ')
  
  for (name_col in column_info){
    data.bed.raw[,name_col] <- as.factor(str_split_fixed(str_split_fixed(str_extract_all(string = data.bed.raw$INFO, pattern = paste(name_col,'.+;', sep = '=')), '=', 2)[,2], ';',2)[,1])
  }
  
  tmp <- str_split_fixed(str_extract_all(string=data.bed.raw$INFO, pattern = 'INV_SPAN=.+|NEAREST_TSS=.+;|INTRONIC=.+|DUP_PARTIAL=.+|COPY_GAIN=.+|DUP_LOF=.+|LOF=.+|UTR=.+'),'=',2)
  data.bed.raw[,'GROUP_SV'] <- tmp[,1]
  data.bed.raw[data.bed.raw$GROUP_SV=='NEAREST_TSS','GROUP_SV'] <- 'INTERGENIC'
  data.bed.raw[data.bed.raw$GROUP_SV=='DUP_LOF','GROUP_SV'] <- 'LOF'
  data.bed.raw[,'GROUP_SV'] <- factor(data.bed.raw[,'GROUP_SV'], levels = c('UTR','DUP_PARTIAL','LOF','INV_SPAN','COPY_GAIN','INTRONIC','INTERGENIC'))
  data.bed.raw[,'Gene_svtk'] <- str_split_fixed(tmp[,2],';',2)[,1]
  
  data.bed.raw <- data.bed.raw[,!str_detect(colnames(data.bed.raw),'INFO')]
  
  data.bed.raw$GROUP_LEN <- '50bp-1kb'
  data.bed.raw[abs(data.bed.raw$SVLEN)>1000,]$GROUP_LEN <- '1kb-10kb'
  data.bed.raw[abs(data.bed.raw$SVLEN)>10000,]$GROUP_LEN <- '10kb-100kb'
  data.bed.raw[abs(data.bed.raw$SVLEN)>100000,]$GROUP_LEN <- '100kb-1Mb'
  data.bed.raw[abs(data.bed.raw$SVLEN)>1000000,]$GROUP_LEN <- 'above 1Mb'
  data.bed.raw$GROUP_LEN <- factor(data.bed.raw$GROUP_LEN, levels = c('50bp-1kb', '1kb-10kb', '10kb-100kb', '100kb-1Mb', 'above 1Mb'))
  
  data.bed.raw <- subset(data.bed.raw, SUPP>0)
  
  data.bed.raw$GROUP_SUPP <- NA
  data.bed.raw[data.bed.raw$SUPP==1,]$GROUP_SUPP <- 'Singleton'
  data.bed.raw[data.bed.raw$SUPP>1,]$GROUP_SUPP <- 'Polymorphic'
  data.bed.raw[data.bed.raw$SUPP>num_ont_samples/2,]$GROUP_SUPP <- 'Major'
  data.bed.raw[data.bed.raw$SUPP==num_ont_samples,]$GROUP_SUPP <- 'Shared'
  
  data.bed.raw$GROUP_SUPP <- factor(data.bed.raw$GROUP_SUPP, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  data.bed.raw$GROUP_SUPP_TIB <- NA
  data.bed.raw[data.bed.raw$SUPP_TIB==1,]$GROUP_SUPP_TIB <- 'Singleton'
  data.bed.raw[data.bed.raw$SUPP_TIB>1,]$GROUP_SUPP_TIB <- 'Polymorphic'
  data.bed.raw[data.bed.raw$SUPP_TIB>num_ont_samples_TIB/2,]$GROUP_SUPP_TIB <- 'Major'
  data.bed.raw[data.bed.raw$SUPP_TIB==num_ont_samples_TIB,]$GROUP_SUPP_TIB <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_TIB <- factor(data.bed.raw$GROUP_SUPP_TIB, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  data.bed.raw$GROUP_SUPP_HAN <- NA
  data.bed.raw[data.bed.raw$SUPP_HAN==1,]$GROUP_SUPP_HAN <- 'Singleton'
  data.bed.raw[data.bed.raw$SUPP_HAN>1,]$GROUP_SUPP_HAN <- 'Polymorphic'
  data.bed.raw[data.bed.raw$SUPP_HAN>num_ont_samples_HAN/2,]$GROUP_SUPP_HAN <- 'Major'
  data.bed.raw[data.bed.raw$SUPP_HAN==num_ont_samples_HAN,]$GROUP_SUPP_HAN <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_HAN <- factor(data.bed.raw$GROUP_SUPP_HAN, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  
  data.bed.raw$GROUP_SUPP_NGS <- NA
  
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS) & data.bed.raw$SUPP_NGS==1,]$GROUP_SUPP_NGS <- 'Singleton'
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS) &data.bed.raw$SUPP_NGS>1,]$GROUP_SUPP_NGS <- 'Polymorphic'
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS) &data.bed.raw$SUPP_NGS>num_ngs_samples/2,]$GROUP_SUPP_NGS <- 'Major'
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS) &data.bed.raw$SUPP_NGS==num_ngs_samples,]$GROUP_SUPP_NGS <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_NGS <- factor(data.bed.raw$GROUP_SUPP_NGS, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  
  data.bed.raw$GROUP_SUPP_TIB_NGS <- NA
  data.bed.raw[!is.na(data.bed.raw$SUPP_TIB_NGS) & data.bed.raw$SUPP_TIB_NGS==1,]$GROUP_SUPP_TIB_NGS <- 'Singleton'
  data.bed.raw[!is.na(data.bed.raw$SUPP_TIB_NGS) &data.bed.raw$SUPP_TIB_NGS>1,]$GROUP_SUPP_TIB_NGS <- 'Polymorphic'
  data.bed.raw[!is.na(data.bed.raw$SUPP_TIB_NGS) &data.bed.raw$SUPP_TIB_NGS>num_ngs_samples_TIB/2,]$GROUP_SUPP_TIB_NGS <- 'Major'
  data.bed.raw[!is.na(data.bed.raw$SUPP_TIB_NGS) &data.bed.raw$SUPP_TIB_NGS==num_ngs_samples_TIB,]$GROUP_SUPP_TIB_NGS <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_TIB_NGS <- factor(data.bed.raw$GROUP_SUPP_TIB_NGS, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  data.bed.raw$GROUP_SUPP_HAN_NGS <- NA
  data.bed.raw[!is.na(data.bed.raw$SUPP_HAN_NGS) & data.bed.raw$SUPP_HAN_NGS==1,]$GROUP_SUPP_HAN_NGS <- 'Singleton'
  data.bed.raw[!is.na(data.bed.raw$SUPP_HAN_NGS) & data.bed.raw$SUPP_HAN_NGS>1,]$GROUP_SUPP_HAN_NGS <- 'Polymorphic'
  data.bed.raw[!is.na(data.bed.raw$SUPP_HAN_NGS) & data.bed.raw$SUPP_HAN_NGS>num_ngs_samples_HAN/2,]$GROUP_SUPP_HAN_NGS <- 'Major'
  data.bed.raw[!is.na(data.bed.raw$SUPP_HAN_NGS) & data.bed.raw$SUPP_HAN_NGS==num_ngs_samples_HAN,]$GROUP_SUPP_HAN_NGS <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_HAN_NGS <- factor(data.bed.raw$GROUP_SUPP_HAN_NGS, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  data.bed.raw$GROUP_SUPP_AFR_NGS <- NA
  data.bed.raw[!is.na(data.bed.raw$SUPP_AFR_NGS) & data.bed.raw$SUPP_AFR_NGS==1,]$GROUP_SUPP_AFR_NGS <- 'Singleton'
  data.bed.raw[!is.na(data.bed.raw$SUPP_AFR_NGS) & data.bed.raw$SUPP_AFR_NGS>1,]$GROUP_SUPP_AFR_NGS <- 'Polymorphic'
  data.bed.raw[!is.na(data.bed.raw$SUPP_AFR_NGS) & data.bed.raw$SUPP_AFR_NGS>num_ngs_samples_AFR/2,]$GROUP_SUPP_AFR_NGS <- 'Major'
  data.bed.raw[!is.na(data.bed.raw$SUPP_AFR_NGS) & data.bed.raw$SUPP_AFR_NGS==num_ngs_samples_AFR,]$GROUP_SUPP_AFR_NGS <- 'Shared'
  
  data.bed.raw$GROUP_SUPP_AFR_NGS <- factor(data.bed.raw$GROUP_SUPP_AFR_NGS, levels = c('Singleton','Polymorphic','Major', 'Shared'))
  
  data.bed.raw$GROUP_POP <- 'TIB'
  data.bed.raw[data.bed.raw$SUPP_HAN>0,]$GROUP_POP <- 'HAN'
  data.bed.raw[data.bed.raw$SUPP_HAN>0&data.bed.raw$SUPP_TIB>0,]$GROUP_POP <- 'Shared'
  data.bed.raw$GROUP_POP <- factor(data.bed.raw$GROUP_POP, levels = c('Shared','HAN','TIB'))

  data.bed.raw$GROUP_POP_NGS <- 'Shared'
  
  for (subgroup  in c('HAN', 'TIB', "AFR")) {
    condition_han <- !is.na(data.bed.raw$SUPP_NGS)&data.bed.raw[, paste('SUPP_',subgroup,'_NGS', sep = '')]==data.bed.raw$SUPP_NGS
    data.bed.raw[condition_han,]$GROUP_POP_NGS <- subgroup
  }
  
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS)&
                 data.bed.raw$SUPP_AFR_NGS>0&
                 data.bed.raw$SUPP_HAN_NGS>0&
                 data.bed.raw$SUPP_TIB_NGS>0,]$GROUP_POP_NGS <- 'All'
  
  data.bed.raw$GROUP_POP_NGS <- factor(data.bed.raw$GROUP_POP_NGS, levels = c('All','Shared','AFR','HAN','TIB'))
  
  data.bed.raw$GROUP_POP_DETAIL_NGS <- 'Shared'
  
  for (subgroup  in c('HANN', 'HANS', 'TIBG', 'TIBL', "AFR")) {
    condition_han <- !is.na(data.bed.raw$SUPP_NGS)&data.bed.raw[, paste('SUPP_',subgroup,'_NGS', sep = '')]==data.bed.raw$SUPP_NGS
    data.bed.raw[condition_han,]$GROUP_POP_DETAIL_NGS <- subgroup
  }
  
  data.bed.raw[!is.na(data.bed.raw$SUPP_NGS)&
                 data.bed.raw$SUPP_AFR_NGS>0&
                 data.bed.raw$SUPP_HANN_NGS>0&
                 data.bed.raw$SUPP_HANS_NGS>0&
                 data.bed.raw$SUPP_TIBG_NGS>0&
                 data.bed.raw$SUPP_TIBL_NGS>0,]$GROUP_POP_DETAIL_NGS <- 'All'

  data.bed.raw$GROUP_POP_DETAIL_NGS <- factor(data.bed.raw$GROUP_POP_DETAIL_NGS, levels = c('All','Shared','AFR','TIBG','TIBL','HANN','HANS'))
  
  data_samples_genotypes <- data.bed.raw[,colnames(data.bed.raw)%in%samplelist$SampleID]
  
  for (label_sample in samplelist$SampleID) {
    
    gt <- str_split_fixed(data.bed.raw[,label_sample],':',3)
    tmp <- rep(0,nrow(gt))
    
    tmp[gt[,1]=='./.'] <- -1
    tmp[gt[,1]=='0/0'] <- 0
    tmp[gt[,1]=='0/1'] <- 1
    tmp[gt[,1]=='1/1'] <- 2
    
    data_samples_genotypes[,label_sample] <- tmp
    
    for (svtype in levels(data.bed.raw$GROUP_SUPP)){
      data_samples_sv[data_samples_sv$id==label_sample&data_samples_sv$Type==svtype,]$Discovery <- length(which(tmp>0&data.bed.raw$GROUP_SUPP==svtype))
    }
    
    data.bed.raw <- data.bed.raw[,!str_detect(colnames(data.bed.raw),label_sample)]
  }
  
  colnames(data.bed.raw)[2] <- 'Chrom'
  colnames(data.bed.raw)[3] <- 'POS'
  colnames(data.bed.raw)[7] <- 'SVID'
  
  rownames(data_samples_genotypes) <- data.bed.raw$SVID
  data.bed.raw$POS <- as.numeric(str_split_fixed(data.bed.raw$SVID,'_',4)[,2])
  data.bed.raw$SVTYPE <- as.character(data.bed.raw$SVTYPE)
  data.bed.raw[data.bed.raw$SVTYPE=='BND',]$SVTYPE <- 'TRA'
  data.bed.raw$SVTYPE <- factor(data.bed.raw$SVTYPE, levels = c('DEL', 'INS', 'DUP', 'INV', 'TRA'))
  
  data.bed.raw$AF_Database <- 0
  data.bed.raw$Database <- NA
  
  order_names <- c('DGV_GAIN_Frequency','DGV_LOSS_Frequency','GD_POPMAX_AF','X1000g_max_AF','DDD_DEL_Frequency','DDD_DUP_Frequency','IMH_AF')
  database_names <- c('DGV','DGV','GnomAD','1000Genome','DDD','DDD','IMH')
  
  data.database.matrix <- data.bed.raw[,order_names]
  index_max <- max.col(data.database.matrix)
  
  data.bed.raw$AF_Database <- data.database.matrix[cbind(1:nrow(data.database.matrix),index_max)] 
  # data.bed.raw$Database <- colnames(data.database.matrix)[index_max]
  data.bed.raw$Database <- database_names[index_max]
  data.bed.raw[data.bed.raw$AF_Database<0,]$Database <- 'None'
  data.bed.raw[data.bed.raw$AF_Database<0,]$AF_Database <- 0
  data.bed.raw$Database <- factor(data.bed.raw$Database)
  
  data.bed.raw[data.bed.raw$AF_Database>1, 'AF_Database'] <- 1
  
  data.bed.raw <- data.bed.raw[,!str_detect(colnames(data.bed.raw),'SV.end')]
  data.bed.raw <- data.bed.raw[,!str_detect(colnames(data.bed.raw),'SV.type')]
  data.bed.raw <- data.bed.raw[,!str_detect(colnames(data.bed.raw),'SV.length')]
  
  data.bed.raw$GCcontent_left <- as.numeric(as.character(data.bed.raw$GCcontent_left))
  data.bed.raw$GCcontent_right <- as.numeric(as.character(data.bed.raw$GCcontent_right))
  
  tryCatch({
    data.bed.raw$RepClass <- as.character(data.bed.raw$RepClass)
    data.bed.raw[data.bed.raw$RepClass==''|is.na(data.bed.raw$RepClass),]$RepClass <- 'None'
    data.bed.raw[str_detect(data.bed.raw$RepClass, 'Unknown'),]$RepClass <- 'Other'
    for(subclass in c( 'Other','RNA', 'RC', 'DNA', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple_repeat', 'Low_complexity' )){
      data.bed.raw[str_detect(data.bed.raw$RepClass, subclass),]$RepClass <- subclass
    }
  },error=function(e){
  })
  
  data.bed.raw$RepClass <-factor(data.bed.raw$RepClass, levels = c('Low_complexity', 'Simple_repeat', 'Satellite','SINE','LINE','LTR','DNA','RNA','RC','Other','None'))
  
  data.bed.raw$RepFamily <- str_split_fixed(data.bed.raw$RepFamily, '/', 2)[,1]
  data.bed.raw[data.bed.raw$RepFamily=='',]$RepFamily <- 'None'
  data.bed.raw$RepFamily <-factor(data.bed.raw$RepFamily)
  
  
  tryCatch({
    data.bed.raw$Repeats_type_left <- as.character(data.bed.raw$Repeats_type_left)
    data.bed.raw[data.bed.raw$Repeats_type_left==''|is.na(data.bed.raw$Repeats_type_left),]$Repeats_type_left <- 'None'
    data.bed.raw[str_detect(data.bed.raw$Repeats_type_left, 'Unknown'),]$Repeats_type_left <- 'Other'
    for(subclass in c( 'Other','RNA', 'RC', 'DNA', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple_repeat', 'Low_complexity' )){
      data.bed.raw[str_detect(data.bed.raw$Repeats_type_left, subclass),]$Repeats_type_left <- subclass
    }
  },error=function(e){
  })
  
  data.bed.raw$Repeats_type_left <-factor(data.bed.raw$Repeats_type_left,levels = c('Low_complexity', 'Simple_repeat', 'Satellite','SINE','LINE','LTR','DNA','RNA','RC','Other','None'))
  
  tryCatch({
    data.bed.raw$Repeats_type_right <- as.character(data.bed.raw$Repeats_type_right)
    data.bed.raw[data.bed.raw$Repeats_type_right==''|is.na(data.bed.raw$Repeats_type_right),]$Repeats_type_right <- 'None'
    data.bed.raw[str_detect(data.bed.raw$Repeats_type_right, 'Unknown'),]$Repeats_type_right <- 'Other'
    for(subclass in c( 'Other','RNA', 'RC', 'DNA', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple_repeat', 'Low_complexity' )){
      data.bed.raw[str_detect(data.bed.raw$Repeats_type_right, subclass),]$Repeats_type_right <- subclass
    }
  },error=function(e){
  })
  
  data.bed.raw$Repeats_type_right <-factor(data.bed.raw$Repeats_type_right,levels = c('Low_complexity', 'Simple_repeat', 'Satellite','SINE','LINE','LTR','DNA','RNA','RC','Other','None'))
  
  data.bed.raw$TADcoordinates <- as.character(data.bed.raw$TADcoordinates)
  data.bed.raw$TADcoordinates2 <- as.character(data.bed.raw$TADcoordinates2)
  data.bed.raw$dbVar_variant <- as.character(data.bed.raw$dbVar_variant)
  data.bed.raw$promoters <- as.character(data.bed.raw$promoters)
  data.bed.raw$IMH_ID_others <- as.character(data.bed.raw$IMH_ID_others)
  return(list(data_sv=data.bed.raw,data_samples=data_samples_sv,data_genotypes=data_samples_genotypes))
}
