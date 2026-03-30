library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)

setwd('/picb/neurosys/chenrenrui/Mouse_ST/A molecularly defined and spatially resolved cell atlas of the whole mouse brain/')

chip = list.files('./data_filter/neocortex/')
chip = str_split(chip,".txt",simplify = T)[,1]

cell_pvalue_dir <- list.files('./subclass_10K/neocortex_sample_subclass/cell_pvalue/')
chip <- intersect(chip,cell_pvalue_dir)

Merfish_brain_cell_type <- read.delim('./Merfish_brain_cell_type_new.txt')


cl <- 62
registerDoParallel(cl)
foreach (chip.1=chip[101:155],.packages = c("data.table","stringr")) %dopar% {
  file_names_list <- list.files(paste('./subclass_10K/neocortex_sample_subclass/cell_pvalue/',chip.1,sep = ''))
  dir.create(paste('./subclass_10K/neocortex_sample_subclass/merge_region/',chip.1,sep = ''))
  for (sample.1 in 1:500) {
    file_names_list_sample <- file_names_list[grep(paste('sample_',sample.1,'_',sep = ''),file_names_list)]
    dir.create(paste('./subclass_10K/neocortex_sample_subclass/merge_region/',chip.1,'/','sample_',sample.1,sep = ''))
    
    for (file_name in file_names_list_sample) {
      mouse_sig_p.total <- read.csv(paste('./subclass_10K/neocortex_sample_subclass/cell_pvalue/',chip.1,'/',file_name,sep = ''))
      mouse_sig_p.total <- mouse_sig_p.total[mouse_sig_p.total$p<0.05,]
      
      if (nrow(mouse_sig_p.total)>0) {
        
        cell_sliding_windows_result_p <- mouse_sig_p.total
        
        
        cell_sliding_windows_result_p$xstart <- as.numeric(str_split(cell_sliding_windows_result_p$loc,"[_]",simplify = T)[,1])
        cell_sliding_windows_result_p$ystart <- as.numeric(str_split(cell_sliding_windows_result_p$loc,"[_]",simplify = T)[,2])
        
        cell_sliding_windows_result_p$xend <- cell_sliding_windows_result_p$xstart + 0.4
        cell_sliding_windows_result_p$yend <- cell_sliding_windows_result_p$ystart + 0.4
        
        
        load(paste('./subclass_10K/neocortex_sample/cell_window_layer/',chip.1,'/',
                   str_split(file_name,"[_]",simplify = T)[,3],'_',
                   str_split(file_name,"[_]",simplify = T)[,4],
                   '.RData',sep = ''))
        rm(file_tmp_layer)
        
        cell_windows_layer_d <- cell_windows_layer_d_layer
        cell_windows_layer_d$loc_1 <- paste(cell_windows_layer_d$layer,cell_windows_layer_d$loc,sep='_')
        
        Neruon_enrich <- cell_sliding_windows_result_p[cell_sliding_windows_result_p$p<0.05,]
        
        Neruon_enrich$loc_1 <- paste(Neruon_enrich$layer,Neruon_enrich$loc,sep = '_')
        
        Neruon_subclass <- unique(Neruon_enrich$subclass)
        
        for (subclass in unique(Neruon_subclass)) {
          temp_Neruon_enrich <- Neruon_enrich[Neruon_enrich$subclass %in% subclass,]
          temp_Neruon_enrich_cell_id <- cell_windows_layer_d[cell_windows_layer_d$loc_1 %in% temp_Neruon_enrich$loc_1,]
          
          temp_Neruon_enrich_cell_list <- list()
          a = 1
          for (i in unique(temp_Neruon_enrich_cell_id$loc_1)) {
            temp1 <- temp_Neruon_enrich_cell_id[loc_1 %in% i]
            temp_Neruon_enrich_cell_list[[a]] <- temp1$cell_label
            names(temp_Neruon_enrich_cell_list)[a] <- unique(paste(temp1$layer,temp1$loc,sep = '_'))
            a = a + 1
          }
          
          
          # 初始化整合后的区域列表
          merged_regions <- list()
          
          # 遍历每个区域
          for (i in 1:length(temp_Neruon_enrich_cell_list)) {
            region <- temp_Neruon_enrich_cell_list[[i]]
            
            # 检查当前区域是否与已整合的区域有重叠
            overlap <- sapply(merged_regions, function(x) any(intersect(x, region)!=0))
            
            overlap.loc <- which(overlap==TRUE)
            
            if (any(overlap)) {
              if (length(overlap.loc)>1) {
                temp1 <- union(unlist(merged_regions[overlap.loc]), region)
                merged_regions <- c(merged_regions, list(temp1))
                merged_regions <- merged_regions[-c(overlap.loc)]
              }else {
                merged_regions[[overlap.loc]] <- union(unlist(merged_regions[overlap.loc]), region)
              }
              # 如果有重叠的区域，则将当前区域与重叠的区域合并
              
            } else {
              # 如果没有重叠的区域，则将当前区域添加到整合后的区域列表中
              merged_regions <- c(merged_regions, list(region))
            }
          }
          
          temp_Neruon_merged_regions <- merged_regions
          
          temp_Neruon_merged_regions_table <- list()
          a = 1
          for (i in 1:length(temp_Neruon_merged_regions)) {
            temp1 <- temp_Neruon_merged_regions[[i]]
            temp1 <- unique(temp_Neruon_enrich_cell_id$loc_1[temp_Neruon_enrich_cell_id$cell_label %in% temp1])
            temp1 <- temp_Neruon_enrich[temp_Neruon_enrich$loc_1 %in% temp1,]
            temp1$merge_regions <- paste('regions',i,sep = '_')
            temp_Neruon_merged_regions_table[[a]] <- temp1
            names(temp_Neruon_merged_regions_table)[a] <- unique(temp1$layer)
            a = a + 1
          }
          temp_Neruon_merged_regions_table <- rbindlist(temp_Neruon_merged_regions_table)
          temp_Neruon_merged_regions_table$chip <- str_split(file_name,"[_]",simplify = T)[,3]
          temp_Neruon_merged_regions_table <- temp_Neruon_merged_regions_table[,c(13,6,1,12,7:11,5)]
          temp_Neruon_merged_regions_table$sample <- paste('sample_',sample.1,sep = '') 
          temp_Neruon_merged_regions_table <- temp_Neruon_merged_regions_table[,c(1,11,2:10)]
          
          
          write.table(temp_Neruon_merged_regions_table,
                      file = paste('./subclass_10K/neocortex_sample_subclass/merge_region/',chip.1,'/',
                                   'sample_',sample.1,'/',
                                   gsub('_cell_sliding_window_result_p.csv','',file_name),'_',
                                   'merged_rgions_loc.txt',sep = ''),
                      quote = F,row.names = F,sep = '\t')
          
          
          file_tmp <- read.csv(paste('./subclass_10K/neocortex_sample_subclass/data_sample/',chip.1,'/',
                                       str_split(file_name,"[_]",simplify = T)[1],'_',
                                       str_split(file_name,"[_]",simplify = T)[2],'_',
                                       str_split(file_name,"[_]",simplify = T)[3],'_',
                                       str_split(file_name,"[_]",simplify = T)[4],
                                       '.txt',sep = ''))
    
          file_region <- read.delim(paste('./data_filter/neocortex/',
                                        str_split(file_name,"[_]",simplify = T)[3],'.txt',sep = ''))
          file_region <- file_region[,c('cell_label','ccf_region_name'),]
          file_region <- file_region[!duplicated(file_region),]
          
          
          file_tmp <- merge(file_tmp,Merfish_brain_cell_type[,c('subclass','cell_Neuron_type')],by='subclass')
          file_tmp <- merge(file_tmp,file_region,by='cell_label')
          
          temp_Neruon_merged_regions_table_cell_id <- list()
          a = 1
          for (i in 1:length(temp_Neruon_merged_regions)) {
            temp1 <- temp_Neruon_merged_regions[[i]]
            temp1 <- file_tmp[file_tmp$cell_label %in% temp1,]
            
            temp1 <- data.frame(slide =  str_split(file_name,"[_]",simplify = T)[3],
                                sample = paste('sample_',sample.1,sep = ''),
                                layer = unique(temp1$layer),
                                region = paste(unique(temp1$ccf_region_name),collapse = ','),
                                merge_regions = paste('regions',i,sep = '_'),
                                total_cell_num = nrow(temp1),
                                subclass = subclass,
                                enrich_subclass_cell_ids_num = length(temp1$cell_label[temp1$subclass==subclass]),
                                enrich_subclass_cell_ids = paste(temp1$cell_label[temp1$subclass==subclass],collapse = ','),
                                Glut_Neruon_cell_ids_num = length(temp1$cell_label[temp1$cell_Neuron_type=='Glut']),
                                Glut_Neruon_cell_ids = paste(temp1$cell_label[temp1$cell_Neuron_type=='Glut'],collapse = ','),
                                GABA_Neruon_cell_ids_num = length(temp1$cell_label[temp1$cell_Neuron_type=='GABA']),
                                GABA_Neruon_cell_ids = paste(temp1$cell_label[temp1$cell_Neuron_type=='GABA'],collapse = ','))
            temp_Neruon_merged_regions_table_cell_id[[a]] <- temp1
            names(temp_Neruon_merged_regions_table_cell_id)[a] <- unique(temp1$layer)
            a = a + 1
          }
          temp_Neruon_merged_regions_table_cell_id <- rbindlist(temp_Neruon_merged_regions_table_cell_id)
          write.table(temp_Neruon_merged_regions_table_cell_id,
                      file = paste('./subclass_10K/neocortex_sample_subclass/merge_region/',chip.1,'/','sample_',sample.1,'/',
                                   gsub('_cell_sliding_window_result_p.csv','',file_name),'_',
                                   'merged_regions_table_cell_id.txt',sep = ''),
                      sep = '\t',quote = F,row.names = F)
          
          print(paste(gsub('_sliding_window_result_p.csv','',file_name),
                      'merge region:',nrow(temp_Neruon_merged_regions_table_cell_id)))
          
          rm(list = ls()[grep('temp',ls())])
          rm(file_tmp)
          rm(file_region)
          rm(Neruon_enrich)
          rm(merged_regions)
          rm(mouse_sig_p.total)
          rm(cell_windows_layer_d)
          rm(cell_windows_layer_d_layer)
          rm(cell_sliding_windows_result_p)
          
        }
      }
      
    }
  }
  

  
}
stopCluster(cl)
