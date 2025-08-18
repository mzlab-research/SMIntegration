
library(tibble)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(tidyr)
library(readr)
args <- commandArgs(T)
projectdir<- args[1]
outdir<-projectdir
input_m<-file.path(projectdir,"metab.rds")
data_m <- readRDS(file = input_m)
input_t<-file.path(projectdir,"trans.rds")
data_t <- readRDS(file = input_t)

t=as.data.frame(ft@assays$Spatial$counts)
m=as.data.frame(fm@assays$Spatial$counts)


t_long <- t %>%
  rownames_to_column(var = "geneID") %>%  # 将行名转换为列
  pivot_longer(
    cols = -geneID,  # 选择需要转换的列，排除geneID列
    names_to = "x_y",  # 新列的名称
    values_to = "value"  # 值列的名称
  )
t_long$x_y <- gsub("^sample:", "", t_long$x_y)
t_long <-t_long %>%
  separate(col = x_y, into = c("x", "y"), sep = "_")
colnames(t_long)=c("geneID","x","y","MIDCount")
t_long$geneID<-as.character(t_long$geneID)
t_long$x<-as.numeric(t_long$x)
t_long$y<-as.numeric(t_long$y)
t_long$MIDCount<-as.numeric(t_long$MIDCount)
write.table(t_long,file = file.path(outdir,"trans_bin_filter.txt"),sep="\t",row.names=F,quote=F)

m_long <- m %>%
  rownames_to_column(var = "geneID") %>%  # 将行名转换为列
  pivot_longer(
    cols = -geneID,  # 选择需要转换的列，排除geneID列
    names_to = "x_y",  # 新列的名称
    values_to = "value"  # 值列的名称
  )
m_long$x_y <- gsub("^sample:", "", m_long$x_y)
m_long <-m_long %>%
  separate(col = x_y, into = c("x", "y"), sep = "_")
colnames(m_long)=c("mz","x","y","Intensity")
m_long$mz<-as.character(m_long$mz)
m_long$x<-as.numeric(m_long$x)
m_long$y<-as.numeric(m_long$y)
m_long$Intensity<-as.numeric(m_long$Intensity)
write.table(m_long,file = file.path(outdir,"metab_bin_filter.txt"),sep="\t",row.names=F,quote=F)