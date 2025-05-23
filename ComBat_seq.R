library(sva)
library(tinyarray)
library(readr)
library(tidyverse)
# 导入文件
disease_exp <- read_delim("GSE125512_expression.xls", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
healthy_exp <- read_delim("GSE202709_expression.xls", delim = "\t", escape_double = FALSE, trim_ws = TRUE) |> 
  select(1:5)

# 表达谱
exp <- healthy_exp |> 
  inner_join(disease_exp,by = c("Gene")) |> 
  column_to_rownames("Gene")

# 分组文件
group <- c(rep("healthy",4),rep("disease",22))
batch <- c(rep("batch1",4),rep("batch2",22))
coldata <- data.frame(group,batch)
rownames(coldata) <- colnames(exp)

# 去除批次效应
expr_combat <- sva::ComBat_seq(as.matrix(exp), batch = coldata$batch)

# PCA
tinyarray::draw_pca(exp = expr_combat,group_list = factor(coldata$group))

# expr_combat <- ComBat(dat = expr_combat, batch = coldata$batch,par.prior=TRUE)
# tinyarray::draw_pca(exp = expr_combat,group_list = factor(coldata$group))

write.table(coldata,"coldata.xls",sep = "\t",quote = F,row.names = T)
write.table(expr_combat,"gene_count_matrix.xls",sep = "\t",quote = F,row.names = T)
