library(GEOquery)
library(tidyverse)
library(limma)
library(edgeR)


GSE <- "GSE255070"
# GEO_file <-  AnnoProbe::geoChina(gse)
GEO_file <- getGEO(GSE,destdir = '.',AnnotGPL = T,getGPL = T)

pd <- pData(GEO_file[[1]])  # 提取数据中的样本临床信息（比如：年龄、性别、是否存活）
exp <- exprs(GEO_file[[1]])  # 提取数据中的样本基因表达矩阵
plate <- fData(GEO_file[[1]])  #提取数据中的平台信息
#plate <- read_delim("GSE53146_family.soft", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 91,show_col_types = FALSE) #导入探针信息



# 整理数据 --------------------------------------------------------------------
# 调整pd的行名顺序与exp列名完全一致（分组需要）
p = identical(rownames(pd),colnames(exp));p
if (!p) {
  exp <- exp[, match(rownames(pd), colnames(exp))]
}


# 探针名替换为基因名 ---------------------------------------------------------------
ids <- plate |> 
  filter(!is.na(miRNA_ID)) |> 
  select(ID,miRNA_ID)

ids$ID <- as.character(ids$ID)

exp <- exp |> 
  data.frame() |> 
  rownames_to_column("ID") |> 
  left_join(ids,by = "ID") |> 
  distinct(miRNA_ID,.keep_all = TRUE) |> 
  filter(!is.na(miRNA_ID) & trimws(miRNA_ID) != "") |> 
  column_to_rownames("miRNA_ID") |> 
  select(-ID)

write.table(exp,paste0(GSE,"_expression.xls"),row.names = T,sep = "\t",quote = F)
# exp <- log2(exp)



# 构建分组 --------------------------------------------------------------------
metadata <- pd |> 
  select(title,sample = geo_accession) |> 
  mutate(condition = case_when(str_detect(title, "control") ~ "Control",
                               TRUE ~ "Treate")) |> 
  select(-title,-sample)

p = identical(rownames(metadata),colnames(exp));p

#过滤掉表达量过低的基因
exp <- exp[rowMeans(exp) > 0, ]


# 差异分析 --------------------------------------------------------------------
group_list = factor(metadata$condition,levels = c("Control","Treate"))
design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(exp)


fit <- lmFit(exp, design) # 拟合模型
constrasts = paste(rev(levels(group_list)),collapse = "-");constrasts
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) # 对比矩阵
fit2=contrasts.fit(fit,cont.matrix) # 应用对比矩阵
fit2=eBayes(fit2)# 贝叶斯平滑
DEG = topTable(fit2, coef=constrasts,adjust="fdr", sort.by="P", n=Inf) # 提取结果
DEG = na.omit(DEG)

# 提取差异显著基因 ----------------------------------------------------------------
psig <- DEG %>% filter(P.Value < 0.05 & abs(logFC) >= 1) %>% mutate(status = case_when(logFC >= 1 ~ "up",logFC <= 1 ~ "down"))
psig_out <- data.frame(gene=rownames(psig), psig[, c("logFC", "P.Value", "adj.P.Val", "status")], exp[rownames(psig),])
write.table(psig_out, paste0(GSE,"_miRNA_difference_psig.xls"), sep="\t", row.names=F, col.names=T, quote=F)

DEG <- DEG |> rownames_to_column("gene")
write.table(DEG, paste0(GSE,"GSE53146_miRNA_difference.xls"), sep="\t", row.names=F, col.names=T, quote=F)

