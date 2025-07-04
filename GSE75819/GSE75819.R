library(GEOquery)
library(tidyverse)
library(limma)
library(edgeR)
library(lumiHumanIDMapping)
library(lumi)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("lumiHumanIDMapping")
BiocManager::install("lumi")

GSE <- "GSE75819"
GEO_file <-  AnnoProbe::geoChina(GSE)

# 样本注释信息
pd <- pData(GEO_file[[1]]) 
row.names(pd) <- pd[,1]

# 探针注释信息
plate <- read_delim("GSE75819_family.soft", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 122,show_col_types = FALSE) 
ids <- plate |> 
  filter(!is.na(`Symbol`)) |> 
  dplyr::select(ID, `Symbol`) |> 
  mutate(
    `Symbol` = str_split(`Symbol`, " /// ", simplify = TRUE)[, 1]
  )

# 表达谱 下载好的格式不需要改变
exp <- read.delim("GSE75819_non-normalized.txt",skip = 5,header = T)
exp <- exp |> 
  dplyr::select(-ProbeID)
write.table(exp,file ='tmp.txt',sep ='\t',quote = F)

# 标准化
x.lumi <- lumiR('tmp.txt',lib.mapping='lumiHumanIDMapping') 
lumi.N.Q <- lumiExpresso(x.lumi)


dat <- exprs(lumi.N.Q)
probeid<-lumi.N.Q@featureData@data
probeid<-probeid[match(rownames(dat),rownames(probeid)),]
identical(rownames(dat),rownames(probeid))
probeid$PROBE_ID->rownames(dat)
dat[1:4,1:4]

#有些探针对应不到基因，去掉symbol为空的
colnames(ids)<-c("probeid","symbol")
ids=ids[ids$symbol!= '',]
dat=dat[rownames(dat) %in% ids$probeid,]
ids=ids[match(rownames(dat),ids$probeid),]
head(ids)  
#在这个示例数据集中probeid的组成有字母有数字，有些只有数字的就要转成character
ids$probeid=as.character(ids$probeid)
identical(rownames(dat),ids$probeid)

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果
dat=dat[ids$probeid,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
range(dat)

# 调整pd的行名顺序与exp列名完全一致（分组需要）
p = identical(rownames(pd),colnames(dat));p
if (!p) {
  dat <- dat[, match(rownames(pd), colnames(dat))]
}

write.table(dat,paste0(GSE,"_expression.xls"),row.names = T,sep = "\t",quote = F)
write.table(pd,paste0(GSE,"_pd.xls"),row.names = F,sep = "\t",quote = F)
