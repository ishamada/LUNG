##############################
# TCGA-HNSC RNA-Seq Analysis Script
# Author: Nora Mohamed
# Purpose: Download clinical and RNA-Seq counts for selected cases,
#          remove metastatic samples, and prepare counts matrix
##############################

##############################
# 1. Load Required Libraries
##############################
library(ggplot2)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(sva)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(123)

##############################
# 2. Load Cases and Prepare Clinical Data
##############################

# Load your selected cases
cases <- read.csv("/home/norah/Downloads/new task/cases.tsv", sep = '\t')

# Query clinical data from GDC (if internet is available)
# If not, clinical can be extracted later from the prepared SE object
# clinical_hnsc <- GDCquery_clinic("TCGA-HNSC")

##############################
# 3. Build Query for RNA-Seq Data (all selected cases)
##############################

# First, build a query for all samples in the cases list
query_TCGA_subset <- GDCquery(
  project = 'TCGA-HNSC',
  data.category = 'Transcriptome Profiling',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  barcode = cases$submitter_id,   # Only selected patient barcodes
  access = 'open'
)

##############################
# 4. Download Data
##############################
GDCdownload(query_TCGA_subset, directory = "/home/norah/Downloads/new task")

##############################
# 5. Prepare SummarizedExperiment Object
##############################
hnsc.tcga.data <- GDCprepare(
  query_TCGA_subset,
  summarizedExperiment = TRUE,
  directory = "/home/norah/Downloads/new task"
)

# Extract clinical data from the prepared object
clinical <- as.data.frame(colData(hnsc.tcga.data))

# Optional: match clinical rows to cases file
final_barcodes <- clinical[clinical$patient %in% cases$submitter_id, ]
barcodes <- final_barcodes$barcode

# Check distribution of sample types
table(final_barcodes$definition)

##############################
# 6. Filter Out Metastatic Samples
##############################
hnsc_no_metastatic <- hnsc.tcga.data[, hnsc.tcga.data$definition != "Metastatic"]
table(hnsc_no_metastatic$definition)  # Verify removal

##############################
# 7. Extract Counts Matrix
##############################
counts_data <- assay(hnsc_no_metastatic)  # genes x samples

# Optional: save the filtered counts matrix for future use
saveRDS(counts_data, "/home/norah/Downloads/new task/counts_data_filtered.rds")
saveRDS(hnsc_no_metastatic, "/home/norah/Downloads/new task/hnsc_no_metastatic.rds")

##############################
# 8. Optional: Explore Clinical Variables
##############################
# For example, check prior treatment distribution
table(final_barcodes$prior_treatment)

##############################
# Script End
##############################

# ===============================================
# DEA for mitochondrial genes in HNSCC
# Using filtered cases metadata
# ===============================================

# 1️⃣ Load libraries
library(DESeq2)
library(sva)
library(biomaRt)
library(ggplot2)
library(pheatmap)

#### hello ####

# 2️⃣ Load data
counts <- read.delim("~/Downloads/new task/TCGA_HNSC_Raw_Counts_Matrix.tsv", row.names = 1, check.names = FALSE)
metadata <- read.csv("~/Downloads/new task/cases_filtered_metadata.csv", row.names = 1)  # rownames = sample IDs

# تأكد أن الأعمدة في counts والـ metadata متطابقة
counts <- counts[, rownames(metadata)]

# 3️⃣ DEA setup: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ origin_of_tissue + sample_type)  # sample_type = Primary/Metastatic

# Exclude metastasis if not already filtered
dds <- dds[, dds$sample_type == "Primary"]

# 4️⃣ Pre-filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 5️⃣ Batch correction using sva
# Create model matrices
mod <- model.matrix(~ sample_type, data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds))
svobj <- sva(counts(dds), mod, mod0)

# Add surrogate variables to design
colData(dds) <- cbind(colData(dds), svobj$sv)
design(dds) <- as.formula(paste("~", paste(c(paste0("sv", 1:ncol(svobj$sv)), "sample_type"), collapse = "+")))

# 6️⃣ Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("sample_type", "Primary", "Metastatic"))  # adjust if needed

# 7️⃣ Convert Ensembl IDs to Gene Symbols using BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- rownames(res)
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = genes,
                  mart = ensembl)

res$gene_symbol <- gene_map$hgnc_symbol[match(rownames(res), gene_map$ensembl_gene_id)]

# 8️⃣ Filter mitochondrial genes (مثلاً قائمة mt_genes موجودة عندك)
# mt_genes <- c("MT-CO1","MT-CO2",...)  # ضع قائمة الجينات الميتوكوندريا المطلوبة
res_mt <- res[res$gene_symbol %in% mt_genes, ]

# 9️⃣ Save results
write.csv(as.data.frame(res_mt), "DEA_mitochondrial_HNSC.csv")

# 🔟 PCA plot for outlier detection
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("sample_type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() +
  ggtitle("PCA - Primary HNSCC samples")

..................................................................

# 1️⃣ اقرأ ملف metadata بعد الفلترة
filtered_cases <- read.csv("~/Downloads/new task/cases_filtered_metadata.csv", row.names = 1)

# 2️⃣ تأكد أن barcodes متوافقة مع الأعمدة
# لو الأعمدة في counts_data هي barcodes، لازم نرتب filtered_cases بنفس الترتيب
filtered_cases <- filtered_cases[match(colnames(counts_data), filtered_cases$barcode), ]

# تحقق من الترتيب
all(colnames(counts_data) == filtered_cases$barcode)





cases_filtered <- read.csv("~/Downloads/new task/cases_filtered_metadata.csv", sep="\t")

# أضيفي condition مباشرة
cases_filtered$condition <- ifelse(
  grepl("Normal", cases_filtered$`barcode.definition.tissue_or_organ_of_origin`, ignore.case = TRUE),
  "Normal",
  "Tumor"
)


table(cases_filtered$condition)  # لازم يطلع Normal + Tumor



#############################
# TCGA-HNSC Filtered Cases Preparation
# Author: Nora Mohamed
# Purpose: Read original CSV, split columns if needed, add condition,
#          remove Metastatic samples, ready for DEA
##############################

library(dplyr)
library(tidyr)

# 1️⃣ اقرأ CSV الأصلي
cases_filtered <- read.csv("~/Downloads/new task/cases_filtered_metadata.csv", sep = ",", stringsAsFactors = FALSE)

# 2️⃣ لو العمود فيه كل المعلومات في خلية واحدة (اختياري)
# لاحظي: لو كل شيء مفصول أصلاً في أعمدة، ممكن تتخطى هذه الخطوة
if(ncol(cases_filtered) == 1){
  cases_filtered <- separate(
    cases_filtered,
    col = 1,
    into = c("barcode", "definition", "origin"),
    sep = ",",
    extra = "merge",
    fill = "right"
  )
}

# 3️⃣ تنظيف النصوص من فراغات إضافية
cases_filtered <- cases_filtered %>%
  mutate(across(c(definition, tissue_or_organ_of_origin, barcode), ~trimws(.)))

# 4️⃣ إضافة condition: Tumor أو Normal
cases_filtered <- cases_filtered %>%
  mutate(condition = ifelse(grepl("Normal", definition, ignore.case = TRUE),
                            "Normal", "Tumor"))

table(cases_filtered$condition)  # للتأكد

# 5️⃣ إزالة أي Metastatic samples (لو موجودة في definition)
cases_filtered <- cases_filtered %>%
  filter(!grepl("Metastatic", definition, ignore.case = TRUE))

# 6️⃣ التأكد من الأعمدة النهائية
head(cases_filtered)
colnames(cases_filtered)

# 7️⃣ حفظ النسخة الجاهزة
write.csv(cases_filtered, "/home/norah/Downloads/new task/cases_filtered_ready.csv", row.names = FALSE)

##############################
# النتيجة: cases_filtered_ready.csv جاهز للتحليل DEA
##############################

library(DESeq2)
rownames(cases_filtered) <- cases_filtered$barcode

# استبدال جميع النقاط (.) بالواصلات (-) في أسماء أعمدة counts_data
colnames(counts_data) <- gsub("\\.", "-", colnames(counts_data))

# قم بالتأكد من أن الأسماء تبدو صحيحة الآن
head(colnames(counts_data), 5)


# 1. تصفية counts_data
counts_data_filtered <- counts_data[, common_samples]

# 2. تصفية وترتيب cases_filtered (ترتيبها ليتطابق مع ترتيب أعمدة counts_data)
cases_filtered_final <- cases_filtered[common_samples, ]

# 3. التحقق الأخير
all.equal(colnames(counts_data_filtered), rownames(cases_filtered_final)) 
# يجب أن تكون TRUE

# 4. بناء كائن DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_data_filtered, 
  colData = cases_filtered_final, 
  design = ~ condition
)
rownames(cases_filtered) <- cases_filtered$barcode
# فلترة الجينات منخفضة العد
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
dim(dds)


vsd <- vst(dds, blind = TRUE)
library(ggplot2)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"%")) +
  ylab(paste0("PC2: ",percentVar[2],"%")) +
  ggtitle("PCA (VST) - قبل تصحيح batch") +
  theme_minimal()


pca_df <- pcaData
# حساب z-score لكل PC
pca_df$z1 <- scale(pca_df$PC1)
pca_df$z2 <- scale(pca_df$PC2)
# تعريف outlier: |z| > 3 في أي PC
outliers <- pca_df[abs(pca_df$z1) > 3 | abs(pca_df$z2) > 3, ]
outliers

# افتراض اسم العمود بالـ clinical: origin_of_tissue أو origin.of.tissue أو origin
# تأكدي من اسمه:
colnames(cases_filtered)


# عرض توزيع origin_of_tissue مع الأخذ في الاعتبار أي NA
table(cases_filtered$tissue_or_organ_of_origin, useNA = "ifany")


library(sva)

# نصيحة: استخدمي counts مقيّمة (log) أو VST لوضعية sva
mat <- assay(vsd)  # log-like transformed

# تصميم المصفوفات للموديل (condition) وموديل الـ null
mod <- model.matrix(~ condition, data = cases_filtered)
mod0 <- model.matrix(~ 1, data = cases_filtered)

# تشغيل svaseq لاستخراج SVs
svobj <- svaseq(as.matrix(mat), mod, mod0)
svs <- svobj$sv   # matrix of surrogate variables
dim(svs)
# أضيفيها للـ colData
for(i in seq_len(ncol(svs))){
  cases_filtered[[paste0("SV", i)]] <- svs[, i]
}
# اعادة بناء dds مع design يتضمن SVs و origin_of_tissue كـ batch (أو أحدهما)
# مثال design مع origin_of_tissue + أول 2 SV:
colData_new <- cases_filtered
dds <- DESeqDataSetFromMatrix(
  countData = counts_data_filtered,
  colData = colData_new,
  design = ~ tissue_or_organ_of_origin + SV1 + SV2 + condition
)

# تأكدي أنtissue_or_organ_of_origin عامل عاملية
dds$tissue_or_organ_of_origin <- factor(dds$tissue_or_organ_of_origin)


dds <- DESeq(dds, test = "Wald")
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
# ترتيب النتائج بحسب adjusted p-value
resOrdered <- res[order(res$padj), ]
summary(resOrdered)
# حفظ النتائج
write.csv(as.data.frame(resOrdered), file = "/home/norah/Downloads/new task/DE_results_all_genes.csv")



library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# افتراض أن rownames(resOrdered) هي Ensembl IDs مثل ENSG00000... 
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(resOrdered),
                  mart = ensembl)

# دمج الخريطة مع النتائج
res_df <- as.data.frame(resOrdered)
res_df$ensembl_gene_id <- rownames(res_df)
res_merged <- merge(res_df, gene_map, by = "ensembl_gene_id", all.x = TRUE)
# افتح الناتج وحفظ
write.csv(res_merged, "/home/norah/Downloads/new task/DE_results_with_symbols.csv", row.names = FALSE)


genes_of_interest <- c("ACADL","ACOX1","ADIPOQ","ADRB3","AGRP","CHKB","CIDEA","COX7A1",
                       "CPT1A","CPT1B","CPT1C","CPT2","DIO2","ELOVL3","FABP3","FABP4",
                       "FGF21","FNDC5","GHRL","GHSR","LEP","LEPR","LPL","MBOAT4","MC3R",
                       "MC4R","NPY","PLIN1","PM20D1","PPARA","PPARD","PPARG","PPARGC1A",
                       "PPARGC1B","PRDM16","PDK4","PNPLA2","RETN","SLC25A14","SLC25A20",
                       "SLC25A27","SLC27A1","SLC2A4","SREBF1","TBX1","TMEM26","UCP1","UCP2","UCP3")

# إذا res_merged يحتوي عمود hgnc_symbol:
goi_res <- res_merged[res_merged$hgnc_symbol %in% genes_of_interest, ]
# أظهر النتائج
goi_res
# حفظ GOI
write.csv(goi_res, "/home/norah/Downloads/new task/DE_results_GOI.csv", row.names = FALSE)


goi_res$significant <- with(goi_res, ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
                                            ifelse(log2FoldChange > 0, "Up", "Down"), "NotSig"))
table(goi_res$significant) 

write.csv(goi_res,
          "/home/norah/Downloads/new task/GOI_results_with_significance.csv",
          row.names = FALSE)


###############################
# 1) VST قبل إزالة الباتش
###############################

library(DESeq2)
library(sva)
library(ggplot2)
library(limma)

# VST pre-correction
vsd_pre <- vst(dds, blind = TRUE)
vsd_pre_mat <- assay(vsd_pre)

################################
# 2) PCA قبل إزالة الباتش (للكشف)
################################

pca_pre <- prcomp(t(vsd_pre_mat))

pca_pre_df <- data.frame(
  PC1 = pca_pre$x[,1],
  PC2 = pca_pre$x[,2],
  condition = colData(dds)$condition,
  tissue = colData(dds)$tissue_or_organ_of_origin
)

ggplot(pca_pre_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw(base_size = 14) +
  ggtitle("PCA Before Batch Correction")

########################################
# 3) SVA – تحديد المتغيرات المختلطة SVs
########################################

# النموذج الأساسي بدون batch
mod <- model.matrix(~ condition, data = colData(dds))

# نموذج الـnull
mod0 <- model.matrix(~ 1, data = colData(dds))

# حساب surrogate variables (SV1, SV2 ...)
svobj <- sva(vsd_pre_mat, mod, mod0)

# إضافة SVs للـcolData
for(i in 1:ncol(svobj$sv)) {
  colData(dds)[[paste0("SV", i)]] <- svobj$sv[, i]
}

########################################
# 4) إزالة الباتش من بيانات VST
########################################

vsd_corrected_mat <- removeBatchEffect(
  x = vsd_pre_mat,
  covariates = svobj$sv,
  design = mod
)

# تحويل المصفوفة مرة أخرى إلى DESeqTransform object
vsd_corrected <- vsd_pre
assay(vsd_corrected) <- vsd_corrected_mat

########################################
# 5) PCA بعد إزالة الباتش (المطلوبة)
########################################

pca_post <- prcomp(t(vsd_corrected_mat))

pca_post_df <- data.frame(
  PC1 = pca_post$x[,1],
  PC2 = pca_post$x[,2],
  condition = colData(dds)$condition,
  tissue = colData(dds)$tissue_or_organ_of_origin
)

ggplot(pca_post_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw(base_size = 14) +
  ggtitle("PCA After Batch Correction (SVA)")

###############################
# نهاية السكريبت
###############################



