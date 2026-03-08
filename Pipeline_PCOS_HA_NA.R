# SCRIPT: Transcriptomic Analysis of PCOS Subtypes (HA vs NA)
# DATASET: GSE137684 (Gene Expression Omnibus)
# DESCRIPTION: Pipeline for Microarray Data Fetching, Pre-processing, 
#              Differential Expression (limma), and Enrichment Analysis

# PART A. PERSIAPAN LIBRARY
# Memuat paket yang dibutuhkan untuk analisis dan visualisasi
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(umap)

# Library tambahan untuk Enrichment Analysis
# Jika belum terinstal, jalankan: BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# PART B. PENGAMBILAN DATA DARI GEO
message("Mengunduh dataset GSE137684...")
gset <- getGEO("GSE137684", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# PART C. PRE-PROCESSING DATA EKSPRESI
ex <- exprs(gset)

# Evaluasi kebutuhan Log2 Transformasi berdasarkan kuartil data
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  message("Log2 Transformasi telah dilakukan pada matriks ekspresi.")
} else {
  message("Data sudah dalam skala Log2.")
}

# PART D. DEFINISI KELOMPOK SAMPEL (DESIGN EXPERIMENT)
# Mengekstrak metadata sampel
sample_info <- pData(gset)$title
groups <- rep(NA, length(sample_info))

# Mendefinisikan grup berdasarkan pola teks pada judul sampel
groups[grep("Normoandrogenic", sample_info, ignore.case = TRUE)] <- "NA_PCOS"
groups[grep("Hyperandrogenic", sample_info, ignore.case = TRUE)] <- "HA_PCOS"
groups[grep("Control", sample_info, ignore.case = TRUE)] <- "Control"

# Menyaring data hanya untuk grup HA dan NA (Eksklusi Kontrol sehat)
keep_samples <- groups %in% c("HA_PCOS", "NA_PCOS")
ex <- ex[, keep_samples]
groups <- groups[keep_samples]

# Konversi grup menjadi faktor
gset_group <- factor(groups)
message("Kelompok yang dianalisis: ", paste(levels(gset_group), collapse = " vs "))

# PART E. DESIGN MATRIX & DIFFERENTIAL EXPRESSION (LIMMA)
# Membangun Design Matrix
design <- model.matrix(~0 + gset_group)
colnames(design) <- levels(gset_group)

# Menentukan formula kontras
contrast_formula <- "HA_PCOS - NA_PCOS"
message("Menjalankan analisis untuk kontras: ", contrast_formula)

# Menjalankan model limma
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Mengekstrak hasil (Threshold: p-value < 0.05)
topTableResults <- topTable(fit2, adjust = "none", p.value = 0.05, number = Inf, sort.by = "P")
message("Jumlah gen signifikan (P < 0.05): ", nrow(topTableResults))

# PART F. ANOTASI GENE SYMBOL
# Mengambil tabel anotasi dari metadata GEO
fData_table <- fData(gset)
probe_ids <- rownames(topTableResults)

# Mencocokkan ID Probe dengan Gene Symbol
indices <- match(probe_ids, rownames(fData_table))
topTableResults$SYMBOL <- fData_table[indices, "GENE_SYMBOL"]

# Data Cleaning: Menangani missing values pada simbol gen
# Jika SYMBOL kosong, fallback menggunakan Probe ID agar tidak terjadi error plot
topTableResults$plot_label <- ifelse(
  is.na(topTableResults$SYMBOL) | topTableResults$SYMBOL == "" | topTableResults$SYMBOL == "NA", 
  probe_ids, 
  topTableResults$SYMBOL
)

# PART G. VISUALISASI KUALITAS DATA (BOXPLOT)
# Mengatur margin bawah agar label tidak terpotong
par(mar = c(10, 4, 4, 2) + 0.1)

# Menetapkan warna berdasarkan grup
custom_colors <- ifelse(gset_group == "HA_PCOS", "salmon", "skyblue")

# Eksekusi Boxplot Normalisasi
boxplot(ex, 
        col = custom_colors, 
        las = 2, 
        main = "Distribusi Ekspresi Normalisasi: HA vs NA")

legend("topright", 
       legend = c("HA_PCOS", "NA_PCOS"), 
       fill = c("salmon", "skyblue"))

# PART H. VISUALISASI DEGs (VOLCANO PLOT)
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  P.Val = topTableResults$P.Value,
  Gene = topTableResults$plot_label
)

# Mengatur threshold: P < 0.05 & |logFC| > 0.5
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0.5 & volcano_data$P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -0.5 & volcano_data$P.Val < 0.05] <- "DOWN"

# Eksekusi Volcano Plot
p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot: DEGs HA vs NA")

print(p_volcano)

# PART I. VISUALISASI HEATMAP (TOP 50 DEGs)
# Ekstraksi 50 gen paling signifikan
top50 <- head(topTableResults[order(topTableResults$P.Value), ], 50)
mat_heatmap <- ex[rownames(top50), ]
rownames(mat_heatmap) <- top50$plot_label

# Label grup pasien
annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(mat_heatmap)

# Eksekusi Heatmap
pheatmap(mat_heatmap,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 50 DEGs (HA vs NA)")

# PART J. ENRICHMENT ANALYSIS (GO & KEGG) - GLOBAL LEVEL
# Menyiapkan daftar gen UP & DOWN regulasi untuk analisis Global
sig_genes <- volcano_data$Gene[volcano_data$status %in% c("UP", "DOWN")]

# Konversi Symbol ke Entrez ID
gene_ids <- bitr(sig_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# A. Gene Ontology (Biological Process) - Analisis Global
go_result <- enrichGO(gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",   
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)

if (!is.null(go_result)) {
  plot_go <- dotplot(go_result, showCategory = 10, title = "Top 10 Gene Ontology (BP)")
  print(plot_go)
} else {
  message("Tidak ada istilah GO yang signifikan.")
}

# B. KEGG Pathway Analysis - Analisis Global
# CATATAN PENELITI: Pada studi dengan ukuran sampel kecil (n=4 per grup), 
# penyesuaian statistik (pAdjust) seringkali membuang semua jalur. 
# Jika hasil di bawah ini kosong, gunakan "Targeted Pathway Mapping" 
# menggunakan web KEGG Mapper dengan data ekspor dari Part K.
kegg_result <- enrichKEGG(gene         = gene_ids$ENTREZID,
                          organism     = 'hsa', 
                          pAdjustMethod= "BH",
                          pvalueCutoff = 0.05)

if (!is.null(kegg_result)) {
  plot_kegg <- dotplot(kegg_result, showCategory = 10, title = "Top 10 KEGG Pathways")
  print(plot_kegg)
} else {
  message("INFO: Tidak ada KEGG Pathway global yang signifikan (p.adj > 0.05). Lanjutkan ke validasi Web KEGG Mapper.")
}

# PART K. EKSPOR DATA & PERSIAPAN VALIDASI WEB (g:Profiler & KEGG Mapper)
# 1. Menyimpan seluruh hasil DEGs lengkap
write.csv(topTableResults, "Hasil_GSE137684_HA_vs_NA.csv", row.names = FALSE)

# 2. Menyiapkan data untuk Validasi Web (Top 100 UP & DOWN)
# Memfilter hanya gen yang signifikan (P.Value < 0.05) dan memiliki nama Gen
sig_results <- topTableResults[topTableResults$P.Value < 0.05, ]
sig_results <- sig_results[!is.na(sig_results$plot_label) & !grepl("A_", sig_results$plot_label), ]

# Mengekstrak Top 100 Gen Up-regulated
up_genes <- sig_results[sig_results$logFC > 0, ]
up_genes <- up_genes[order(up_genes$logFC, decreasing = TRUE), ]
top100_up <- head(up_genes$plot_label, 100)

# Mengekstrak Top 100 Gen Down-regulated
down_genes <- sig_results[sig_results$logFC < 0, ]
down_genes <- down_genes[order(down_genes$logFC, decreasing = FALSE), ]
top100_down <- head(down_genes$plot_label, 100)

# Menyimpan sebagai file teks (.txt) agar siap Copy-Paste ke Web
write.table(top100_up, "Top100_UP_Genes_WebInput.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(top100_down, "Top100_DOWN_Genes_WebInput.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Pipeline Selesai! Data utama disimpan sebagai Hasil_GSE137684_HA_vs_NA.csv")
message("File input untuk g:Profiler dan KEGG Mapper telah di-generate (.txt).")