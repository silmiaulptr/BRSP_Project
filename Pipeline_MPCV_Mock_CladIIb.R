# SCRIPT: Transcriptomic Analysis of Mpox (Mock vs Infected Clade IIb)
# DATASET: GSE219036 (RNA-seq / TPM Data)
# DESCRIPTION: Pipeline for Local Data Import, Pre-processing, 
#              Differential Expression (limma), and Enrichment Analysis

# PART A. PERSIAPAN LIBRARY
library(limma)
library(ggplot2)
library(pheatmap)
library(umap)

# Library tambahan untuk Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# PART B. PENGAMBILAN DATA (DARI FILE LOKAL)
file_path <- "C:/Users/SILMI AULIA PUTRI/Downloads/BRSP/GSE219036_TPM_for_GEO_upload_skin.txt.gz"
raw_data <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
message("Data lokal berhasil dimuat!")

# PART C. PRE-PROCESSING DATA EKSPRESI
target_columns <- c(2, 3, 4, 8, 9, 10)
ex <- raw_data[, target_columns]
ex <- log2(ex + 1)
message("Transformasi log2(x+1) telah dilakukan.")

# PART D. DEFINISI KELOMPOK SAMPEL
group_list <- c("Mock", "Mock", "Mock", "Infected", "Infected", "Infected")
group_factor <- factor(group_list, levels = c("Mock", "Infected"))
message("Grup yang dianalisis: ", paste(levels(group_factor), collapse = " vs "))

# PART E. DESIGN MATRIX & DIFFERENTIAL EXPRESSION (LIMMA)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
contrast_formula <- "Infected - Mock"
message("Menjalankan analisis untuk kontras: ", contrast_formula)

fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

# PART F. ANOTASI GENE SYMBOL (Konversi Ensembl ID ke Symbol)
# Memastikan library database manusia sudah aktif
library(org.Hs.eg.db)

# 1. Membersihkan Ensembl ID dari nomor versi 
# (Memotong karakter apa pun setelah tanda titik, misal ENSG...16 menjadi ENSG...)
clean_ensembl <- gsub("\\..*", "", rownames(topTableResults))

# 2. Melakukan pemetaan (Mapping) dari kamus ENSEMBL ke SYMBOL
mapped_symbols <- mapIds(org.Hs.eg.db,
                         keys = clean_ensembl,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# 3. Memasukkan hasil terjemahan ke dalam kolom SYMBOL
# Jika ada gen yang tidak memiliki nama (NA), kita tetap gunakan ID aslinya agar data tidak hilang
topTableResults$SYMBOL <- ifelse(is.na(mapped_symbols), rownames(topTableResults), mapped_symbols)

message("Konversi Ensembl ID ke Gene Symbol berhasil dilakukan!")

# PART G. VISUALISASI KUALITAS DATA (STANDAR SCOPUS)
# 1. Boxplot Distribusi
# Mengatur margin bawah agar label tidak terpotong
par(mar = c(10, 4, 4, 2) + 0.1)

# Menggunakan palet ramah buta warna (Okabe-Ito)
group_colors <- ifelse(group_factor == "Mock", "#0072B2", "#D55E00")
boxplot(ex, 
        col = group_colors, 
        las = 2, outline = FALSE,
        main = "Boxplot Distribusi Ekspresi Log2(TPM+1)", 
        ylab = "Expression Value")
legend("topright", legend = levels(group_factor), fill = c("#0072B2", "#D55E00"), cex = 0.8)

# 2. Density Plot (6 Garis, 6 Warna Gradasi Publikasi)
expr_long <- data.frame(
  Expression = as.vector(as.matrix(ex)),
  Sample = rep(colnames(ex), each = nrow(ex)),
  Group = rep(group_factor, each = nrow(ex))
)

# Membuat 6 warna khusus: 3 nuansa Biru (Mock) dan 3 nuansa Oranye/Merah (Infected)
# Warna ini sangat direkomendasikan untuk standar publikasi dan aman bagi buta warna
gradasi_warna <- c(
  "#56B4E9", "#0072B2", "#003366", # 3 Warna untuk Mock
  "#E69F00", "#D55E00", "#990000"  # 3 Warna untuk Infected Clade IIb
)

# Memetakan warna ke nama sampel secara otomatis
nama_sampel <- unique(expr_long$Sample)
warna_sampel <- setNames(gradasi_warna, nama_sampel)

p_density <- ggplot(expr_long, aes(x = Expression, color = Sample)) + 
  geom_density(linewidth = 1.2, alpha = 0.8) + 
  scale_color_manual(values = warna_sampel) +
  theme_minimal(base_size = 14) + 
  labs(title = "Distribusi Ekspresi Gen per Sampel (MPXV)", 
       x = "Expression Value (log2 TPM+1)", 
       y = "Density") +
  theme(legend.position = "right")

print(p_density)

# 3. PCA Plot (Sebagai Pengganti UMAP untuk Sampel Kecil)
# Menghitung varians untuk setiap gen (baris)
gene_vars <- apply(ex, 1, var)

# Membuang gen yang nilainya konstan (varians = 0) di semua sampel
ex_filtered <- ex[gene_vars > 0, ]

# Menghitung PCA pada data yang sudah bersih dari varians nol
pca_result <- prcomp(t(ex_filtered), scale. = TRUE)

# Mengekstrak 2 komponen utama (PC1 dan PC2) untuk plot
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group_factor
)

# Visualisasi PCA Standar Scopus
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) + 
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(values = c("Mock" = "#0072B2", "Infected" = "#D55E00")) +
  theme_minimal(base_size = 14) + 
  labs(title = "PCA Plot: Mock vs Infected Clade IIb",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme(legend.position = "right")

print(p_pca)

# PART H. VISUALISASI DEGs (VOLCANO & HEATMAP STANDAR SCOPUS)
# 1. Volcano Plot
volcano_data <- data.frame(
  logFC = topTableResults$logFC, 
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.7, size = 1.5) +
  # Skema warna elegan: Merah bata (UP), Biru laut (DOWN), Abu-abu redup (NO)
  scale_color_manual(values = c("DOWN" = "#0072B2", "NO" = "#CCCCCC", "UP" = "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") +
  ggtitle("Volcano Plot: Infected vs Mock")
print(p_volcano)

# i. STRATEGI PELABELAN (Anotasi Top 5 Gen UP dan Top 5 Gen DOWN)
volcano_data$label_teks <- NA

# Mengambil HANYA 5 gen UP dan 5 gen DOWN dengan p-value paling ekstrem
top5_up <- head(volcano_data[volcano_data$status == "UP", ][order(volcano_data[volcano_data$status == "UP", ]$adj.P.Val), "Gene"], 5)
top5_down <- head(volcano_data[volcano_data$status == "DOWN", ][order(volcano_data[volcano_data$status == "DOWN", ]$adj.P.Val), "Gene"], 5)

# Fungsi cerdas untuk memotong teks panjang menjadi maksimal 15 karakter
potong_teks <- function(x) {
  ifelse(nchar(x) > 15, paste0(substr(x, 1, 12), "..."), x)
}

# Memasukkan nama gen yang SUDAH DIPENDEKAN ke kolom label teks
volcano_data$label_teks <- ifelse(
  volcano_data$Gene %in% c(top5_up, top5_down), 
  potong_teks(as.character(volcano_data$Gene)), 
  NA
)

# ii. Eksekusi Plot dengan geom_label_repel yang Lebih Rapi
p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("DOWN" = "#0072B2", "NO" = "#E0E0E0", "UP" = "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  geom_label_repel(aes(label = label_teks), 
                   size = 3.5, 
                   fontface = "bold", 
                   color = "black",    
                   fill = "white",     
                   box.padding = 1.0,  # Jarak antar kotak diperlebar
                   point.padding = 0.5, 
                   segment.color = "black", 
                   segment.size = 0.5,
                   max.overlaps = Inf,
                   force = 10,         # Kekuatan menolak antar kotak ditingkatkan (agar tidak numpuk)
                   na.rm = TRUE,       
                   show.legend = FALSE) + 
  
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") +
  ggtitle("Volcano Plot: Infected vs Mock (Top 10 DEGs Labeled)")

print(p_volcano)

# 2. Heatmap (Top 50 DEGs)
top50_genes <- head(topTableResults[order(topTableResults$adj.P.Val), ], 50)
mat_heatmap <- ex[rownames(top50_genes), ]

# MENGGANTI NAMA BARIS HEATMAP MENJADI GENE SYMBOL
rownames(mat_heatmap) <- top50_genes$SYMBOL

annotation_col <- data.frame(Group = group_factor)
rownames(annotation_col) <- colnames(mat_heatmap)

ann_colors <- list(Group = c(Mock = "#0072B2", Infected = "#D55E00"))
heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

pheatmap(mat_heatmap, 
         scale = "row", 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = heat_colors,
         show_colnames = FALSE, 
         show_rownames = TRUE,  # Sekarang teks di kanan akan berupa Gene Symbol!
         fontsize_row = 7,
         clustering_method = "complete", 
         main = "Top 50 DEGs: Infected vs Mock")


# PART I. ENRICHMENT ANALYSIS (GO & KEGG) - GLOBAL LEVEL
sig_genes <- volcano_data$Gene[volcano_data$status %in% c("UP", "DOWN")]
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_result <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
if (!is.null(go_result)) {
  print(dotplot(go_result, showCategory = 10, title = "Top 10 Gene Ontology (BP)"))
} else { message("Tidak ada istilah GO yang signifikan.") }

kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', 
                          pAdjustMethod= "BH", pvalueCutoff = 0.05)
if (!is.null(kegg_result)) {
  print(dotplot(kegg_result, showCategory = 10, title = "Top 10 KEGG Pathways"))
} else { message("INFO: Jika KEGG Global kosong. Lanjutkan ke validasi web.") }

# PART J. EKSPOR DATA & PERSIAPAN VALIDASI WEB
write.csv(topTableResults, "Hasil_DEGs_Mpox_Infected_vs_Mock.csv", row.names = FALSE)

sig_results <- topTableResults[topTableResults$adj.P.Val < 0.05 & abs(topTableResults$logFC) > 1, ]
sig_results <- sig_results[!is.na(sig_results$SYMBOL) & sig_results$SYMBOL != "", ]

up_genes <- sig_results[sig_results$logFC > 0, ]
up_genes <- up_genes[order(up_genes$logFC, decreasing = TRUE), ]
top100_up <- head(up_genes$SYMBOL, 100)

down_genes <- sig_results[sig_results$logFC < 0, ]
down_genes <- down_genes[order(down_genes$logFC, decreasing = FALSE), ]
top100_down <- head(down_genes$SYMBOL, 100)

write.table(top100_up, "Top100_UP_Mpox.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(top100_down, "Top100_DOWN_Mpox.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

message(paste("Pipeline Selesai! File tersimpan di direktori:", getwd()))