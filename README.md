# 🧬 Portofolio Bioinformatika: Analisis Transkriptomik & Differentially Expressed Genes (DEGs)

[![OmicsLite BRSP](https://img.shields.io/badge/Program-BRSP_OmicsLite-purple.svg)](#)
[![R Programming](https://img.shields.io/badge/Language-R-blue.svg)](#)
[![Status](https://img.shields.io/badge/Status-Completed-success.svg)](#)

## Tentang Repositori Ini
Repositori ini berisi kumpulan tugas, mini-proyek, dan proyek akhir (*Capstone Project*) yang saya kerjakan selama mengikuti **Bioinformatics Research Starter Program (BRSP)** dari OmicsLite. Program pelatihan 5 minggu ini menggunakan pendekatan *Project-Based Learning* (PBL) untuk mempelajari analisis data *transcriptomics* (khususnya *Differentially Expressed Gene* / DEG) secara runtut dan aplikatif. Dalam repositori ini, saya mendemonstrasikan alur riset bioinformatika secara utuh, mulai dari literasi ilmiah, eksplorasi basis data Gene Expression Omnibus (GEO), pemrosesan data dengan R, hingga interpretasi biologis (*Enrichment Analysis*).

---

## 🏆 *Capstone Project* (Proyek Akhir)
Sebagai penugasan akhir program ini, saya melakukan analisis komparatif ekspresi gen secara mandiri. 

* **Topik Riset:** Analisis Transkriptomik MPXV (Mpox Virus) Clade IIb (Mock vs MPXV).
* **Dataset yang Digunakan:** `GSE219036` (Data TPM RNA-seq sampel kulit).
* **Hasil & Laporan Utama:** Silakan baca detail analisis, visualisasi, dan interpretasi biologisnya di dokumen **[Analisis Transkriptomik MPXV CladIIb](Analisis_Transkriptomik_MPXV_CladIIb.md)**.
* **Skrip R:** [`Pipeline_MPCV_Mock_CladIIb.R`](Pipeline_MPCV_Mock_CladIIb.R)
* **Dataset Mentah:** [`Dataset_GSE219036_TPM_for_GEO_upload_skin.txt.gz`](Dataset_GSE219036_TPM_for_GEO_upload_skin.txt.gz)

---

## Rekam Jejak Pembelajaran (Minggu 1 - 3)
Selain proyek akhir, repositori ini juga mendokumentasikan proses belajar saya secara bertahap selama program berlangsung:

### 🔹 Minggu 1: Konseptual dan Literasi Ilmiah
Mempelajari dasar *transcriptomics* dan berlatih melakukan telaah literatur ilmiah.
* 📄 **[Review Jurnal: Efek Glabridin pada Biofilm Staphylococcus aureus](Review-Jurnal-Efek-Glabridin-pada-Biofilm-Staphylococcus-aureus.md)**

### 🔹 Minggu 2: Pengenalan Data & Analisis Berbasis Web (GEO2R)
[cite_start]Eksplorasi data microarray dari NCBI GEO (Dataset `GSE137684` - Sel Granulosa PCOS) dan melakukan analisis DEG tahap awal.
* 📄 **[Analisis Patogenesis PCOS DEGs pada Sel Granulosa](Analisis%20Patogenesis%20PCOS%20DEGs%20pada%20Sel%20Granulosa.md)**

### 🔹 Minggu 3: *Workflow* DEG Analysis menggunakan R
Menerapkan analisis reproduktif (*reproducible analysis*) menggunakan bahasa pemrograman R untuk membandingkan pasien *Normoandrogenic* (NA) dan *Hyperandrogenic* (HA) PCOS.
* 📄 **Laporan Analisis:** **[Analisis Transkriptomik Normoandrogenic dan Hyperandrogenic PCOS](Analisis%20Transkiptomik%20Normoandrogenic%20dan%20Hyperandrogenic%20PCOS.md)**
* **Skrip Analisis:** [`Pipeline_PCOS_HA_NA.R`](Pipeline_PCOS_HA_NA.R)
* **Hasil Output:** [`Hasil DEGs_GSE137684_HA_vs_NA.csv`](Hasil%20DEGs_GSE137684_HA_vs_NA.csv)
* **Dataset Lokal:** [`Dataset_GSE137684_series_matrix.txt.gz`](Dataset_GSE137684_series_matrix.txt.gz)

---

## Tools & Packages yang Digunakan

Proyek ini menggunakan bahasa pemrograman **R / RStudio** dengan memanfaatkan berbagai *packages* dari ekosistem CRAN dan Bioconductor:

* **Pengambilan & Manajemen Data:**
  * `GEOquery`: Untuk mengunduh dataset microarray (GSE137684) beserta tabel anotasinya secara langsung dari *database* NCBI GEO.
  * `dplyr`: Digunakan untuk manipulasi dan transformasi struktur data (data wrangling).
  * *Base R* (`read.table`): Digunakan untuk proses ekstraksi dan pembacaan *raw expression matrix* dari file lokal berekstensi `.txt.gz`.

* **Analisis Statistik & *Differential Expression* (DEG):**
  * `limma`: *Package* utama yang digunakan untuk pemodelan linear eksperimen (*design matrix*, *contrasts*), perhitungan *fold-change*, dan penapisan signifikansi statistik (p-value / FDR).

* **Anotasi Gen & *Enrichment Analysis*:**
  * `org.Hs.eg.db`: Basis data anotasi genom manusia yang digunakan untuk pemetaan ID (*mapping*), seperti mengonversi *Ensembl ID* atau *Probe ID* menjadi *Gene Symbol* dan *Entrez ID*.
  * `clusterProfiler`: Digunakan untuk mengeksekusi analisis pengayaan fungsi biologis (Gene Ontology - Biological Process) dan jalur metabolisme (KEGG Pathways).
  * `enrichplot`: Digunakan secara spesifik untuk memvisualisasikan hasil dari analisis *clusterProfiler* ke dalam bentuk *Dotplot*.

* **Visualisasi Data Eksekusi Tinggi:**
  * `ggplot2`: Digunakan secara luas untuk membuat visualisasi utama seperti *Volcano Plot*, *PCA Plot*, *Boxplot*, dan *Density Plot*.
  * `pheatmap`: Untuk membangun *Heatmap* hierarkis dari Top 50 DEGs, lengkap dengan anotasi warna pengelompokan sampel pasien.
  * `umap` / `prcomp` (PCA): Digunakan pada tahap analisis komponen utama (reduksi dimensi) untuk memvalidasi pemisahan pola varians antar grup sampel.
  * `ggrepel` (implisit pada `geom_label_repel`): Berperan dalam anotasi teks *Volcano Plot* agar label spesifik gen (Top DEGs) otomatis terhindar dari tumpang tindih (*overlapping*).

## 💡 Cara Mereproduksi Analisis
Jika Anda ingin mengeksekusi skrip R di komputer Anda:
1. *Clone* repositori ini menggunakan Git: 
   ```bash
   git clone [https://github.com/silmiaulptr/BRSP_Project.git](https://github.com/silmiaulptr/BRSP_Project.git)
