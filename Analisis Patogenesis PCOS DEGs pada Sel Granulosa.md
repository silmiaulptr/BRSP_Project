# Analisis Patogenesis PCOS Melalui Ekspresi Gen Diferensial (DEGs) pada Sel Granulosa 

### Latar Belakang
*Polycystic Ovary Syndrome* (PCOS) adalah gangguan endokrin dan metabolik yang menjadi penyebab utama infertilitas anovulatorik pada wanita usia reproduksi. Salah satu karakteristik utama PCOS adalah terhentinya perkembangan folikel ovarium (*follicular arrest*). Secara patofisiologis, sel granulosa sel somatik penyokong folikel memainkan peran krusial dalam patogenesis penyakit ini. Analisis komparatif ini bertujuan untuk mengevaluasi profil *Differentially Expressed Genes* (DEGs) pada sel granulosa pasien PCOS dibandingkan dengan wanita sehat (Kontrol), guna mengidentifikasi gen-gen penggerak (*driver genes*) yang merusak fungsi reproduksi.

---

### Metodologi Analisis

![Flowchart Metodologi](Flowchart_Metodologi_Granulosa.png)
*Gambar 1. Alur kerja analisis transkriptomik dari persiapan data hingga interpretasi hasil.*

Seperti yang diilustrasikan pada **Gambar 1**, tahapan analisis dirancang secara sistematis. Data matriks ekspresi diekstraksi dari *Gene Expression Omnibus* (GEO) publik dengan aksesi **GSE155489**. Analisis komparatif antara kelompok Kontrol (Normal) dan PCOS dieksekusi menggunakan perangkat lunak statistik berbasis web, GEO2R. Kontrol kualitas yang ketat diterapkan menggunakan metode koreksi *False Discovery Rate* (FDR) dengan ambang batas *Adjusted p-value* < 0.05 untuk menyaring gen yang valid. Output data DEG yang lolos seleksi kemudian diekstraksi untuk divisualisasikan.

---

### Hasil dan Interpretasi Biologis

#### 1. Evaluasi Kualitas Data (Boxplot)
Sebelum mengidentifikasi gen-gen spesifik, sangat penting untuk memastikan data terbebas dari bias teknis (*batch effect*).

![Boxplot Normalisasi](Boxplot_Granulosa.jpg)

*Gambar 2. Distribusi nilai median ekspresi antar-sampel yang telah dinormalisasi.*

Berdasarkan visualisasi *Boxplot* pada **Gambar 2**, nilai median intensitas ekspresi (Log2) pada seluruh sampel berada pada garis horizontal yang sejajar. Kesejajaran ekuilibrium ini mengonfirmasi bahwa data hitungan ekspresi dari repositori GEO telah ternormalisasi dengan sangat baik. Variasi yang dianalisis dipastikan murni akibat kondisi patologis PCOS, bukan akibat kesalahan instrumen.

#### 2. Pemisahan Identitas Molekuler (UMAP)
Untuk memvalidasi apakah penyakit ini secara fundamental merombak profil genetik, dilakukan analisis proyeksi dimensi reduksi spasial.

![UMAP Plot](UMAP_Granulosa.jpg)

*Gambar 3. UMAP Plot menunjukkan pemisahan koordinat spasial yang tegas berdasarkan kondisi penyakit.*

Berdasarkan hasil visualisasi pada **Gambar 3**, *UMAP plot* menampilkan klasterisasi spasial yang sangat tajam dan tanpa *overlap* antara kelompok sampel Kontrol dan PCOS. Jarak yang merentang jauh antar kedua klaster ini membuktikan bahwa sel granulosa penderita PCOS telah kehilangan identitas molekuler normalnya akibat paparan penyakit.

#### 3. Asimetri Transkripsional Global (Volcano Plot)
Setelah validasi jarak molekuler dipastikan, distribusi keseluruhan gen yang mengalami perubahan aktivitas divisualisasikan.

![Volcano Plot](Volcano_Granulosa.jpg)

*Gambar 4. Volcano Plot menunjukkan distribusi DEGs yang asimetris.*

Analisis diferensial pada **Gambar 4 (Volcano Plot)** mengungkap fenomena asimetri ekstrem. Lanskap grafik sangat didominasi oleh kepadatan titik pada area *LogFC* positif, yang merepresentasikan gen-gen *up-regulation* secara signifikan. Dominasi tak lazim ini mengindikasikan terjadinya **hiperaktivitas transkripsional** di dalam sel granulosa PCOS. Hal ini merupakan respons kompensasi seluler terhadap stres kronis di lingkungan ovarium, yang memaksa sel terus-menerus memproduksi mRNA secara berlebihan.

#### 4. Karakterisasi Molekuler Penggerak Utama (Top DEGs)
Untuk memahami kerusakan spesifik yang terjadi, data difokuskan pada gen-gen penyandi fungsional dengan nilai signifikansi tertinggi.

![Tabel Top DEGs](Tabel_DEGs_Granulosa.png)

*Gambar 5. Daftar gen teratas (Top DEGs) beserta metrik signifikansinya.*

Merujuk pada metrik di **Gambar 5**, disfungsi sel granulosa pada PCOS terbukti didalangi oleh malfungsi beberapa jalur mekanisme:
1. Gen **MMD** (*Monocyte to macrophage differentiation-associated*) menempati urutan signifikansi teratas secara absolut (*padj* = 1.11E-69, Log2FC = 2.264), diikuti oleh peningkatan gen **NFKBIA** (Log2FC = 1.643). Sinergi keduanya mengonfirmasi aktivasi makrofag folikuler agresif yang menciptakan badai sitokin di dalam ovarium yang merusak oosit.
2. Lonjakan ekspresi drastis terjadi pada gen **CLDN11** (Log2FC = 1.86) yang menyandikan protein *tight junction*. *Over expression* ini mengubah *Blood Follicle Barrier* menjadi kaku dan fibrotik, memblokir lalu lintas nutrisi dan membuat oosit kelaparan.
3. *Down regulation* pada gen **KCNK1** (Log2FC = -3.699) merusak fungsi kanal kalium, menghancurkan stabilitas "kelistrikan" seluler, dan membuat folikel menjadi kebal terhadap sinyal ovulasi dari tubuh.

---

### Kesimpulan
*Pipeline* komputasi transkriptomik yang diterapkan pada dataset GSE155489 secara komprehensif membuktikan bahwa PCOS memicu pergeseran molekuler yang masif, asimetris, dan destruktif pada sel granulosa. Disregulasi konstelasi gen kunci yakni hiperaktivasi **MMD** dan **NFKBIA** (inflamasi fokal), **CLDN11** (fibrosis barier folikel), serta represi **KCNK1** (kegagalan ion seluler) mentransformasi lingkungan ovarium sehat menjadi lingkungan patologis yang memicu berhentinya pematangan sel telur. 
