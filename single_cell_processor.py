# single_cell_processor.py
import scanpy as sc

class SingleCellProcessor:
    def __init__(self, file_path):
        print("SingleCellProcessor başlatılıyor...", flush=True)
        self.adata = self.load_data(file_path)  # veri dosyasının yolunu alır ve yükler

    def load_data(self, file_path):
        print(f"Veri yükleniyor: {file_path}", flush=True)
        adata = sc.read(file_path)
        print("Veri başarıyla yüklendi!", flush=True)

        print("Veri boyutu:", adata.shape)
        print("Veri başı (obs):", adata.obs.head())
        print("Veri başı (var):", adata.var.head())

        if adata.var.empty:
            print("Gen bilgileri (var) bulunamadı, alternatif alanlar kontrol ediliyor...", flush=True)
            print("Veri başı (uns):", adata.uns.keys())
            print("Veri başı (obsm):", adata.obsm.keys())
            print("Veri başı (layers):", adata.layers.keys())
        else:
            print("Gen bilgileri mevcut.")
        
        return adata

    def quality_control(self):
        try:
            print("Kalite kontrol işlemi başlıyor...", flush=True)
            sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)
            print(f"Veri boyutu kalite kontrol sonrası: {self.adata.shape}", flush=True)
            print(f"Veri başı (obs) kalite kontrol sonrası:\n{self.adata.obs.head()}", flush=True)
            print(f"Veri başı (var) kalite kontrol sonrası:\n{self.adata.var.head()}", flush=True)
            print("Kalite kontrol tamamlandı!", flush=True)
        except Exception as e:
            print(f"Kalite kontrol hatası: {e}", flush=True)

    def filter_data(self):
        try:
            print("Filtreleme işlemi başlıyor...", flush=True)
            sc.pp.filter_cells(self.adata, min_genes=200)
            sc.pp.filter_cells(self.adata, max_genes=8000)
            sc.pp.filter_genes(self.adata, min_cells=3)
            print("Filtreleme tamamlandı!", flush=True)
        except Exception as e:
            print(f"Filtreleme hatası: {e}", flush=True)

    def normalize_and_scale(self):
        print("Normalize ve ölçekleme işlemi başlıyor...", flush=True)
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)  # yüksek varyanslı genleri seçer 
        sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        self.adata = self.adata[:, self.adata.var.highly_variable]
        sc.pp.scale(self.adata, max_value=10)
        print("Normalize ve ölçekleme tamamlandı!", flush=True)

    def compute_neighbors_and_umap(self):
        print("UMAP hesaplama işlemi başlıyor...", flush=True)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(self.adata)
        print("UMAP hesaplama tamamlandı!", flush=True)

    def plot_umap(self):
        print("UMAP çizimi başlıyor...", flush=True)
        return sc.pl.umap(self.adata, color=['Cell type'], show=False)
