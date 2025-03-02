# main.py
import logging
logging.basicConfig(level=logging.DEBUG)  # Hata ayıklama bilgilerini kaydetmek

import scanpy as sc
import anndata
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shiny import App, ui, render
from single_cell_processor import SingleCellProcessor  # Sınıfı import ettik

# Kullanıcı arayüzü UI kısmı: Dosya yükleme
app_ui = ui.page_fluid(
    ui.input_file("file1", "Veri dosyasını seçin"),  # Dosya yüklemek için
    ui.output_plot("umap_plot")  # UMAP grafiğini göstermek için
)

# Server kısmı yüklenen dosyayı işler ve umap oluşturur
def server(input, output, session):
    print("Server başlatıldı!", flush=True)

    @output
    @render.plot
    def umap_plot():
        file = input.file1()  # Dosyayı al
        if file is not None:  # Dosya yüklendiyse
            if isinstance(file, list):  # Eğer file bir liste ise
                print(f"Yüklenen dosyanın listesi: {file}", flush=True)
                file_path = file[0]['datapath']  # İlk dosyanın yolunu al
            else:
                file_path = file['datapath']  # Eğer dict türünde ise
            print(f"Yüklenen dosya yolu: {file_path}", flush=True)  # Debug: Dosya yolunu yazdır
            processor = SingleCellProcessor(file_path)
            processor.quality_control()
            processor.filter_data()
            processor.normalize_and_scale()
            processor.compute_neighbors_and_umap()
            return processor.plot_umap()
        else:
            fig, ax = plt.subplots()
            ax.set_title("Lütfen bir dosya yükleyin.") 
            return fig

app = App(app_ui, server)  # Uygulamayı başlat

# Sunucu başlatma komutu
app.run()
