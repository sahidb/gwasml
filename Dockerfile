FROM continuumio/miniconda3:24.1.2-0

WORKDIR /app

# Salin file environment dan requirement
COPY environment.yml .
RUN conda env create -f environment.yml

# Aktifkan env default
SHELL ["conda", "run", "-n", "gwasml", "/bin/bash", "-c"]

# Salin kode pipeline dan folder data
COPY . .

# Pastikan data/image dan data/ exist
RUN mkdir -p data/image

# Set default entrypoint (gunakan argumen untuk override input/output filenya)
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "gwasml", "python", "gwas_pipeline.py"]
