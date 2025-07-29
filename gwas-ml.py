import pandas as pd
import numpy as np
from sklearn.linear_model import Lasso, Ridge, ElasticNet
from sklearn.feature_selection import mutual_info_regression
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from xgboost import XGBRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, mean_absolute_error
import pybiomart
import os
import subprocess  # Untuk menjalankan Plink dari Python
import matplotlib.pyplot as plt  # Untuk plotting
from dotenv import load_dotenv  # Tambahan untuk memuat .env

# Muat .env jika ada
load_dotenv()

# Fungsi 1: Quality Control (QC) menggunakan Plink2
def perform_qc(input_file_prefix, output_file_prefix):
    """
    Jalankan QC menggunakan Plink2 untuk membersihkan data GWAS.
    Threshold diambil dari .env atau default.
    """
    geno = float(os.getenv('QC_GENO_THRESHOLD', 0))
    mind = float(os.getenv('QC_MIND_THRESHOLD', 0.05))
    maf = float(os.getenv('QC_MAF_THRESHOLD', 0.01))
    ld_window = int(os.getenv('QC_LD_WINDOW', 50))
    ld_step = int(os.getenv('QC_LD_STEP', 5))
    ld_r2 = float(os.getenv('QC_LD_R2_THRESHOLD', 0.3))
    hwe = float(os.getenv('QC_HWE_PVALUE', 1e-5))
    
    cmd = [
        'plink2',
        '--bfile', input_file_prefix,
        '--geno', str(geno),
        '--mind', str(mind),
        '--maf', str(maf),
        '--indep-pairwise', str(ld_window), str(ld_step), str(ld_r2),
        '--hwe', str(hwe),
        '--make-bed',
        '--out', output_file_prefix
    ]
    subprocess.run(cmd, check=True)
    print(f"QC selesai: Output di {output_file_prefix}")
    return output_file_prefix

# Fungsi 2: Feature Selection (Seleksi SNP)
def select_features(X, y, method='elastic_net'):
    """
    Pilih SNP relevan menggunakan metode ML.
    Parameter diambil dari .env atau default.
    """
    n_features = int(os.getenv('FEATURE_SELECTION_N_FEATURES', 5000))
    
    if method == 'lasso':
        alpha = float(os.getenv('FEATURE_SELECTION_LASSO_ALPHA', 0.00045))
        model = Lasso(alpha=alpha)
        model.fit(X, y)
        selected = X.columns[np.abs(model.coef_) > 0][:n_features]
    elif method == 'ridge':
        alpha = float(os.getenv('FEATURE_SELECTION_RIDGE_ALPHA', 0.005))
        model = Ridge(alpha=alpha)
        model.fit(X, y)
        coefs = np.abs(model.coef_)
        selected = X.columns[np.argsort(coefs)[-n_features:]]
    elif method == 'elastic_net':
        alpha = float(os.getenv('FEATURE_SELECTION_ELASTIC_NET_ALPHA', 0.0033))
        l1_ratio = float(os.getenv('FEATURE_SELECTION_ELASTIC_NET_L1_RATIO', 0.5))
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
        model.fit(X, y)
        selected = X.columns[np.abs(model.coef_) > 0][:n_features]
    elif method == 'mutual_info':
        mutual_top = int(os.getenv('FEATURE_SELECTION_MUTUAL_INFO_TOP', 5000))
        scores = mutual_info_regression(X, y)
        selected = X.columns[np.argsort(scores)[-mutual_top:]]
    else:
        raise ValueError("Metode tidak dikenal")
    
    print(f"{len(selected)} SNP terpilih dengan {method}")
    return selected

# Fungsi Tambahan: Plot Performance
def plot_performance(r2, mae, model_name, save_dir="data/image"):
    """
    Buat dan simpan bar chart untuk metrik R² dan MAE.
    """
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    metrics = ['R2', 'MAE']
    scores = [r2, mae]
    colors = ['steelblue', 'darkorange']
    
    plt.figure(figsize=(6, 4))
    bars = plt.bar(metrics, scores, color=colors)
    plt.title(f"Performance Metrics for {model_name}")
    plt.ylim(0, max(scores) * 1.2)
    plt.ylabel('Score')
    
    for bar, score in zip(bars, scores):
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + max(scores)*0.02, f'{score:.2f}', ha='center', va='bottom')
    
    plt.savefig(os.path.join(save_dir, f"{model_name}_performance.png"))
    plt.close()
    print(f"Grafik performa disimpan: {save_dir}/{model_name}_performance.png")

# Fungsi 3: Association Analysis
def perform_association(X, y, method='svr', save_dir='data/image'):
    """
    Analisis asosiasi SNP dengan fenotipe menggunakan model ML.
    Parameter diambil dari .env atau default.
    """
    n_top = int(os.getenv('ASSOCIATION_TOP_SNPS', 100))
    
    if method == 'lr':
        model = LinearRegression()
    elif method == 'rf':
        n_estimators = int(os.getenv('ASSOCIATION_RF_N_ESTIMATORS', 100))
        max_depth = int(os.getenv('ASSOCIATION_RF_MAX_DEPTH', 20))
        model = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth)
    elif method == 'svr':
        c = float(os.getenv('ASSOCIATION_SVR_C', 100))
        epsilon = float(os.getenv('ASSOCIATION_SVR_EPSILON', 0.01))
        gamma = os.getenv('ASSOCIATION_SVR_GAMMA', 'scale')
        model = SVR(kernel='sigmoid', C=c, epsilon=epsilon, gamma=gamma)
    elif method == 'xgboost':
        learning_rate = float(os.getenv('ASSOCIATION_XGBOOST_LEARNING_RATE', 0.1))
        max_depth = int(os.getenv('ASSOCIATION_XGBOOST_MAX_DEPTH', 3))
        n_estimators = int(os.getenv('ASSOCIATION_XGBOOST_N_ESTIMATORS', 300))
        model = XGBRegressor(learning_rate=learning_rate, max_depth=max_depth, n_estimators=n_estimators)
    else:
        raise ValueError("Metode tidak dikenal")
    
    model.fit(X, y)
    preds = model.predict(X)
    r2 = r2_score(y, preds)
    mae = mean_absolute_error(y, preds)
    
    # Permutation importance untuk top SNP
    importance = permutation_importance(model, X, y, n_repeats=10)
    top_indices = np.argsort(importance.importances_mean)[-n_top:]
    top_snps = X.columns[top_indices]
    
    # Panggil fungsi plot
    plot_performance(r2, mae, method, save_dir)
    
    results = {'r2': r2, 'mae': mae, 'top_snps': top_snps}
    print(f"Metrik {method}: R²={r2:.2f}, MAE={mae:.2f}")
    return results

# Fungsi 4: SNP Enrichment (Analisis Pasca-GWAS)
def perform_enrichment(snps, output_csv='enrichment.csv'):
    """
    Validasi biologis SNP menggunakan pybiomart.
    """
    dataset = pybiomart.Dataset(name='hsapiens_snp_variation', host='http://www.ensembl.org')
    filters = {'snp_ids': snps}
    attributes = ['refsnp_id', 'chrom_start', 'associated_gene', 'phenotype_description']
    results = dataset.query(attributes=attributes, filters=filters)
    
    results.to_csv(output_csv, index=False)
    print(f"Enrichment selesai: Output di {output_csv}")
    return results

# Fungsi Tambahan: Untuk Imputed Data
def impute_data(input_file, output_file):
    """
    Imputasi data menggunakan tools eksternal.
    """
    cmd = ['plink2', '--bfile', input_file, '--make-bed', '--out', output_file]  # Placeholder
    subprocess.run(cmd, check=True)
    print(f"Imputasi selesai: {output_file}")
    return output_file

# Fungsi Tambahan: Untuk Rare Variants
def analyze_rare_variants(input_file):
    """
    Analisis varian langka.
    """
    maf_threshold = float(os.getenv('RARE_VARIANTS_MAF_THRESHOLD', 0.01))
    data = pd.read_csv(input_file + '.bim', sep='\t', header=None)
    rare = data[data[4] < maf_threshold]
    print(f"{len(rare)} varian langka ditemukan")
    return rare

# Fungsi Main untuk Jalankan Seluruh Pipeline
def main(input_prefix='data/input', phenotype_file='data/phenotype.csv'):
    qc_output = perform_qc(input_prefix, 'data/qc_output')
    
    # Placeholder load data
    X = pd.DataFrame(np.random.rand(100, 10000))
    y = pd.Series(np.random.rand(100))
    
    selected_snps = select_features(X, y, method='elastic_net')
    X_selected = X[selected_snps]
    
    methods = ['lr', 'rf', 'svr', 'xgboost']
    for method in methods:
        assoc_results = perform_association(X_selected, y, method=method)
    
    perform_enrichment(assoc_results['top_snps'])
    
    imputed = impute_data(qc_output, 'data/imputed')
    analyze_rare_variants(imputed)

if __name__ == "__main__":
    # Parsing manual, bisa diganti argparse jika mau lebih aman
    import sys
    input_prefix = sys.argv[1] if len(sys.argv) > 1 else 'data/input'
    phenotype_file = sys.argv[2] if len(sys.argv) > 2 else 'data/phenotype.csv'
    main(input_prefix, phenotype_file)
