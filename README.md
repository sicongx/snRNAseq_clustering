# snRNA-seq median+5SD QC, SCTransform / LogNormalize, clustree

两个 snRNA-seq 项目的完整分析流程：**scn1lab** 与 **dyrk1a**。从 10x CellRanger 矩阵读入，经 QC、过滤、SCTransform（失败时自动 LogNormalize）、Harmony、clustree 与聚类，全部在 Rmd 中按步骤运行并重新生成结果。

本仓库**仅含代码与 Rmd**；原始 10x 矩阵体积过大（单个 `matrix.mtx.gz` 通常 >100 MB），需单独获取后放到 `data/` 再运行。

## 目录结构

```
snRNAseq_median5sd_clustering/
├── scripts/
├── scn1lab_snRNAseq/
│   ├── scn1lab_full_pipeline_median5sd_clustree.Rmd
│   ├── data/          # 自行放入 10x 矩阵（见下）
│   └── out/           # 运行后生成
└── dyrk1a_snRNAseq/
    ├── dyrk1a_full_pipeline_median5sd_clustree.Rmd
    ├── data/
    └── out/
```

路径均相对仓库根目录 `snRNAseq_median5sd_clustering/`。

## 数据准备

从 CellRanger 输出的 `filtered_feature_bc_matrix/`（或等价目录）复制三个文件到对应样本文件夹：

```
<project>/data/<folder>/
  barcodes.tsv.gz
  features.tsv.gz
  matrix.mtx.gz
```

**scn1lab**（`scn1lab_snRNAseq/data/`）：

| 文件夹 | 样本 ID |
|--------|---------|
| `scnWT1` | scn1_wt_rep1 |
| `scnWT2` | scn1_wt_rep2 |
| `scnWT3` | scn1_wt_rep3 |
| `scnHET1` | scn1_het_rep1 |
| `scnHET2` | scn1_het_rep2 |
| `scnHET3` | scn1_het_rep3 |
| `scnHOM1` | scn1_hom_rep1 |
| `scnHOM2` | scn1_hom_rep2 |
| `scnHOM3` | scn1_hom_rep3 |

**dyrk1a**（`dyrk1a_snRNAseq/data/`）：

| 文件夹 | 样本 ID |
|--------|---------|
| `dyrk1a_WT1_DRT` | dyrk1a_wt_rep1 |
| `dyrk1a_WT2_DRT` | dyrk1a_wt_rep2 |
| `dyrk1a_WT3_DRT` | dyrk1a_wt_rep3 |
| `dyrk1a_HET1_DRT` | dyrk1a_het_rep1 |
| `dyrk1a_HET2_DRT` | dyrk1a_het_rep2 |
| `dyrk1a_HET3_DRT` | dyrk1a_het_rep3 |
| `dyrk1a_HOM1_DRT` | dyrk1a_hom_rep1 |
| `dyrk1a_HOM2_DRT` | dyrk1a_hom_rep2 |
| `dyrk1a_HOM3_DRT` | dyrk1a_hom_rep3 |

放好后，打开 Rmd 运行 `paths_check` chunk，应列出 9 个样本且 `matrix_dir` 均存在。

## 环境

- R ≥ 4.2，RStudio
- 建议内存 ≥ 64 GB（Step 6 合并 SCT 占用较大；内存不足时会自动改用 LogNormalize）

```r
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "tibble", "patchwork",
  "harmony", "future", "clustree", "rmarkdown"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scDblFinder", "glmGamPoi"),
                     update = FALSE, ask = FALSE)
```

## 如何运行

1. 完成上述数据放置。
2. 在 RStudio 中将工作目录设为 `snRNAseq_median5sd_clustering/`。
3. 打开 Rmd 并 **Run All**（或按 chunk 逐步运行）：
   - `scn1lab_snRNAseq/scn1lab_full_pipeline_median5sd_clustree.Rmd`
   - `dyrk1a_snRNAseq/dyrk1a_full_pipeline_median5sd_clustree.Rmd`
4. **Step 6 分三步**（`pipeline-options` 中开关）：
   - **6a** `run_step6a_sctransform <- TRUE`
   - **6b** 6a 完成后 `run_step6b_clustree <- TRUE`
   - **6c** 根据 clustree 修改 `clustering_resolutions`，`run_step6c_clustering <- TRUE`

固定 QC：`nFeature_min = 200`，`nFeature_max = 每样本 median + 5×SD`，`mt_max = 3`，仅 scDblFinder 单细胞。

## 产出结果

| 路径 | 步骤 | 内容 |
|------|------|------|
| `out/qc_plots/` | 4 | 合并前 QC 图 |
| `out/threshold_diagnostics/` | 4a–4b | 阈值诊断 |
| `out/qc_filter_plots/` | 5 | 过滤前后对比 |
| `out/*_merged_raw.rds` | 3 | 合并 raw 对象 |
| `out/*_merged_filtered_median5sd.rds` | 5 | 过滤后合并对象 |
| `out/clustering/nFmin200_nFmax_median5SD_mtmax3/` | 6 | 聚类主输出 |

Step 6 目录内：`clustree.png`、`*_harmony_prepped.rds`、`umap_*.png`、`*_harmony_res*.rds`、`clustering_summary.csv`。

样本 ID 与基因型见各 Rmd 中的 `samples` 表。
