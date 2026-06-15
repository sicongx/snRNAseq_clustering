# snRNA-seq median+5SD QC, SCTransform / LogNormalize, clustree

两个 snRNA-seq 项目的完整分析流程：**scn1lab** 与 **dyrk1a**。从 10x CellRanger 矩阵读入，经 QC、过滤、SCTransform（失败时自动 LogNormalize）、Harmony、clustree 与聚类，全部在 Rmd 中按步骤运行并重新生成结果。

## 目录结构

```
snRNAseq_median5sd_clustering/
├── scripts/                          # 共用 R 脚本
├── scn1lab_snRNAseq/
│   ├── scn1lab_full_pipeline_median5sd_clustree.Rmd
│   ├── data/                         # 10x 矩阵（每样本一文件夹）
│   └── out/                          # 运行后生成的输出（初始为空）
└── dyrk1a_snRNAseq/
    ├── dyrk1a_full_pipeline_median5sd_clustree.Rmd
    ├── data/
    └── out/
```

路径均相对仓库根目录 `snRNAseq_median5sd_clustering/`，无需修改本机绝对路径。

## 环境

- R ≥ 4.2，RStudio
- 建议内存 ≥ 64 GB（Step 6 合并 SCT 占用较大；内存不足时会自动改用 LogNormalize）

安装依赖（在 R 中运行一次）：

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

1. 在 RStudio 中将工作目录设为 `snRNAseq_median5sd_clustering/`（包含 `scripts/` 的文件夹）。
2. 打开要运行的 Rmd：
   - `scn1lab_snRNAseq/scn1lab_full_pipeline_median5sd_clustree.Rmd`
   - `dyrk1a_snRNAseq/dyrk1a_full_pipeline_median5sd_clustree.Rmd`
3. 默认从 Step 1 完整跑通（`run_step1_3_build_merge <- TRUE`）。可先 **Run All**，或按 chunk 逐步运行。
4. **Step 6 分三步**（Rmd 内 `pipeline-options` 可开关）：
   - **6a** `run_step6a_sctransform <- TRUE`：SCT + Harmony + UMAP
   - **6b** 6a 完成后设 `run_step6b_clustree <- TRUE`：生成 clustree 图
   - **6c** 根据 clustree 修改 `clustering_resolutions`，设 `run_step6c_clustering <- TRUE`：聚类与 UMAP 图

固定 QC 参数：`nFeature_min = 200`，`nFeature_max = 每样本 median + 5×SD`，`mt_max = 3`，仅保留 scDblFinder 单细胞。

## 产出结果

各项目输出在对应 `out/` 下：

| 路径 | 步骤 | 内容 |
|------|------|------|
| `out/qc_plots/` | 4 | 合并前 QC 图 |
| `out/threshold_diagnostics/` | 4a–4b | 阈值诊断 |
| `out/qc_filter_plots/` | 5 | 过滤前后对比 |
| `out/scn1_merged_raw.rds` 或 `out/dyrk1a_merged_raw.rds` | 3 | 合并 raw 对象 |
| `out/*_merged_filtered_median5sd.rds` | 5 | 过滤后合并对象 |
| `out/clustering/nFmin200_nFmax_median5SD_mtmax3/` | 6 | 聚类主输出 |

Step 6 聚类目录内主要包括：

- `clustree.png` — 分辨率 clustree
- `*_harmony_prepped.rds` — Harmony 后、聚类前对象
- `umap_*.png` — UMAP 图
- `*_harmony_res*.rds` — 各 resolution 的 Seurat 对象
- `clustering_summary.csv` — 聚类摘要

样本 ID 与基因型见各 Rmd 中的 `samples` 表。
