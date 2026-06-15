# threshold_diagnostics — note

**Location:** `<project>/out/threshold_diagnostics/`  
**When:** Step 4a–4b, before formal filtering (Step 5)  
**Input:** Unfiltered merged raw metadata (`nFeature_RNA`, `orig.ident`)

---

## Why these tables exist

Choosing `nFeature_min` and `nFeature_max` is a key QC decision. These CSVs summarize how nFeature is distributed and how many cells would be kept under different candidate cutoffs—so you can pick thresholds manually. They do **not** auto-recommend or rank options.

**Important:** The three main tables use **nFeature only**. Mitochondrial % and singlet/doublet filters are applied in Step 5. Retention % here will therefore be **higher** than after the full filter.

**Retention rule (same as Step 5 for nFeature):** keep cells with `nFeature > min` and `nFeature < max`.

Each table includes **overall** (all cells) and **per_sample** (`orig.ident`) rows where relevant.

---

## 1. `nfeature_distribution.csv`

**Purpose:** Describe the shape of the nFeature distribution before filtering.

| Column | Meaning |
|--------|---------|
| `scope` | `overall` or `per_sample` |
| `sample` | `all` or sample ID |
| `n_cells_total` | Cell count in that scope |
| `nFeature_median` | Median detected genes per cell |
| `nFeature_Q25`, `Q75`, `Q90` | 25th / 75th / 90th percentile |
| `nFeature_sd`, `nFeature_mad` | Standard deviation; median absolute deviation |
| `n_below_200`, `pct_below_200` | Cells with nFeature &lt; 200 |
| `n_at_or_above_Q75`, `pct_at_or_above_Q75` | Cells at or above the 75th percentile (~25% overall) |

**Use:** Check central tendency, spread, and how many low-nFeature cells exist. Compare `per_sample` rows—large differences mean one global cutoff may affect replicates unevenly.

---

## 2. `nfeature_sd_mad_bounds.csv`

**Purpose:** Data-driven candidate windows based on **median ± k × spread** (k = 2–5; spread = SD or MAD).

**Floor (200):** The lower bound is never below 200: `raw_lower = max(200, median − k × spread)`. Many rows show `raw_lower = 200` because the formula was truncated by this floor.

**Window:** The nFeature interval used to count retained cells—strict open bounds: `> lower` and `< upper`.

- **Raw window:** `raw_lower`, `raw_upper` (formula values; may be decimals)
- **Round window:** `round_min`, `round_max` (rounded to practical integer cutoffs, e.g. 200 and 2500)

| Column | Meaning |
|--------|---------|
| `method` | `median -/+ k * SD`, `median -/+ k * MAD`, or `count nFeature below floor` |
| `multiplier` | k (2–5); NA for the floor row |
| `spread` | SD or MAD used |
| `raw_lower`, `raw_upper` | Unrounded bounds |
| `round_min`, `round_max` | Rounded bounds |
| `n_cells_raw_window`, `pct_raw_window` | Cells / % in the raw window |
| `n_cells_round_window`, `pct_round_window` | Cells / % in the round window |

The **`count nFeature below floor`** row counts cells with nFeature &lt; 200 (not cells inside a retention window).

**Use:** See what spread-based intervals look like; compare raw vs rounded retention; inspect per-sample rows if replicates differ.

---

## 3. `nfeature_round_grid.csv`

**Purpose:** Lookup table for preset **integer** min/max pairs (e.g. 200/4000, 400/2500).

Default grid: min ∈ {200, 300, …, 1000}, max ∈ {2000, 2500, …, 6000}; only pairs with min &lt; max.

| Column | Meaning |
|--------|---------|
| `nFeature_min`, `nFeature_max` | Candidate lower / upper bounds |
| `n_cells_overall` | Cells in window (all samples) |
| `pct_retained_overall` | % retained overall |
| `min_pct_retained_any_sample` | Lowest % retained among samples |
| `max_pct_retained_any_sample` | Highest % retained among samples |

**Use:** Quickly compare round-number cutoffs. Pay attention to `min_pct_retained_any_sample` so no replicate is disproportionately lost.

---

## Optional: `user_filter_impact.csv`

Written only if `nFeature_min`, `nFeature_max`, and `mt_max` are set in the Rmd before Step 4a. Shows impact of **your chosen** thresholds, including mt and singlet where applicable.

---

## Quick workflow

1. **distribution** — How does nFeature look? Any bad replicates?
2. **sd_mad_bounds** — What do spread-based windows suggest?
3. **round_grid** — Pick integer min/max; check overall and worst-sample retention.
4. Set thresholds in the Rmd → **Step 5** applies nFeature + mt + singlet and writes `filter_summary/`.

**Regenerate:** Step 4a in the full pipeline, or `Rscript C:/Yale/scripts/run_threshold_diagnostics.R`
