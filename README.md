# LAWS-HiC

**L**ocally **A**daptive **W**eighting and **S**creening for **Hi-C** interaction calling.

LAWS-HiC is a computational tool designed to enhance the power of standard Hi-C peak callers by adjusting their initial results. Using a Locally Adaptive Weighting and Screening (LAWS) approach [1], LAWS-HiC adjusts p-values from standard Hi-C peak callers by accounting for the local spatial dependency structure of data, alleviating the common assumption of statistical independency among interactions made by most standard peak callers. This method enables more accurate identification of significant chromatin interactions.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("lb-zhou/LAWSHiC")
```

Dependencies: `data.table`, `foreach`, `doParallel` (installed automatically).

---

## Input files

### 1. Results from any standard Hi-C peak caller

A tab-separated file with a header row containing at least the following four columns (additional columns are preserved in the output):

| Column | Type | Description |
|--------|------|-------------|
| `chr` | integer | Chromosome number (e.g. `19` for chr19) |
| `fragmentMid1` | integer | Midpoint of the upstream fragment (bp) |
| `fragmentMid2` | integer | Midpoint of the downstream fragment (bp) |
| `p-value` | numeric | Interaction p-value |

**Example (mESC chr19, 10 kb resolution):**

```
chr  fragmentMid1  fragmentMid2  p-value
19   3285000       3265000       8.166582e-06
19   3345000       3325000       1.692066e-13
19   3355000       3325000       8.166582e-06
19   3355000       3335000       1.989453e-05
19   3365000       3345000       3.258155e-06
```

### 2. TAD file

A tab-separated or comma-separated file with a header row and three columns:

| Column | Type | Description |
|--------|------|-------------|
| `chr` | integer | Chromosome number |
| `start` | integer | TAD start coordinate (bp) |
| `end` | integer | TAD end coordinate (bp) |

**Example (mESC chr19):**

```
chr,start,end
19,0,4445000
19,4445000,5435000
19,5435000,6425000
19,6425000,7715000
19,7715000,8515000
19,8515000,9255000
```

---

## Usage

`laws_hic()` processes one chromosome at a time. Example data for mESC chr19 (4% downsampled, TADs 2–5 covering 4.4–8.5 Mb) are bundled with the package:

```r
library(LAWSHiC)

hic_file <- system.file("extdata", "example_hic_chr19.txt.gz", package = "LAWSHiC")
tad_file <- system.file("extdata", "example_tads_chr19.txt",   package = "LAWSHiC")

result <- laws_hic(
  chromosome_number = 19,
  input_file        = hic_file,
  tad_file          = tad_file,
  output_prefix     = "output/chr19",   # optional: write results to files
  ncores            = 4L
)

# Results are also returned as R objects
head(result$results)
head(result$diagnostics)
```

---

## Output

### Results table (`*_laws_results.txt`)

Contains all original input columns plus five new columns:

| Column | Description |
|--------|-------------|
| `i` | Bin index of `fragmentMid1` |
| `j` | Bin index of `fragmentMid2` |
| `LAWS_pvalue` | LAWS-adjusted p-value |
| `LAWS_pi` | Estimated local sparsity level. `NA` for bin pairs in skipped TADs. |
| `tad_skipped` | `TRUE` if the TAD was too short or had too few bin pairs; raw `p-value` is passed through unchanged. |

**Example output (GM12878 chr19, 1B reads, 4 rows):**

```
chr  fragmentMid1  fragmentMid2  p-value       i     j     LAWS_pvalue   LAWS_pi  tad_skipped
19   1585000       1435000       2.342e-01     159   144   1.000e+00     0.0500   FALSE
19   52485000      52465000      1.662e-04     5249  5247  9.409e-04     0.1501   FALSE
19   19075000      18915000      3.163e-01     1908  1892  1.000e+00     0.2200   FALSE
19   36595000      36445000      3.097e-01     3660  3645  5.052e-01     0.3800   FALSE
```

### Diagnostics table (`*_tad_diagnostics.txt`)

One row per (TAD × tau-grid point), recording the tau threshold, screening set size, and kernel weight for each of the 10 tau values.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `chromosome_number` | — | Chromosome number to process |
| `input_file` | — | Path to Hi-C p-value file |
| `tad_file` | — | Path to TAD file |
| `output_prefix` | `NULL` | If set, write results to `<prefix>_laws_results.txt` and `<prefix>_tad_diagnostics.txt` |
| `ncores` | `1` | Parallel cores for the bin-pair loop |
| `bin_size_bp` | `10000` | Hi-C resolution in bp |
| `min_tad_length_bp` | `200000` | Minimum TAD length to apply LAWS |
| `h_proportion` | `0.1` | Kernel bandwidth as a fraction of TAD span in bins |
| `tau_lower_limit` | `0.20` | Lower quantile of the p-value distribution used as the finest tau threshold |
| `tau_upper_limit` | `0.95` | Upper quantile of the p-value distribution used as the coarsest tau threshold |

---

## Down-sampling Hi-C data

The `scripts/` directory contains two R scripts for down-sampling raw Hi-C contact counts to simulate lower sequencing depths.

### Step 1 — Down-sample contact counts (`scripts/downsample.R`)

Takes long-format contact count files containing five columns (`chr1`, `mid1`, `chr2`, `mid2`, and `count`) and down-samples the total read count using multinomial sampling.

```bash
Rscript scripts/downsample.R \
  <INFDIR>            \  # directory containing per-chromosome input files
  <PREFIX>            \  # filename prefix for input (e.g. "" or "MAPQGE30_")
  <SUFFIX>            \  # filename suffix for input (e.g. "_10kb_RAWobserved.adj.fhic2.gz")
  <N_CHROMS>          \  # number of chromosomes (19 for mouse, 22 for human)
  <DRATIO>            \  # down-sampling fraction relative to original (e.g. 0.04)
  <DRATIO_LABEL>      \  # label used in output filename (e.g. "0.04")
  <OUTDIR>               # output directory
```

### Step 2 (optional) — Build FitHiC fragment files (`scripts/downsample_fragmentfile.R`)

Aggregates the down-sampled contact counts per fragment to create the marginal contact count files required by FitHiC2. This step is needed only when FitHiC2 is used for peak calling.

```bash
Rscript scripts/downsample_fragmentfile.R \
  <INFDIR>      \  # directory containing down-sampled files from Step 1
  <PREFIX>      \  # filename prefix (e.g. "chr")
  <SUFFIX>      \  # suffix of the down-sampled files (e.g. "_fhic2.downsampled_0.04.txt.gz")
  <N_CHROMS>    \  # number of chromosomes
  <OUTSUFFIX>      # suffix for output fragment files (e.g. "_fhic1.downsampled_0.04")
```

After these two steps, run FitHiC2 on the down-sampled inputs to obtain p-values, then supply the FitHiC2 output to `laws_hic()` as described above.

---

## References

[1] Cai, T. Tony, Wenguang Sun, and Yin Xia. “LAWS: A Locally Adaptive Weighting and Screening Approach to Spatial Multiple Testing.” Journal of the American Statistical Association 117, no. 539 (2021): 1370–83. doi:10.1080/01621459.2020.1859379.

## Contact

For errors, questions, or feedback, please contact Lingbo Zhou at  
[lzhou1@unc.edu](mailto:lzhou1@unc.edu).
