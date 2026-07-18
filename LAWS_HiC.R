#!/usr/bin/env Rscript
# =============================================================================
# LAWS-HiC: Local Adaptive Weighted Smoothing for Hi-C
#
# Estimates the local null proportion pi(i,j) for each bin pair using a
# multi-tau Gaussian kernel smoother, then applies weighted BH FDR control.
#
# Input:
#   --input_file   TSV with header, columns: chr  fragmentMid1  fragmentMid2  p-value
#   --tad_file     TSV with header, columns: chr  start  end
#
# Output (two files written to --output_prefix):
#   <prefix>_laws_results.txt         per-bin-pair results
#   <prefix>_tad_diagnostics.txt      per-TAD tau-grid diagnostics
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(foreach)
  library(doParallel)
})

# ── Helper functions ──────────────────────────────────────────────────────────

fdp_with_qvalue <- function(pv, pi, q) {
  # Weighted BH FDR control with q-value computation.
  # pv : LAWS-adjusted p-values
  # pi : local null proportion estimates (same length as pv)
  # q  : nominal FDR level
  m      <- length(pv)
  st.pv  <- sort(pv)
  ord    <- order(pv)

  # Adjusted statistic: (sum of (1-pi_k)) * p_(i) / i
  pvi    <- sum(1 - pi) * st.pv / seq_len(m)

  # Q-values: expected FDR at each sorted rank
  qvals  <- numeric(m)
  for (i in seq_len(m)) {
    qvals[i] <- min(1, sum(1 - pi) * st.pv[i] / i)
  }
  # Enforce monotonicity (non-decreasing from most to least significant)
  for (i in (m - 1):1) qvals[i] <- min(qvals[i], qvals[i + 1])

  # Map back to original order
  qvals_ordered        <- numeric(m)
  qvals_ordered[ord]   <- qvals

  # FDR decision
  de <- integer(m)
  if (sum(pvi <= q) == 0) {
    k  <- 0L
    pk <- 1
  } else {
    k  <- max(which(pvi <= q))
    pk <- st.pv[k]
    de[pv <= pk] <- 1L
  }

  list(nr = k, th = pk, de = de, qvalue = qvals_ordered)
}

# ── Command-line options ──────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-c", "--chromosome_number"), type = "integer", default = NA,
              help = "Chromosome number to process (integer, e.g. 1 for chr1)."),

  make_option("--input_file", type = "character", default = NA,
              help = "Path to Hi-C input file (TSV, header required).
                      Expected columns: chr  fragmentMid1  fragmentMid2  p-value"),

  make_option("--tad_file", type = "character", default = NA,
              help = "Path to TAD list file (TSV, header required).
                      Expected columns: chr  start  end"),

  make_option("--output_prefix", type = "character", default = NA,
              help = "Output file prefix (directory + base name, no extension).
                      Two files are written: <prefix>_laws_results.txt and
                      <prefix>_tad_diagnostics.txt"),

  make_option("--ncores", type = "integer", default = 1L,
              help = "Number of parallel cores for the bin-pair loop [default: 1]."),

  make_option("--bin_size_bp", type = "integer", default = 10000L,
              help = "Resolution of the Hi-C data in base pairs [default: 10000].
                      Used to convert fragmentMid coordinates to bin indices."),

  make_option("--min_tad_length_bp", type = "integer", default = 200000L,
              help = "Minimum TAD length (bp) to apply LAWS. Smaller TADs and
                      TADs with fewer than 20 bin pairs are skipped; raw p-values
                      are passed through unchanged [default: 200000]."),

  make_option("--h_proportion", type = "numeric", default = 0.1,
              help = "Kernel bandwidth as a fraction of the TAD span in bins.
                      h = max(2, tad_span_bins * h_proportion) [default: 0.1]."),

  make_option("--tau_lower_limit", type = "numeric", default = 0.05,
              help = "Lower quantile of the p-value distribution used to define
                      the coarsest tau threshold [default: 0.05]."),

  make_option("--tau_upper_limit", type = "numeric", default = 0.80,
              help = "Upper quantile of the p-value distribution used to define
                      the finest tau threshold [default: 0.80].")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Validate required arguments ───────────────────────────────────────────────

stopifnot(
  "chromosome_number is required" = !is.na(opt$chromosome_number),
  "input_file is required"        = !is.na(opt$input_file),
  "tad_file is required"          = !is.na(opt$tad_file),
  "output_prefix is required"     = !is.na(opt$output_prefix)
)

cat(
  "=== LAWS-HiC ===\n",
  "chromosome:      ", opt$chromosome_number, "\n",
  "input_file:      ", opt$input_file,        "\n",
  "tad_file:        ", opt$tad_file,          "\n",
  "output_prefix:   ", opt$output_prefix,     "\n",
  "ncores:          ", opt$ncores,            "\n",
  "bin_size_bp:     ", opt$bin_size_bp,       "\n",
  "min_tad_length:  ", opt$min_tad_length_bp, "\n",
  "h_proportion:    ", opt$h_proportion,      "\n",
  "tau_lower_limit: ", opt$tau_lower_limit,   "\n",
  "tau_upper_limit: ", opt$tau_upper_limit,   "\n",
  sep = ""
)

chr          <- opt$chromosome_number
bin_size     <- opt$bin_size_bp
min_tad_len  <- opt$min_tad_length_bp
h_proportion <- opt$h_proportion
tau_q_low    <- opt$tau_lower_limit
tau_q_high   <- opt$tau_upper_limit
K_tau        <- 10L           # number of tau values in the grid
q_fdr        <- 0.05          # nominal FDR level
gc_frequency <- 100L          # GC every N TADs

registerDoParallel(cores = opt$ncores)
cat("Cores detected: ", detectCores(), " | Using: ", opt$ncores, "\n", sep = "")

# ── Load input ────────────────────────────────────────────────────────────────

input <- fread(opt$input_file, header = TRUE)
stopifnot(
  "input must have columns: chr, fragmentMid1, fragmentMid2, p-value" =
    all(c("chr", "fragmentMid1", "fragmentMid2", "p-value") %in% names(input))
)

# Filter to requested chromosome
input <- input[chr == opt$chromosome_number]
cat("Bin pairs on chr", chr, ":", nrow(input), "\n")

# Compute bin indices (bin centre = fragmentMid)
half_bin <- bin_size / 2L
input[, i   := (fragmentMid1 + half_bin) / bin_size]
input[, j   := (fragmentMid2 + half_bin) / bin_size]
input[, bin := paste(chr, fragmentMid1, fragmentMid2, sep = ":")]

# ── Load TADs ─────────────────────────────────────────────────────────────────

TAD_all <- fread(opt$tad_file, header = TRUE)
stopifnot(
  "tad_file must have columns: chr, start, end" =
    all(c("chr", "start", "end") %in% names(TAD_all))
)
TAD_chr <- TAD_all[chr == opt$chromosome_number]
n_TAD   <- nrow(TAD_chr)
cat("TADs on chr", chr, ":", n_TAD, "\n")

# ── Pre-allocate results ──────────────────────────────────────────────────────

result_chr <- data.table(matrix(ncol = ncol(input) + 3L, nrow = 0L))
setnames(result_chr, c(names(input), "LAWS_qvalue", "LAWS_pvalue", "LAWS_pi"))

tad_tau_log <- vector("list", n_TAD)

# =============================================================================
# Main TAD loop
# =============================================================================

for (tad_i in seq_len(n_TAD)) {

  tad_start     <- TAD_chr$start[tad_i]
  tad_end       <- TAD_chr$end[tad_i]
  tad_length_bp <- tad_end - tad_start

  # Bin pairs whose span overlaps this TAD
  tad_ind    <- which(input$fragmentMid2 >= tad_start & input$fragmentMid1 <= tad_end)
  s_set_full <- input[tad_ind]
  s_set      <- s_set_full[, .(i, j, `p-value`)]

  # ── Small-TAD fallback ──────────────────────────────────────────────────────
  if (tad_length_bp < min_tad_len || nrow(s_set) < 20L) {
    cat(sprintf("Skip TAD %d/%d  [%d, %d]  length=%d bp  n_pairs=%d\n",
                tad_i, n_TAD, tad_start, tad_end, tad_length_bp, nrow(s_set)))

    s_set_full[, `:=`(
      LAWS_pvalue = `p-value`,
      LAWS_pi     = NA_real_,
      LAWS_qvalue = NA_real_,
      tad_skipped = TRUE
    )]
    for (kt in seq_len(K_tau))
      s_set_full[, paste0("LAWS_pi_tau_", kt) := NA_real_]

    result_chr <- rbindlist(list(result_chr, s_set_full), fill = TRUE)
    next
  }

  s_set_full[, tad_skipped := FALSE]

  # ── Step 1: guard p=0, compute tau grid on linear p-value scale ─────────────
  # p=0 is replaced by the smallest positive double to avoid log(0) or /0.
  # p=1 pairs are excluded from the quantile computation because they carry no
  # signal about the bulk null and would inflate the tau thresholds.
  # The tau grid is derived directly from quantiles of p (no -log10 round-trip):
  #   original logic:  tau_star  = quantile(-log10(p_valid), seq(low, high))
  #                    tau_linear = 10^(-tau_star)
  #   equivalent:      tau_linear = quantile(p_valid, seq(1-high, 1-low))
  # (-log10 is strictly decreasing, so probability order inverts.)
  pvals       <- pmax(s_set$`p-value`, .Machine$double.xmin)
  pvals_valid <- pvals[pvals < 1]

  tau_linear_vec <- quantile(
    pvals_valid,
    probs  = seq(1 - tau_q_high, 1 - tau_q_low, length.out = K_tau),
    names  = FALSE
  )

  # ── Step 2: compute T_k sizes and omega weights ─────────────────────────────
  # T(tau_k) = { (i,j) : p_ij > tau_k }  (the "null-looking" screening set)
  T_k_sizes <- vapply(tau_linear_vec,
                      function(tau_k) sum(pvals > tau_k),
                      integer(1L))

  # omega_k ∝ sqrt(|T_k|)  (inverse-variance weighting approximation)
  omega_raw             <- sqrt(T_k_sizes)
  omega_raw[T_k_sizes == 0L] <- 0
  omega_sum             <- sum(omega_raw)
  omega <- if (omega_sum == 0) rep(1 / K_tau, K_tau) else omega_raw / omega_sum

  # ── TAD-level diagnostics ───────────────────────────────────────────────────
  tad_tau_log[[tad_i]] <- data.frame(
    tad_idx       = tad_i,
    tad_start     = tad_start,
    tad_end       = tad_end,
    tad_length_bp = tad_length_bp,
    n_bin_pairs   = nrow(s_set),
    tau_k_idx     = seq_len(K_tau),
    tau_linear    = tau_linear_vec,
    T_k_size      = T_k_sizes,
    omega         = omega
  )

  # ── Step 3: bandwidth h ─────────────────────────────────────────────────────
  tad_span_bins <- (tad_end - tad_start) / bin_size
  h             <- max(2, tad_span_bins * h_proportion)

  m_vec <- s_set$i
  n_vec <- s_set$j

  cat(sprintf(
    "TAD %d/%d  [%d, %d]  length=%d bp  n_pairs=%d  h=%.3f  tau=[%.2e, %.2e]\n",
    tad_i, n_TAD, tad_start, tad_end, tad_length_bp, nrow(s_set),
    h, tau_linear_vec[1], tau_linear_vec[K_tau]
  ))

  # Pre-compute per-tau screening index vectors (shared across workers)
  tau_index_list <- lapply(tau_linear_vec, function(tau_k) pvals > tau_k)

  # ── Step 4: parallel pi estimation over bin pairs ────────────────────────────
  temp <- foreach(
    k        = seq_len(nrow(s_set)),
    .combine = "rbind",
    .export  = c("m_vec", "n_vec", "pvals",
                 "tau_index_list", "tau_linear_vec", "omega", "h", "K_tau")
  ) %dopar% {

    # Gaussian kernel weights centred on bin pair k
    d2   <- (m_vec - m_vec[k])^2 + (n_vec - n_vec[k])^2
    v_h  <- exp(-d2 / (2 * h^2))
    sv_h <- sum(v_h)

    # Per-tau local null proportion estimate
    pi_per_tau <- numeric(K_tau)
    for (kt in seq_len(K_tau)) {
      t_idx <- tau_index_list[[kt]]
      if (sum(t_idx) == 0L || sv_h == 0) {
        pi_per_tau[kt] <- 1e-5          # conservative fallback
      } else {
        pi_raw         <- 1 - sum(v_h[t_idx]) / ((1 - tau_linear_vec[kt]) * sv_h)
        pi_per_tau[kt] <- max(1e-5, min(1 - 1e-5, pi_raw))
      }
    }

    # Weighted average across tau grid
    pi_hat <- max(1e-5, min(1 - 1e-5, sum(omega * pi_per_tau)))

    # LAWS-adjusted p-value
    p_laws <- min(1, pvals[k] * (1 - pi_hat) / pi_hat)

    out <- data.frame(LAWS_pvalue = p_laws, LAWS_pi = pi_hat,
                      stringsAsFactors = FALSE)
    for (kt in seq_len(K_tau))
      out[[paste0("pi_tau_", kt)]] <- pi_per_tau[kt]
    out
  }

  # ── Step 5: TAD-level FDR control ───────────────────────────────────────────
  fdp_res <- fdp_with_qvalue(temp$LAWS_pvalue, temp$LAWS_pi, q_fdr)

  s_set_full[, `:=`(
    LAWS_qvalue = fdp_res$qvalue,
    LAWS_pvalue = temp$LAWS_pvalue,
    LAWS_pi     = temp$LAWS_pi
  )]
  for (kt in seq_len(K_tau))
    s_set_full[, paste0("LAWS_pi_tau_", kt) := temp[[paste0("pi_tau_", kt)]]]

  result_chr <- rbindlist(list(result_chr, s_set_full), fill = TRUE)

  if (tad_i %% gc_frequency == 0L) gc()

}  # end TAD loop

# ── Write outputs ─────────────────────────────────────────────────────────────

tad_diagnostics <- rbindlist(Filter(Negate(is.null), tad_tau_log))

fwrite(result_chr,      paste0(opt$output_prefix, "_laws_results.txt"),
       sep = "\t", row.names = FALSE)
fwrite(tad_diagnostics, paste0(opt$output_prefix, "_tad_diagnostics.txt"),
       sep = "\t", row.names = FALSE)

cat("\nDone.\n")
cat("Results:     ", paste0(opt$output_prefix, "_laws_results.txt"),     "\n")
cat("Diagnostics: ", paste0(opt$output_prefix, "_tad_diagnostics.txt"),  "\n")
cat("Bin pairs processed:", nrow(result_chr), "\n")
