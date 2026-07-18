#' Run LAWS-HiC on a single chromosome
#'
#' @description
#' Estimates the local null proportion \eqn{\hat{\pi}(i,j)} for each Hi-C bin
#' pair using a multi-tau Gaussian kernel smoother, then returns a
#' LAWS-adjusted p-value and local null proportion estimate for each pair.
#'
#' The tau grid is defined on the linear p-value scale as quantiles of the
#' empirical p-value distribution (excluding p = 1), which is mathematically
#' equivalent to the -log10 quantile approach used in the original implementation.
#'
#' @param chromosome_number Integer. Chromosome number to process
#'   (e.g. `19` for chr19).
#' @param input_file Character. Path to a tab-separated Hi-C file with a header
#'   row containing at least these four columns (in any order): `chr`,
#'   `fragmentMid1`, `fragmentMid2`, `p-value`. Additional columns are
#'   preserved in the output.
#' @param tad_file Character. Path to a tab-separated (or comma-separated) TAD
#'   file with a header row and three columns: `chr`, `start`, `end`.
#' @param output_prefix Character or `NULL` (default `NULL`). If provided, two
#'   tab-separated files are written:
#'   `<output_prefix>_laws_results.txt` and
#'   `<output_prefix>_tad_diagnostics.txt`.
#' @param ncores Integer. Number of parallel cores for the bin-pair loop
#'   (default `1`).
#' @param bin_size_bp Integer. Hi-C resolution in base pairs (default `10000`).
#'   Used to convert `fragmentMid` coordinates to bin indices.
#' @param min_tad_length_bp Integer. Minimum TAD length in bp required to apply
#'   LAWS. TADs shorter than this value, or with fewer than 20 bin pairs, are
#'   skipped; their raw p-values are passed through unchanged (default `200000`).
#' @param h_proportion Numeric. Gaussian kernel bandwidth as a fraction of the
#'   TAD span in bins. `h = max(2, tad_span_bins * h_proportion)` (default
#'   `0.1`).
#' @param tau_lower_limit Numeric. Lower quantile of the p-value distribution
#'   used as the finest (most stringent) tau threshold in the grid (default
#'   `0.20`). Only p-values below this quantile are considered potential
#'   signals when estimating pi at the finest scale.
#' @param tau_upper_limit Numeric. Upper quantile of the p-value distribution
#'   used as the coarsest tau threshold in the grid (default `0.95`). Only
#'   p-values below this quantile are considered potential signals when
#'   estimating pi at the coarsest scale.
#'
#' @return Invisibly returns a named list with two `data.table` objects:
#' \describe{
#'   \item{`results`}{Per-bin-pair table. Contains all original input columns
#'     plus five new columns: `i` (bin index of `fragmentMid1`), `j` (bin index
#'     of `fragmentMid2`), `LAWS_pvalue` (LAWS-adjusted p-value), `LAWS_pi`
#'     (estimated local null proportion), and `tad_skipped` (logical). For
#'     skipped TADs, `LAWS_pvalue` equals the raw `p-value`, `LAWS_pi` is `NA`,
#'     and `tad_skipped` is `TRUE`.}
#'   \item{`diagnostics`}{Per-TAD tau-grid diagnostics with columns
#'     `tad_idx`, `tad_start`, `tad_end`, `tad_length_bp`, `n_bin_pairs`,
#'     `tau_k_idx`, `tau_linear`, `T_k_size`, `omega`.}
#' }
#'
#' @export
#' @import data.table
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' result <- laws_hic(
#'   chromosome_number = 19,
#'   input_file        = "hic_chr19.txt",
#'   tad_file          = "tads.txt",
#'   output_prefix     = "output/chr19",
#'   ncores            = 4L
#' )
#' head(result$results)
#' head(result$diagnostics)
#' }
laws_hic <- function(
    chromosome_number,
    input_file,
    tad_file,
    output_prefix     = NULL,
    ncores            = 1L,
    bin_size_bp       = 10000L,
    min_tad_length_bp = 200000L,
    h_proportion      = 0.1,
    tau_lower_limit   = 0.20,
    tau_upper_limit   = 0.95
) {
  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "chromosome_number must be a single integer" =
      is.numeric(chromosome_number) && length(chromosome_number) == 1L,
    "input_file must be a character string" = is.character(input_file),
    "input_file does not exist"             = file.exists(input_file),
    "tad_file must be a character string"   = is.character(tad_file),
    "tad_file does not exist"               = file.exists(tad_file),
    "output_prefix must be NULL or a character" =
      is.null(output_prefix) || is.character(output_prefix)
  )

  # ── Constants ────────────────────────────────────────────────────────────────
  K_tau        <- 10L   # number of tau values in the grid
  gc_frequency <- 100L  # run gc() every N TADs

  chr         <- chromosome_number
  bin_size    <- bin_size_bp
  min_tad_len <- min_tad_length_bp

  registerDoParallel(cores = ncores)
  message("Cores detected: ", detectCores(), " | Using: ", ncores)

  # ── Load input ───────────────────────────────────────────────────────────────
  input <- fread(input_file, header = TRUE)
  stopifnot(
    "input_file must have columns: chr, fragmentMid1, fragmentMid2, p-value" =
      all(c("chr", "fragmentMid1", "fragmentMid2", "p-value") %in% names(input))
  )
  input <- input[chr == chromosome_number]
  message("Bin pairs on chr", chr, ": ", nrow(input))

  # Bin indices computed as vectors; not added to input to keep output clean
  half_bin  <- bin_size / 2L
  bin_i_all <- (input$fragmentMid1 + half_bin) / bin_size
  bin_j_all <- (input$fragmentMid2 + half_bin) / bin_size

  # ── Load TADs ────────────────────────────────────────────────────────────────
  TAD_all <- fread(tad_file, header = TRUE)
  stopifnot(
    "tad_file must have columns: chr, start, end" =
      all(c("chr", "start", "end") %in% names(TAD_all))
  )
  TAD_chr <- TAD_all[chr == chromosome_number]
  n_TAD   <- nrow(TAD_chr)
  message("TADs on chr", chr, ": ", n_TAD)

  # ── Pre-allocate result containers ───────────────────────────────────────────
  result_list <- vector("list", n_TAD)
  tad_tau_log <- vector("list", n_TAD)

  # ── Main TAD loop ────────────────────────────────────────────────────────────
  for (tad_i in seq_len(n_TAD)) {

    tad_start     <- TAD_chr$start[tad_i]
    tad_end       <- TAD_chr$end[tad_i]
    tad_length_bp <- tad_end - tad_start

    tad_ind    <- which(input$fragmentMid2 >= tad_start &
                          input$fragmentMid1 <= tad_end)
    s_set_full <- input[tad_ind]
    s_set_full[, i := bin_i_all[tad_ind]]
    s_set_full[, j := bin_j_all[tad_ind]]
    m_vec      <- s_set_full$i
    n_vec      <- s_set_full$j
    pvals      <- pmax(s_set_full$`p-value`, .Machine$double.xmin)
    n_pairs    <- nrow(s_set_full)

    # ── Small-TAD fallback: pass through raw p-values ─────────────────────────
    if (tad_length_bp < min_tad_len || n_pairs < 20L) {
      message(sprintf("Skip TAD %d/%d  [%d, %d]  length=%d bp  n_pairs=%d",
                      tad_i, n_TAD, tad_start, tad_end, tad_length_bp, n_pairs))
      s_set_full[, LAWS_pvalue := `p-value`]
      s_set_full[, LAWS_pi     := NA_real_]
      s_set_full[, tad_skipped := TRUE]
      result_list[[tad_i]] <- s_set_full
      next
    }

    # ── Step 1: tau grid on the p-value scale ────────────────────────────────
    # p = 0: replaced by smallest positive double to avoid division by zero.
    # p = 1: excluded from quantile computation (would inflate thresholds).
    pvals_valid    <- pvals[pvals < 1]
    tau_linear_vec <- quantile(
      pvals_valid,
      probs  = seq(tau_lower_limit, tau_upper_limit, length.out = K_tau),
      names  = FALSE
    )

    # ── Step 2: T_k sizes and omega weights ───────────────────────────────────
    T_k_sizes <- vapply(tau_linear_vec,
                        function(tau_k) sum(pvals > tau_k),
                        integer(1L))

    omega_raw              <- sqrt(T_k_sizes)
    omega_raw[T_k_sizes == 0L] <- 0
    omega_sum              <- sum(omega_raw)
    omega <- if (omega_sum == 0) rep(1 / K_tau, K_tau) else omega_raw / omega_sum

    # ── TAD-level diagnostics ─────────────────────────────────────────────────
    tad_tau_log[[tad_i]] <- data.frame(
      tad_idx       = tad_i,
      tad_start     = tad_start,
      tad_end       = tad_end,
      tad_length_bp = tad_length_bp,
      n_bin_pairs   = n_pairs,
      tau_k_idx     = seq_len(K_tau),
      tau_linear    = tau_linear_vec,
      T_k_size      = T_k_sizes,
      omega         = omega
    )

    # ── Step 3: kernel bandwidth ──────────────────────────────────────────────
    h <- max(2, (tad_end - tad_start) / bin_size * h_proportion)

    message(sprintf(
      "TAD %d/%d  [%d, %d]  length=%d bp  n_pairs=%d  h=%.3f  tau=[%.2e, %.2e]",
      tad_i, n_TAD, tad_start, tad_end, tad_length_bp, n_pairs,
      h, tau_linear_vec[1], tau_linear_vec[K_tau]
    ))

    tau_index_list <- lapply(tau_linear_vec, function(tau_k) pvals > tau_k)

    # ── Step 4: parallel pi estimation over bin pairs ─────────────────────────
    temp <- foreach(
      k        = seq_len(n_pairs),
      .combine = "rbind",
      .export  = c("m_vec", "n_vec", "pvals",
                   "tau_index_list", "tau_linear_vec", "omega", "h", "K_tau")
    ) %dopar% {

      d2   <- (m_vec - m_vec[k])^2 + (n_vec - n_vec[k])^2
      v_h  <- exp(-d2 / (2 * h^2))
      sv_h <- sum(v_h)

      pi_per_tau <- numeric(K_tau)
      for (kt in seq_len(K_tau)) {
        t_idx <- tau_index_list[[kt]]
        if (sum(t_idx) == 0L || sv_h == 0) {
          pi_per_tau[kt] <- 1e-5
        } else {
          pi_raw         <- 1 - sum(v_h[t_idx]) / ((1 - tau_linear_vec[kt]) * sv_h)
          pi_per_tau[kt] <- max(1e-5, min(1 - 1e-5, pi_raw))
        }
      }

      pi_hat <- max(1e-5, min(1 - 1e-5, sum(omega * pi_per_tau)))
      p_laws <- min(1, pvals[k] * (1 - pi_hat) / pi_hat)

      data.frame(LAWS_pvalue = p_laws, LAWS_pi = pi_hat, stringsAsFactors = FALSE)
    }

    s_set_full[, LAWS_pvalue := temp$LAWS_pvalue]
    s_set_full[, LAWS_pi     := temp$LAWS_pi]
    s_set_full[, tad_skipped := FALSE]
    result_list[[tad_i]] <- s_set_full

    if (tad_i %% gc_frequency == 0L) gc()

  }  # end TAD loop

  # ── Collect results ───────────────────────────────────────────────────────────
  result_chr      <- rbindlist(Filter(Negate(is.null), result_list))
  tad_diagnostics <- rbindlist(Filter(Negate(is.null), tad_tau_log))

  # ── Write output files (optional) ─────────────────────────────────────────────
  if (!is.null(output_prefix)) {
    f_results <- paste0(output_prefix, "_laws_results.txt")
    f_diag    <- paste0(output_prefix, "_tad_diagnostics.txt")
    fwrite(result_chr,      f_results, sep = "\t")
    fwrite(tad_diagnostics, f_diag,    sep = "\t")
    message("\nResults:     ", f_results)
    message("Diagnostics: ", f_diag)
  }

  message("Bin pairs processed: ", nrow(result_chr))

  invisible(list(results = result_chr, diagnostics = tad_diagnostics))
}
