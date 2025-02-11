#' @title LAWS-Adjusted 3D Peak Calling P-Values
#' @description This function take peak caller results and topologically associating domains (TADs) list as input, adjusting the p-values by LAWS procedure then returns the LAWS-adjusted p-values and the estimated sparsity level.
#' @param input A data.table containing Hi-C interactions with columns: `chr`, `fragmentMid1`, `fragmentMid2`, and `p_value`.
#' @param domain_input A data.table containing TAD boundary information with columns: `chr`, `x1`, and `x2`.
#' @param chr A character string or interger indicating the chromosome to process, please make sure it is consistent with your input.
#' @param resolution A numeric value representing the resolution of Hi-C data in kilobases (kb).
#' @return A `data.table` containing the original peak calling results with two additional columns: `p_laws` (adjusted p-values) and `pi` (estimated sparsity level).
#' @import data.table
#' @import kedd
#' @export
laws_hic <- function(input, domain_input, chr, resolution) {

  library(kedd)
  library(data.table)

  kb <- as.numeric(resolution)

  # load the peak caller results:
  colnames(input) <- c("chr", "fragmentMid1", "fragmentMid2", "p_value")
  # assign bin pair identifiers:
  setDT(input)[, i:=(fragmentMid1+(kb*1e3/2))/(kb*1e3)][, j:=(fragmentMid2+(kb*1e3/2))/(kb*1e3)]

  # load TAD information:
  colnames(domain_input) <- c("chr", "x1", "x2")
  TAD_region_chr <- domain_input[domain_input$chr == chr, ]
  n_TAD <- nrow(TAD_region_chr)

  # define the data to store result:
  result_chr <- data.table(matrix(ncol = ncol(input) + 2, nrow = 0))
  setnames(result_chr, c(names(input), "p_laws", "pi"))

  # Iterate over each TAD region
  for(i in 1:n_TAD){

    # Identify Hi-C interactions within the current TAD region
    tad_ind <- which((input$fragmentMid2 >= TAD_region_chr$x1[i]) &
                       (input$fragmentMid1 <= TAD_region_chr$x2[i]))
    s_set_full <- input[tad_ind, ]
    s_set <- input[tad_ind, .(i, j, `p_value`)]  # Select minimal columns

    # Derive tau and bandwidth h for the current set
    q_tau <- 0.9
    tau <- bh.func(s_set$`p_value`, q_tau)$th

    m <- s_set$i
    n <- s_set$j

    h_est <- kedd::h.ccv(sqrt((m[1] - m)^2 + (n[1] - n)^2))$h
    h <- h_est

    cat("Processed TAD boundary ", i, "/", n_TAD,
        " with tau = ", tau, " and h = ", h, "\n")

    # Initialize a matrix to store results for each interaction
    temp <- matrix(ncol = 2, nrow = nrow(s_set))

    # Iterate over each interaction within the TAD
    for(k in 1:nrow(s_set)){

      bin_pair_i <- s_set$i[k]
      bin_pair_j <- s_set$j[k]

      # Calculate squared distances
      d2 <- (m - bin_pair_i)^2 + (n - bin_pair_j)^2

      # Compute the spatial weighting
      v_h <- exp(-(d2 / (h^2)) / 2)

      # Determine which p-values exceed the tau threshold
      tau_index <- s_set$`p_value` > tau

      # Calculate the estimated sparsity level (pi)
      estimated_pis_raw <- 1 - sum(v_h[tau_index]) / ((1 - tau) * sum(v_h))
      estimated_pis <- max(min(1 - 1e-5, estimated_pis_raw), 1e-5)

      # Compute the weight for the current interaction
      weight_s <- estimated_pis / (1 - estimated_pis)

      # Adjust the p-value based on the weight
      p_s <- (s_set$`p_value`[k]) / weight_s

      # Store the results, ensuring p_laws does not exceed 1
      temp[k, ] <- c(min(1, p_s), estimated_pis)
    }

    # Assign the computed p_laws and pi to the full set
    s_set_full[, `:=`(`p_laws` = temp[, 1], `pi` = temp[, 2])]

    # Combine the results
    result_chr <- rbindlist(list(result_chr, s_set_full))

  }

  return(result_chr)

}
