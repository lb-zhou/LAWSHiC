library(testthat)
library(LAWSHiC)

test_that("laws_hic_paral works correctly", {
  # Prepare mock data
  data <- data.table::fread("/Users/waltchen/github-repo/LAWSHiC/example_data.tsv")
  domain <- data.table::fread("/Users/waltchen/github-repo/LAWSHiC/example_domain_list.tsv")
  resolution <- 10
  chr <- 21

  # Run the function
  result <- laws_hic_paral(data, domain, chr, resolution)
  
  # Check the result
  expected_col_names <- c("chr", "fragmentMid1", "fragmentMid2", "p_value", "i", "j", "p_laws", "pi")
  result_range <- range(result$fragmentMid1, result$fragmentMid2)
  domain_range <- range(domain$x1, domain$x2)

  expect_true(is.data.table(result))
  expect_true(all(expected_col_names %in% colnames(result)))
  
  # Check if results are within the TAD region using range
  expect_true(result_range[2] <= domain_range[2] & result_range[1] >= domain_range[1], 
              info = "Result bins should be within the TAD region")
              
})
