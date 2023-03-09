#' Get all genes
#'
#' @param files a list of all TWAS results files
#'
#' @return a dataframe with all gene symbols (IDs)
#' @export
#'
#' @examples
#' gene_list <= get_all_genes(files)
get_all_genes = function(files) {
  all_gene = data.frame()
  for (file in files) {
    all_gene = rbind(all_gene,
                     read.table(file, header = T) %>% dplyr::select(ID)) %>%
      distinct(ID)
  }
  return(all_gene)
}

#' Get the test statistics
#'
#' @param Z z-values transformed from TWAS p-values
#'
#' @return a list with test statistics (MaxS), corresponding selected subset (idx), and all possible statistics
#' @export
#'
#' @examples
#' z <- 1:3
#' names(z) <- c("a", "b", "c")
#' test_stats <- get_sum_stats(z)
get_sum_stats <- function(Z){
  names_list <- names(Z)
  Z <- as.numeric(Z)
  names(Z) <- names_list
  K <- length(Z)
  Z1 <- sort(Z, decreasing = TRUE)
  S <- rep(0,K)
  for(i in 1:K){
    S[i] <- sum(Z1[1:i]) / sqrt(i)
  }
  id_max <- which(S == max(S, na.rm = T))
  res <- list(maxS = S[id_max], idx = names(Z1[1:id_max]), S = S)
  return(res)
}

#' Simulate expression z-values from null gene regulated expression distribution across tissues
#'
#' @param gene The gene name
#' @param cov_matrix The null gene regulated expression covariance matrix across tissues for all genes in the reference panel
#' @param n Number of simulations
#'
#' @return A nxm matrix, where n is the number of samples (simulated replicated z-values) for all expression-activated tissues, and m is the tissue number
#' @export
#'
#' @examples
#' x <- sim_null("APOE", cov_matrix, 10000)
sim_null <- function(gene, cov_matrix, n) {
  cov_gene <- cov_matrix[[gene]]
  return(mvrnorm(n, rep(0, nrow(cov_gene)), cov_gene))
}

#' Small p-values estimation with GDP mle estimators
#'
#' @param Nexc Exceedances threshold (default 250)
#' @param list_stats List of statistics from simulations
#' @param sub_out Result from get_sum_stats function
#'
#' @return An estimated p-value
#' @export
#'
#' @examples
#' list_stats <- c(0.25, 1.69, 1.86, 2.3, 0.01)
#' p_value <- estp_gdp(250, list_stats, sub_out)
estp_gdp <- function(Nexc = 250, list_stats, sub_out) {
  n <- length(list_stats)
  y_star <- sort(list_stats, decreasing = T)
  t <- (y_star[Nexc] + y_star[Nexc + 1]) / 2
  z <- y_star[y_star - t > 0] - t
  x_bar <- mean(z)
  x_var <- var(z)
  # mle estimators of GPD distribution
  alpha <- x_bar * (x_bar ^ 2 / x_var + 1) / 2
  k <- (x_bar ^ 2 / x_var - 1) / 2
  # adjusting for unsolvable situation
  while (sub_out$maxS > alpha / k & k > 0) {
    z <- c(z, sub_out$maxS)
    x_bar <- mean(z)
    x_var <- var(z)
    alpha <- x_bar * (x_bar ^ 2 / x_var + 1) / 2
    k <- (x_bar ^ 2 / x_var - 1) / 2
    cat(sub_out$maxS, ">= a/k")
  }
  if (k == 0) {
    omFz <- exp(-(sub_out$maxS - t) / alpha)
  }
  else {
    omFz <- (1 - k * (sub_out$maxS - t) / alpha) ^ (1 / k)
  }
  p_value <- Nexc / n * omFz
  return(p_value)
}

#' Convert gene IDs
#'
#' @param twas a dataframe of TWAS results with Ensemble ID
#'
#' @return a dataframe of TWAS results with gene symble
#' @export
#'
#' @examples
#' twas_z = data.frame(ID = "ENSG00000261456.5", Z = "0.15066")
#' twas_z <= convert_genes_ids(twas_z)
convert_genes_ids = function(twas) {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_list <- gsub("\\..", "", twas$ID)
  dat_out <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   values = gene_list, mart = mart) %>%
    rename(ID = ensembl_gene_id,
           Gene = hgnc_symbol) %>%
    left_join(twas %>%
                mutate(ID = gsub("\\..", "", ID)),
              by = "ID") %>%
    mutate(ID = Gene) %>%
    dplyr::select(-Gene) %>%
    distinct(ID, .keep_all = T)
  return(dat_out)
}
