#' Run Subset-based Cross-tissue TWAS
#'
#' @param path A string of the direction for TWAS results.
#' @param pattern A string of the file pattern for TWAS results (default ".alldat").
#' @param cov_matrix A list of matrix of the gene expression covariance matrix across tissues from the reference panel (default using cov_matrix_GRCh37, can also change to cov_matrix_GRCh38 or make your own matrix list).
#' @param percent_act_tissue A decimal of the minimum percent of activated tissues for each gene regulated expression.
#' @param gene_list An array of the list of interested genes (default NULL; if NULL, it will go over all genes in the TWAS results; if not NULL, percent_act_tissue will be ignored).
#' @param n_more Simulation times for small p-values (default 1e+04; Caution: a very large number may lead to long calculation time; a very small number may lead to inaccurate p-value estimation).
#'
#' @return A dataframe for the Subset-based Cross-tissue TWAS results.
#' @export
#'
#' @examples
#' res_SCTWAS = run_SCTWAS("path_to_TWAS_resutls", cov_matrix)
run_SCTWAS = function(path,
                      cov_matrix = cov_matrix_GRCh37,
                      percent_act_tissue = 0.7,
                      n_more = 1e+04,
                      gene_list = NULL,
                      pattern = ".alldat") {
  # get a list of all TWAS results files
  files = list.files(path, pattern, full.names = TRUE)
  # select = dplyr::select
  # get a dataframe with all genes in the reference panel
  cat("Creating gene library......\n")
  all_gene = get_all_genes(files)
  cat("Done!\n")
  cat(nrow(all_gene), "genes identified. \n")
  # get TWAS z-values and p-values across all tissues
  cat("Loading TWAS results......\n")
  twas_z = all_gene
  twas_p = all_gene
  tissue_names = c("ID")
  # Add the prograss bar
  pb_tissue <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                total = length(files),
                                complete = "=",   # Completion bar character
                                incomplete = "-", # Incomplete bar character
                                current = ">",    # Current bar character
                                clear = FALSE,    # If TRUE, clears the bar when finish
                                width = 100)      # Width of the progress bar
  for (file in files) {
    pb_tissue$tick() # Update the prograss status
    name =  strsplit(basename(file), "[.]")[[1]][2] # extract the tissue name
    tissue_names = c(tissue_names, name)
    tem_dat = read.table(file, header = T) %>%
      filter(!is.na(TWAS.P) & TWAS.P != 0) %>%
      group_by(ID) %>%
      arrange(TWAS.P, by_group = T) %>%
      filter(n() == 1) %>%
      mutate(Z = qnorm(TWAS.P, lower.tail = F)) %>%
      dplyr::select(ID, Z, TWAS.P)
    twas_z = twas_z %>%
      left_join(tem_dat %>% dplyr::select(-TWAS.P), by = "ID")
    twas_z = twas_z %>% # select genes that as least expressed in certain percent of tissues
      filter(rowMeans(is.na(.)) < 1 - percent_act_tissue)
    twas_p = twas_p %>%
      left_join(tem_dat %>% dplyr::select(-Z), by = "ID")
    twas_p = twas_p %>% # select genes that as least expressed in certain percent of tissues
      filter(rowMeans(is.na(.)) < 1 - percent_act_tissue)
    # cat(which(file == files), "/", length(files), "\n")
  }
  names(twas_z) = tissue_names
  names(twas_p) = tissue_names
  cat("Done!\n")
  # performing subset-based test
  cat("Subset-based testing......\n")
  res = data.frame()
  if (!is.null(gene_list)) {
    genes = gene_list
  }
  else {
    genes = twas_z$ID
  }
  for (gene in genes) {
    cat("Testing gene:", gene, "......\n")
    if (! gene %in% twas_z$ID) {
      message("Gene ", gene, " is skipped due to too few activated tissues. Maybe select a lower cutoff of the minimum percent of activated tissues for gene regulated expression [percent_act_tissue]. \n")
      next
    }
    asset_z = data.frame(twas_z) %>%
      filter(ID == gene) %>%
      select_if(~ !is.na(.)) %>%
      dplyr::select(-"ID")
    asset_p = data.frame(twas_p) %>%
      filter(ID == gene) %>%
      select_if(~ !is.na(.)) %>%
      dplyr::select(-"ID")
    sub_out = get_sum_stats(asset_z)
    cat("Simulating p-value for gene:", gene, "......\n")
    n = 1e+03 # first doing 1000 simulations
    x = sim_null(gene, cov_matrix, n)
    list_stats = c()
    for (i in 1:n) {
      tem = get_sum_stats(x[i, ])
      list_stats = c(list_stats, tem$maxS)
    }
    # check whether p-value is small
    M = sum(list_stats >= sub_out$maxS)
    if (M < 10) {
      cat("Small p-value, simulating more replicates......\n")
      n = n_more # then simulating more times
      # Add the prograss bar
      pb_simu <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining (for current gene): :eta]",
                                  total = n,
                                  complete = "=",   # Completion bar character
                                  incomplete = "-", # Incomplete bar character
                                  current = ">",    # Current bar character
                                  clear = FALSE,    # If TRUE, clears the bar when finish
                                  width = 100)      # Width of the progress bar
      x = sim_null(gene, cov_matrix, n)
      list_stats = c()
      for (i in 1:n) {
        pb_simu$tick() # Update the prograss status
        tem = get_sum_stats(x[i, ])
        list_stats = c(list_stats, tem$maxS)
      }
      M = sum(list_stats >= sub_out$maxS)
    }
    if (M >= 10) {
      cat("Generating emperical p-value......\n")
      p_value = M/n
    }
    # tail approximation of p-value if p-value is small
    else {
      cat("Generating p-value with GPD estimation......\n")
      p_value = estp_gdp(250, list_stats, sub_out)
    }
    #sig_p = asset_p[, asset_p < 2.5e-6/length(asset_p)] # adjusted by tissue number
    tem = data.frame(Gene = gene,
                      Subset_Tissue = paste(sub_out$idx, collapse = ", "),
                      Number_of_Tissues = length(sub_out$idx),
                      P_value = p_value)
    res = rbind(res, tem)
    cat("Finished test for gene:", gene, "!\n")
    cat(which(gene == genes), "/", length(genes), "genes finished. \n")
  }
  cat("Job Done!\n")
  return(res)
}


