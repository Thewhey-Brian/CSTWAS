#' Run CSTWAS
#'
#' @param path A string of the direction for TWAS results.
#' @param pattern A string of the file pattern for TWAS results (default ".alldat").
#' @param cov_matrix A string indicating list of matrix of the gene expression covariance matrix across tissues from the reference panel (default using cov_matrix_GRCh37, can also change to cov_matrix_GRCh38 or use your own matrix list). This parameter is omitted if cov_matrix_path is specified.
#' @param percent_act_tissue A decimal of the minimum percent of activated tissues for each gene regulated expression.
#' @param gene_list An array of the list of interested genes (default NULL; if NULL, it will go over all genes in the TWAS results; if not NULL, percent_act_tissue will be ignored).
#' @param n_more Simulation times for small p-values (default 1e+04; Caution: a very large number may lead to long calculation time; a very small number may lead to inaccurate p-value estimation).
#' @param cov_matrix_path Path for downloaded reference gene expression covariance matrix across tissues (need to be named as "cov_matrix") (the reference matrix can be downloaded from: https://github.com/Thewhey-Brian/CSTWAS) If NULL, the function will automatically download the reference panel indicated by cov_matrix.
#'
#' @return cstwas_res: A dataframe for the CSTWAS results.
#'         meta_data: A dataframe for the tissue-specific TWAS results across multiple tissues.
#' @export
#'
#' @examples
#' res_CSTWAS = run_CSTWAS("path_to_TWAS_resutls", cov_matrix = "cov_matrix_GRCh37")
run_CSTWAS = function(path,
                      cov_matrix = "cov_matrix_GRCh37",
                      cov_matrix_path = NULL,
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
  cat(nrow(all_gene), "genes identified. These genes are expressed in one or multiple tissue(s) from TWAS resutls. \n")
  # get TWAS z-values and p-values across all tissues
  cat("Loading TWAS results......\n")
  twas_z = all_gene
  twas_p = all_gene
  twas_meta = data.frame(matrix(vector(), 0, 8))
  names(twas_meta) = c("ID", "Z", "TWAS.P", "CHR", "BP", "Tissue", "P0", "P1")
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
    pb_tissue$tick() # Update the progress status
    name =  strsplit(basename(file), "[.]")[[1]][2] # extract the tissue name
    tissue_names = c(tissue_names, name)
    tem_dat = read.table(file, header = T) %>%
      filter(!is.na(TWAS.P) & TWAS.P != 0) %>%
      group_by(ID) %>%
      arrange(TWAS.P, by_group = T) %>%
      filter(n() == 1) %>%
      mutate(Z = qnorm(TWAS.P, lower.tail = F),
             BP = (P0 + P1) / 2,
             Tissue = name) %>%
      dplyr::select(ID, Z, TWAS.P, CHR, BP, Tissue, P0, P1)
    twas_z = twas_z %>%
      left_join(tem_dat %>% dplyr::select(c(ID, Z)), by = "ID")
    twas_p = twas_p %>%
      left_join(tem_dat %>% dplyr::select(c(ID, TWAS.P)), by = "ID")
    twas_meta = rbind(twas_meta, tem_dat)
  }
  twas_z = twas_z %>% # select genes that as least expressed in certain percent of tissues
    filter(rowMeans(is.na(.)) < 1 - percent_act_tissue)
  twas_p = twas_p %>% # select genes that as least expressed in certain percent of tissues
    filter(rowMeans(is.na(.)) < 1 - percent_act_tissue)
  names(twas_z) = tissue_names
  names(twas_p) = tissue_names
  cat(nrow(twas_z), "genes identified to be activated in more than", percent_act_tissue * 100, "percent of tissues from the TWAS results.\n")
  cat("Done!\n")
  # loading reference gene expression covariance matrix across tissues
  if (is.null(cov_matrix_path)) {
    cat("Downloading reference gene expression covariance matrix", cov_matrix, "......\n")
    matrix_name = paste0(cov_matrix, ".rda")
    matrix_path = paste0("https://github.com/Thewhey-Brian/CSTWAS/blob/main/", cov_matrix, ".rda?raw=true")
    download.file(matrix_path, matrix_name)
    cat("Done!\n")
    cat("Loading reference gene expression covariance matrix", cov_matrix, "......\n")
    load(file.path(getwd(), matrix_name))
    cov_matrix = get(cov_matrix)
    cat("Done!\n")
  }
  else {
    cat("Loading reference gene expression covariance matrix ......\n")
    cov_matrix = get(load(cov_matrix_path))
    cat("Done!\n")
  }
  res = data.frame()
  if (!is.null(gene_list)) {
    genes = gene_list
  }
  else {
    genes = twas_z$ID
  }
  # find overlapped genes between TWAS results and the reference panel
  genes = genes[genes %in% names(cov_matrix)]
  cat("Among", nrow(twas_z), "activated genes", length(genes),"genes have reference information. \n")
  # performing subset-based test
  cat("Subset-based testing......\n")
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
  output = list(cstwas_res = res,
                meta_data = twas_meta)
  return(output)
}

#' Manhattan Plot For TWAS Results
#'
#' @param meta_data meta_data from run_CSTWAS results.
#' @param anot_index An integer indicating how significant results are to be annotated. (-log10(TWAS.P) > anot_index) This parameter will be ignored if anno_gene is not NULL.
#' @param ceiling_ctf An integer indicating how significant results are to be cut by the ceiling. (-log10(TWAS.P) > ceiling_ctf). If is NULL, it will automatically adjust based on the data.
#' @param floor_ctf An integer indicating how insignificant results are to be cut by the floor (-log10(TWAS.P) < floor_ctf). Default 0.
#' @param pts_size An integer indicating the point size.
#' @param anno_gene A list of genes that need to be annotated.
#' @param path Path for saving the plot.
#'
#' @return A Manhattan plot
#' @export
#'
#' @examples
#' mhp_twas(test$meta_data)
mhp_twas <- function(meta_data,
                     anot_index = 15,
                     ceiling_ctf = NULL,
                     floor_ctf = 0,
                     pts_size = 3.6,
                     anno_gene = NULL,
                     path = NULL) {
  dat_all = meta_data
  if(is.null(ceiling_ctf)) {
    ceiling_ctf = ceiling(-log10(min(dat_all$TWAS.P)))
  }
  if(is.null(anno_gene)) {
    dat_plot = dat_all %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%
      left_join(dat_all, ., by=c("CHR"="CHR")) %>%
      arrange(CHR, BP) %>%
      mutate(BPcum = BP + tot) %>%
      mutate(ceiling = ifelse(-log10(TWAS.P) > ceiling_ctf, "yes", "no")) %>%
      mutate(TWAS.P = ifelse(-log10(TWAS.P) > ceiling_ctf, 1e-1^ceiling_ctf, TWAS.P)) %>%
      mutate(is_annotate = ifelse(-log10(TWAS.P) > anot_index, "yes", "no"))
  }
  else{
    dat_plot <- dat_all %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%
      left_join(dat_all, ., by=c("CHR"="CHR")) %>%
      arrange(CHR, BP) %>%
      mutate(BPcum = BP + tot) %>%
      mutate(ceiling = ifelse(-log10(TWAS.P) > ceiling_ctf, "yes", "no")) %>%
      mutate(TWAS.P = ifelse(-log10(TWAS.P) > ceiling_ctf, 1e-1^ceiling_ctf, TWAS.P)) %>%
      mutate(is_annotate = ifelse(ID %in% anno_gene, "yes", "no"))
  }
  axisdf <- dat_plot %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  p = ggplot(dat_plot, aes(x=BPcum, y=-log10(TWAS.P))) +
    geom_point(aes(color=Tissue, shape = ceiling), alpha=0.8, size=pts_size) +
    scale_color_manual(values = sample(col_vector, length(unique(dat_all$Tissue)))) +
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(floor_ctf, ceiling_ctf)) +
    labs(x = "Chromosome", title = paste("Tissue-specific TWAS Manhattan Plot")) +
    geom_label_repel(data=subset(dat_plot, is_annotate=="yes"), aes(label=ID), size=5.2) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      legend.position="top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  if (!is.null(path)) {
    ggsave(path, width = 24, height = 12)
  }
  else {
    p
  }
}

#' Manhattan Plot For the CSTWAS Resutls
#'
#' @param meta_data meta_data from run_CSTWAS results.
#' @param cstwas_res cstwas_res from run_CSTWAS results.
#' @param anot_index An integer indicating how significant results are to be annotated. (-log10(TWAS.P) > anot_index) This parameter will be ignored if anno_gene is not NULL.
#' @param ceiling_ctf An integer indicating how significant results are to be cut by the ceiling. (-log10(TWAS.P) > ceiling_ctf). If is NULL, it will automatically adjust based on the data.
#' @param floor_ctf An integer indicating how insignificant results are to be cut by the floor (-log10(TWAS.P) < floor_ctf). Default 0.
#' @param anno_gene A list of genes that need to be annotated.
#' @param path Path for saving the plot.
#' @param pts_size An integer indicating the point size.
#'
#' @return A Manhattan plot
#' @export
#'
#' @examples
#' mhp_cstwas(test$meta_data, test$cstwas_res)
mhp_cstwas <- function(meta_data, cstwas_res,
                       anot_index = 15,
                       pts_size = 2,
                       ceiling_ctf = NULL,
                       floor_ctf = 0,
                       anno_gene = NULL,
                       path = NULL) {
  dat_all = meta_data
  dat_all = dat_all %>%
    ungroup() %>%
    mutate(Gene = ID) %>%
    select(-c(ID, Z, TWAS.P, Tissue)) %>%
    distinct() %>%
    left_join(cstwas_res %>% select(Gene, P_value), by = "Gene") %>%
    filter(!is.na(P_value)) %>%
    dplyr::rename(P = P_value)
  if(is.null(ceiling_ctf)) {
    ceiling_ctf = ceiling(-log10(min(dat_all$P)))
  }
  if(is.null(anno_gene)) {
    dat_plot <- dat_all %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%
      left_join(dat_all, ., by=c("CHR"="CHR")) %>%
      arrange(CHR, BP) %>%
      mutate(BPcum = BP + tot, CHR = as.factor(CHR)) %>%
      mutate(ceiling = ifelse(-log10(P) > ceiling_ctf, "yes", "no")) %>%
      mutate(P = ifelse(-log10(P) > ceiling_ctf, 1e-1^ceiling_ctf, P)) %>%
      mutate(is_annotate = ifelse(-log10(P) > anot_index, "yes", "no")) %>%
      distinct()
  }
  else{
    dat_plot <- dat_all %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%
      left_join(dat_all, ., by=c("CHR"="CHR")) %>%
      arrange(CHR, BP) %>%
      mutate(BPcum = BP + tot, CHR = as.factor(CHR)) %>%
      mutate(ceiling = ifelse(-log10(P) > ceiling_ctf, "yes", "no")) %>%
      mutate(P = ifelse(-log10(P) > ceiling_ctf, 1e-1^ceiling_ctf, P)) %>%
      mutate(is_annotate = ifelse(Gene %in% anno_gene, "yes", "no")) %>%
      distinct()
  }
  axisdf <- dat_plot %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  p = ggplot(dat_plot, aes(x=BPcum, y=-log10(P))) +
    geom_point(aes(color = CHR, shape = ceiling), alpha=0.8, size=pts_size) +
    scale_color_manual(values = sample(col_vector, 22)) +
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(floor_ctf, ceiling_ctf)) +
    labs(x = "Chromosome", title = paste("CSTWAS Manhattan Plot")) +
    geom_label_repel(data=subset(dat_plot, is_annotate=="yes"), aes(label=Gene), size=5.2) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
      legend.position="top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  if (!is.null(path)) {
    ggsave(path, width = 24, height = 12)
  }
  else {
    p
  }
}


#' Venn diagram for significant GReX associations
#'
#' @param meta_data meta_data from run_CSTWAS results.
#' @param cstwas_res cstwas_res from run_CSTWAS results.
#' @param path Path for saving the plot.
#' @param merge_range An integer indicating how wide (in base pairs) should be considered to merge nearby genes. Default +/- 1000bp.
#'
#' @return A Venn diagram
#' @export
#'
#' @examples
#' venn_diagram(test$meta_data, test$cstwas_res)
venn_diagram = function(meta_data,
                        cstwas_res,
                        merge_range = 1000,
                        path = NULL) {
  pos_data = meta_data %>%
    group_by(ID) %>%
    arrange(TWAS.P) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(Gene = ID,
           P0 = P0 - merge_range,
           P1 = P1 + merge_range,
           range = paste0(P0, "-", P1)) %>%
    select(Gene, P0, P1, TWAS.P, range)
  sc_list = cstwas_res %>%
    filter(P_value <= 2.5e-6) %>%
    left_join(pos_data, by = "Gene")
  ts_list = pos_data %>%
    filter(TWAS.P <= 2.5e-6 / 48) # adjusted by tissue number
  sc_list_ir = IRanges(start = sc_list$P0, end = sc_list$P1)
  ts_list_ir = IRanges(start = ts_list$P0, end = ts_list$P1)
  # reduce gene based on loci
  sc_list_gene = sc_list$Gene[
    sapply(mcols(reduce(sc_list_ir, with.revmap = T))$revmap,
                        "[[",
                        1)
    ]
  ts_list_gene = ts_list$Gene[
    sapply(mcols(reduce(ts_list_ir, with.revmap = T))$revmap,
           "[[",
           1)
  ]
  gene_vd_list = list("CSTWAS" = sc_list_gene,
                      "Tissue-specific TWAS" = ts_list_gene)
  p = ggVennDiagram(gene_vd_list, label = "count", set_size = 4) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    labs(title = "Venn Diagram of Significant GReX Associations") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  if (!is.null(path)) {
    ggsave(path, width = 24, height = 12)
  }
  else {
    p
  }
}

