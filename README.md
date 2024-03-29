# CSTWAS
Transcriptome-wide association study (TWAS) is introduced to identify 
    significant expression-trait associations through imputations. It has been widely 
    used to analyze tissue-specific associations with the reference expression 
    quantitative trait loci (eQTL) panel. To increase the statistical power of TWAS results, 
    meta-analysis methods aggregating TWAS results across multiple tissues are developed. 
    However, most existing meta-analysis methods lose interpretation of disease etiology and 
    have limited power to identify weaker associations when only a few tissues are weakly 
    activated. Therefore, we developed the cross-tissue subset-based meta-analysis method, 
    also called cross-tissue subset-based transcriptome-wide association study (CSTWAS). 
    In this package, we aggregate the TWAS results across tissues and perform meta-analysis 
    through the subset-based test. R functions are provided for researchers to integrate 
    TWAS results across multiple tissues and visualize the result.

# Installation

Use the following codes to install the CSTWAS package
```R
library(devtools)
install_github("Thewhey-Brian/CSTWAS")
```
For details about how to install a R package directly from GitHub: https://rdrr.io/cran/remotes/man/install_github.html.

# Usage

## Prerequisites
In order to integrate **TWAS results** across multiple tissues, the ideal format of the TWAS results shown as following:

- TWAS results for each tissue should contain the following variables: `ID`, `CHR`, `P0`, `P1`, `TWAS.Z`, `TWAS.P`. The additional variables will be ignored.
  
  |   ID   |  CHR  |     P0    |    P1     |  TWAS.Z   |  TWAS.P  | 
  | ------ | :---: | :-------- | :-------- | --------- | -------- | 
  | ALOX5  |  10   | 45869661  | 45941561  | 0.573537  | 0.566281 |
  | COX15  |  10   | 101471601 | 101491857 | 1.539744  | 0.123623 |
  | ZDHHC6 |  10   | 114190058 | 114206672 | 0.834107  | 0.404221 |
  | ABCC2  |  10   | 101542489 | 101611949 | 1.540492  | 0.123440 |
  | VIM.   |  10   | 17270258  | 17279592  | 0.484861  | 0.627775 |
  | FAS    |  10   | 90750414  | 90775542  | -0.925617 | 0.354645 |
  
  -  `ID`: Feature/gene identifier
  -  `CHR`: Chromosome
  -  `P0`: Gene start
  -  `P1`: Gene end
  -  `TWAS.Z`: TWAS Z-score 
  -  `TWAS.P`: TWAS P-value
  
- TWAS results for each tissue are stored in a single file, named as "all.tissue_name.alldat". And all tissue-specific TWAS results files are saved in the same folder/path. 

  For example, TWAS results for Adipose Coronary should be saved as "all.Adipose_Coronary.alldat". 
  
  The resutls files structure looks like:
  
  ![Example of results files structure](/Plots/Files_structure_example.png)

### TWAS results generated from the [FUSION](http://gusevlab.org/projects/fusion/) software. 
Here is an example of formating TWAS results from the FUSION software. Other TWAS software such as [PrediXcan](https://github.com/hakyimlab/PrediXcan) can also be used. 

**Inputs:**
1. GWAS summary statistics (file: ***dise.sumstats***)
   - `SNP`: SNP identifier (rsID)
   - `A1`: first allele (effect allele)
   - `A2`: second allele (other allele)
   - `Z`: Z-scores, sign with respect to `A1`
2. Expression weights (stored in: ***./WEIGHTS/***) ([Download](http://gusevlab.org/projects/fusion/#download-pre-computed-predictive-models))
3. LD reference data--1000 Genomes (stored in: ***./LDREF/1000G.EUR.***) ([Download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2))
4. Tissue list (file: ***tissue_list.txt***)
   - A list of tissue names from expression weights (one tissue name per line)
   > Adipose_Subcutaneous
   > 
   > Brain_Amygdala
   > 
   > Esophagus_Mucosa
   > 
   > Heart_Left_Ventricle
   > 
   > Lung
   > 
   > Nerve_Tibial

**Performing TWAS with FUSION:**
```bash
dise=$1 # Disease name
tiss=$2 # File stores a list of tissues

module load R

values=$(cat $tiss)
mkdir /outcomes/$dise # create outcomes path
for tissue in $values # loop through all tissues
do
	for i in $(seq 1 22) # loop through all 22 chromosomes
	do # create .sh files for all combinations
		echo "module load R;
		      Rscript /Tools/fusion_twas-master/FUSION.assoc_test.R \ # the location of installed FUSION software
		      --sumstats /data/$dise.sumstats \ # the GWAS summary statistics
		      --weights /reference/WEIGHTS/$tissue.pos \ # tissue-specific reference expression weight
		      --weights_dir /reference/WEIGHTS/ \ # the location of reference expression weight
		      --ref_ld_chr /reference/LDREF/1000G.EUR. \ # LD reference data
		      --chr $i \ # specify chromosome
		      --out /outcomes/$dise/$dise_${tissue}_${i}.dat" > $dise_${tissue}_${i}.sh
		qsub -cwd -l mem_free=20G,h_vmem=30G $dise_${tissue}_${i}.sh # submit jobs
		sleep 1
	done
done
```
**Formatting Outpus:**

Since FUSION can only perform chromosome-specific analysis, we need to combine the results from all 22 chromosomes for each tissue.

```bash
out_dir=$1 # Output direction
file=$2 # File stores a list of tissues

tissues=$(cat $file)
for tissue in $tissues
do
	head -n 1 ${out_dir}/${tissue}_1.dat | awk '{print $1, $3, $4, $5, $6, $19, $20}' > ${out_dir}/all.${tissue}.alldat
	tail -q -n +2 ${out_dir}/${tissue}*.dat | awk '{print $1, $3, $4, $5, $6, $19, $20}' >> ${out_dir}/all.${tissue}.alldat
	echo "Finish ${tissue}, length: $(wc -l ${out_dir}/all.${tissue}.alldat)"
done
```

## Run CSTWAS

With the formatted TWAS results, we are ready to perform the CSTWAS anslysis.

!Note: as weights labeled in [FUSION TWAS website](http://gusevlab.org/projects/fusion/#gtex-v8-multi-tissue-expression), cov_matrix_GTEx_v7 uses gene symble while cov_matrix_GTEx_v8 uses Ensembl ID. Please be careful about the ID consistency for TWAS results and gene_list (if needed.) 

If using cov_matrix_GTEx_v8 as reference panel and want gene ID in the CSTWAS results as gene symble, simply set `run_CSTWAS(..., gene_symble = TRUE)`.

### run_CSTWAS: Run CSTWAS

***Inputs:***
- `path`: A string of the direction for TWAS results.
- `cov_matrix`: A string indicating list of matrix of the gene expression covariance matrix across tissues from the reference panel (default using cov_matrix_GTEx_v7, can also change to cov_matrix_GTEx_v8 or use your own matrix list.) This parameter is omitted if cov_matrix_path is specified. 
- `cov_matrix_path`: Path for downloaded reference gene expression covariance matrix across tissues (the reference matrix can be downloaded from: https://github.com/Thewhey-Brian/CSTWAS) If NULL, the function will automatically download the reference panel indicated by cov_matrix from the GitHub repository.
- `percent_act_tissue`: A decimal of the minimum percent of activated tissues for each gene regulated expression.
- `n_more`: Simulation times for small p-values (default 1e+04; Caution: a very large number may lead to long calculation time; a very small number may lead to inaccurate p-value estimation).
- `gene_list`: An array of the list of interested genes (default NULL; if NULL, it will go over all genes in the TWAS results; if not NULL, percent_act_tissue will be ignored).
- `pattern`: A string of the file pattern for TWAS results (default ".alldat").

***Outpus:***

A list containing following dataframes: 
- `cstwas_res`: A dataframe for the CSTWAS results.
	- `Gene` : Feature/gene identifier
	- `Subset_Tissue`: Set of potential gene-expression-specific activated tissues
	- `Number_of_Tissues`: Number of potential gene-expression-specific activated tissues
	- `P_value`: P-value for the cross-tissue subset-based test
- `meta_data`: A dataframe for the tissue-specific TWAS results across multiple tissues.
	- `ID`: Feature/gene identifier
	- `Z`: Z-value transformed from TWAS p-value (qnorm(TWAS.P, lower.tail = F))
	- `TWAS.P`: TWAS p-value
	- `CHR`: Chromosome
	- `BP`: Middle point of gene ((Gene_start + Gene_end)/2)
	- `Tissue`: Tissue of specific gene expression

***Examples:***

```R
res_CSTWAS = run_CSTWAS("path_to_TWAS_resutls", 
                        cov_matrix = "cov_matrix_GTEx_v7", 
                        percent_act_tissue = 0.6, 
                        n_more = 1e+03)
```

- `cstwas_res`:

|  Gene  | Subset_Tissue | Number_of_Tissues | P_value | 
| ------ | ------------- | ----------------- | ---------|
| ALOX5  | Brain_Hippocampus, Small_Intestine_Terminal_Ileum, Spleen, Esophagus_Mucosa | 4 | 0.995 |
| ZDHHC6 | Vagina, Brain_Hypothalamus, Brain_Nucleus_accumbens_basal_ganglia, Brain_Cerebellum, Spleen, Esophagus_Mucosa, Whole_Blood, Cells_Transformed_fibroblasts, Brain_Frontal_Cortex_BA9, Adrenal_Gland | 10 | 0.409 |
| ABCC2  | Brain_Substantia_nigra, Vagina, Cells_EBV.transformed_lymphocytes, Minor_Salivary_Gland, Adipose_Subcutaneous, Heart_Left_Ventricle, Cells_Transformed_fibroblasts, Brain_Amygdala, Testis, Brain_Hippocampus, Artery_Tibial, Whole_Blood | 12 | 0.526 |

- `meta_data`: 

|    ID    |     Z     |  TWAS.P  | CHR |    BP    |       Tissue         |
| -------- | --------- | -------- | --- | -------- | -------------------- | 
| EXOC3L2  | 17.699108 | 2.13e-70 | 19  | 45726674 | Adipose_Subcutaneous |
| PVRL2    | 17.060110 | 1.47e-65 | 19  | 45370958 | Adipose_Subcutaneous |
| CEACAM19 | 8.166297  | 1.59e-16 | 19  | 45178929 | Adipose_Subcutaneous |


### mhp_twas: Manhattan Plot For TWAS Results

***Inputs:***
- `meta_data`: meta_data from run_CSTWAS results.
- `anot_index`: An integer indicating how significant results are to be annotated. (-log10(TWAS.P) > anot_index) This parameter will be ignored if anno_gene is not NULL.
- `ceiling_ctf`: An integer indicating how significant results are to be cut by the ceiling. (-log10(TWAS.P) > ceiling_ctf). If is NULL, it will automatically adjust based on the data.
- `floor_ctf`: An integer indicating how insignificant results are to be cut by the floor (-log10(TWAS.P) < floor_ctf). Default 0.
- `pts_size`: An integer indicating the point size.
- `anno_gene`: A list of genes that need to be annotated.
- `path`: Path for saving the plot.

***Outpus:***

A Manhattan plot of tissue-specific TWAS resutls.

***Examples:***

```R
mhp_twas(res_CSTWAS$meta_data, ceiling_ctf = 30)
```

![Example of Manhattan Plot For TWAS Results](/Plots/example_TWAS_mhp.png)

### mhp_cstwas: Manhattan Plot For the CSTWAS Results

***Inputs:***
- `meta_data`: meta_data from run_CSTWAS results.
- `cstwas_res`: cstwas_res from run_CSTWAS results.
- `anot_index`: An integer indicating how significant results are to be annotated. (-log10(TWAS.P) > anot_index) This parameter will be ignored if anno_gene is not NULL.
- `ceiling_ctf`: An integer indicating how significant results are to be cut by the ceiling. (-log10(TWAS.P) > ceiling_ctf). If is NULL, it will automatically adjust based on the data.
- `floor_ctf`: An integer indicating how insignificant results are to be cut by the floor (-log10(TWAS.P) < floor_ctf). Default 0.
- `pts_size`: An integer indicating the point size.
- `anno_gene`: A list of genes that need to be annotated.
- `path`: Path for saving the plot.

***Outpus:***

A Manhattan plot of the CSTWAS results.

***Examples:***

```R
mhp_cstwas(res_CSTWAS$meta_data, res_CSTWAS$cstwas_res, anot_index = 6)
```

![Example of Manhattan Plot For CSTWAS Results](/Plots/example_csTWAS_mhp.png)


 ### venn_diagram: Venn diagram for significant GReX associations

 ***Inputs:***
 - `meta_data`: meta_data from run_SCTWAS results.
 - `sctwas_res`: sctwas_res from run_SCTWAS results.
 - `merge_range`: An integer indicating how wide (in base pairs) should be considered to merge nearby genes. Default +/- 1000bp..
 - `path`: Path for saving the plot.

 ***Outpus:***

 A Venn diagram showing overlapping conditions of GReX called between Subset-based Cross-tissue TWAS and tissue-specific TWAS.

 ***Examples:***

 ```R
 venn_diagram(res_CSTWAS$meta_data, res_CSTWAS$sctwas_res)
 ```

 ![Example of Venn diagram for significant GReX associations](/Plots/example_Venn_Diagram.png)
