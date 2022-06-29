# SCTWAS
Transcriptome-wide association study (TWAS) is introduced to identify 
    significant expression-trait associations through imputations. It has been widely 
    used to analyze tissue-specific associations with the reference expression 
    quantitative trait loci (eQTL) panel. To increase the statistical power of TWAS results, 
    meta-analysis methods aggregating TWAS results across multiple tissues are developed. 
    However, most existing meta-analysis methods lose interpretation of disease etiology and 
    have limited power to identify weaker associations when only a few tissues are weakly 
    activated. Therefore, we developed the subset-based cross-tissue meta-analysis method, 
    also called subset-based cross-tissue transcriptome-wide association study (SCTWAS). 
    In this package, we aggregate the TWAS results across tissues and perform meta-analysis 
    through the subset-based test. R functions are provided for researchers to integrate 
    TWAS results across multiple tissues and visualize the result.

# Installation

Use the following codes to install the ADTools package
```R
library(devtools)
install_github("Thewhey-Brian/SCTWAS")
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
  
  - `ID`: Feature/gene identifier
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
{
   echo ""
   echo "Usage: $0 -d Disease -t Tissue list"
   echo -e "\t-d Disease"
   echo -e "\t-t Tissue list"
   exit 1 # Exit script after printing help
}

while getopts "d:t:" opt
do
   case "$opt" in
      d ) dise="$OPTARG" ;;
      t ) tiss="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$dise" ] || [ -z "$tiss" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

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
helpFunction()
{
   echo ""
   echo "Usage: $0 -o Output direction -f File of tissue list"
   echo -e "\t-o Output direction"
   echo -e "\t-f File of tissue list"
   exit 1 # Exit script after printing help
}

while getopts "o:f:" opt
do
   case "$opt" in
      o ) out_dir="$OPTARG" ;;
      f ) file="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$out_dir" ] || [ -z "$file" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
tissues=$(cat $file)
for tissue in $tissues
do
	head -n 1 ${out_dir}/${tissue}_1.dat | awk '{print $1, $3, $4, $5, $6, $19, $20}' > ${out_dir}/all.${tissue}.alldat
	tail -q -n +2 ${out_dir}/${tissue}*.dat | awk '{print $1, $3, $4, $5, $6, $19, $20}' >> ${out_dir}/all.${tissue}.alldat
	echo "Finish ${tissue}, length: $(wc -l ${out_dir}/all.${tissue}.alldat)"
done
```

## Run SCTWAS

With the formatted TWAS results, we are ready to perform the subset-based cross-tissue TWAS anslysis.

### run_SCTWAS: Run Subset-based Cross-tissue TWAS

***Inputs:***
- `path`: A string of the direction for TWAS results.
- `cov_matrix`: A string indicating list of matrix of the gene expression covariance matrix across tissues from the reference panel (default using "cov_matrix_GRCh37", can also change to "cov_matrix_GRCh38" or use your own matrix list). This parameter is omitted if cov_matrix_path is specified.
- `cov_matrix_path`: Path for downloaded reference gene expression covariance matrix across tissues (need to be named as "cov_matrix") (the reference matrix can be downloaded from: https://github.com/Thewhey-Brian/SCTWAS) If NULL, the function will automatically download the reference panel indicated by cov_matrix from the GitHub repository.
- `percent_act_tissue`: A decimal of the minimum percent of activated tissues for each gene regulated expression.
- `n_more`: Simulation times for small p-values (default 1e+04; Caution: a very large number may lead to long calculation time; a very small number may lead to inaccurate p-value estimation).
- `gene_list`: An array of the list of interested genes (default NULL; if NULL, it will go over all genes in the TWAS results; if not NULL, percent_act_tissue will be ignored).
- `pattern`: A string of the file pattern for TWAS results (default ".alldat").

***Outpus:***
- `sctwas_res`: A dataframe for the Subset-based Cross-tissue TWAS results.
- `meta_data`: A dataframe for the tissue-specific TWAS results across multiple tissues.

***Examples:***

```R
res_SCTWAS = run_SCTWAS("path_to_TWAS_resutls", 
                        cov_matrix = "cov_matrix_GRCh37", 
                        percent_act_tissue = 0.6, 
                        n_more = 1e+03)
```

- `sctwas_res`:

| Gene | 

- `meta_data`: 


### mhp_twas: Manhattan Plot For TWAS Results



### mhp_sctwas: Manhattan Plot For the Subset-based Cross-tissue TWAS Results






