#!/bin/R

library(data.table)
library(DESeq2)

## DEFINE SEX #######################################################################################################
.sex = 'all' # should be "all", "male", or "female"
#####################################################################################################################

# get release tissue code 
library(MotrpacBicQC)
bic_animal_tissue_code = data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code[grepl('Powder',bic_tissue_name), my_tissue := gsub(' powder','',tolower(bic_tissue_name))]
bic_animal_tissue_code[motrpac_tissue_code == 'T30 Rat PaxGene RNA', my_tissue := 'paxgene rna']
tissue_codes = bic_animal_tissue_code[!is.na(my_tissue)]

## SPECIFY A TISSUE FROM ONE OF THE FOLLOWING STRINGS ###############################################################
print(tissue_codes[,my_tissue])
.tissue = 'gastrocnemius'
#####################################################################################################################

tissue_codes[my_tissue == .tissue]
tissue_release_name = tissue_codes[my_tissue == .tissue, tissue_name_release]
tissue_release_name

# download and read in a file from gsutil 
dl_read_gcp <- function(path,sep='\t',tmpdir='/tmp/motr',GSUTIL_PATH='~/google-cloud-sdk/bin/gsutil'){
  # GSUTIL_PATH needs to point to the software location on your computer
  # tmpdir needs to point to a temporary data path where files from Google Cloud can be downloaded locally and read in 
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
  system(cmd,ignore.stdout = T,ignore.stderr = T)
  # read in the data
  new_path = sprintf('%s/%s',tmpdir,basename(path))
  dt <- fread(new_path,sep=sep,header=T)
  return(dt)
}

dmaqc_metadata = 'gs://motrpac-internal-release2-results/pass1b-06/phenotype/motrpac_pass1b-06_pheno_viallabel-data.txt'
dmaqc_dict = 'gs://motrpac-internal-release2-results/pass1b-06/phenotype/motrpac_pass1b_pheno_merged-dictionary.txt'

# download and format phenotypic data 
dmaqc_metadata = dl_read_gcp(dmaqc_metadata)
cols = dl_read_gcp(dmaqc_dict)
old_cols = colnames(dmaqc_metadata)
new_cols = tolower(cols[match(old_cols, BICUniqueID), FullName])
colnames(dmaqc_metadata) = new_cols

# make some variables human-readable
for (var in c('key.protocol','key.agegroup','key.intervention','key.sacrificetime','registration.sex')){
  d = cols[Field.Name == gsub('.*\\.','',var)]
  keys=unname(unlist(strsplit(d[,Categorical.Values],'\\|')))
  values=unname(unlist(strsplit(d[,Categorical.Definitions],'\\|')))
  names(values) = keys
  # match keys to values; create new column 
  new_var = gsub(".*\\.","",var)
  dmaqc_metadata[,(new_var) := unname(values)[match(get(var), names(values))]]
}

# format to be compatible with existing functions
dmaqc_metadata[,sacrificetime := sapply(sacrificetime, function(x) gsub(' week.*','w',x))]

# clean up 'intervention'
dmaqc_metadata[intervention == 'Control', intervention := 'control']
dmaqc_metadata[grepl('Training',intervention), intervention := 'trained']

# make "group" - "1w", "2w", "4w", "8w", "control8w"
dmaqc_metadata[,group := sacrificetime]
dmaqc_metadata[intervention == 'control', group := paste0('control',group)]

# load counts 
counts = dl_read_gcp(sprintf('gs://motrpac-internal-release2-results/pass1b-06/transcriptomics/%s/transcript-rna-seq/results/motrpac_pass1b-06_%s_transcript-rna-seq_rsem-genes-count.txt', tissue_release_name, tissue_release_name))
counts = data.frame(counts)
rownames(counts) = counts$gene_id
counts$gene_id = NULL

# filter by genes in normalized data 
tmm = dl_read_gcp(sprintf('gs://motrpac-internal-release2-results/pass1b-06/transcriptomics/%s/transcript-rna-seq/results/motrpac_pass1b-06_%s_transcript-rna-seq_normalized-log-cpm.txt', tissue_release_name, tissue_release_name))
# remove pid and bid rows
tmm = tmm[3:nrow(tmm)]
counts = counts[tmm[,viallabel],]

# subset metadata 
meta_data = dmaqc_metadata[as.character(viallabel) %in% colnames(tmm)]
# remove reference standards from counts 
counts = counts[,paste0('X', meta_data[,viallabel])]
dim(counts)

# subset by sex
if(.sex != 'all'){
  meta_data = meta_data[tolower(sex) == .sex]
  counts = counts[,paste0('X', meta_data[,viallabel])]
}

# coerce counts to integer (RSEM has fractional values for counts sometimes)
counts_round = data.frame(apply(counts, c(1,2), as.integer)) 

# run DESeq

if(.sex == 'all'){
  contrast = '~ sex + group'
}else if(.sex %in% c('male','female')){
  contrast = '~ group'
}

# set it up - no interaction term 
dds = DESeqDataSetFromMatrix(countData = counts_round,
                            colData = meta_data,
                            design = eval(parse(text=contrast)))
dds = DESeq(dds)
w1 = results(dds, c('group','1w','control8w'))
