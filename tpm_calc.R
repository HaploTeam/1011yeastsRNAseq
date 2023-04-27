#! /usr/bin/env Rscript

# Compute TPM from raw counts
# 
# arg1 : table of counts

# Transcript Sample1 Sample2 Sample N
# T1 15 30 50       
# T2 0 0 10       
# T3 100 50 100 


# arg2: list of transcript lengths
#
# Transc1 1500
# Transc2 300
# Transc3 2000

args = commandArgs(trailingOnly=TRUE)

# Read raw counts
#raw_counts = read.table("/Users/fabien/counts_table_duplsubset_ORFsubset.txt", header=T, check.names=F, stringsAsFactors=F, row.names=1)
raw_counts = read.table(args[1], header=T, check.names=F, stringsAsFactors=F, row.names=1)

# Read transcript lengths and effective lengths
#lengths = read.table("/Users/fabien/RNAseq/transcripts_lengths_SGDannot_list_withpangenomeORFs.txt", header=F, stringsAsFactors=F, row.names=1)
lengths = read.table(args[2], header=F, stringsAsFactors=F, row.names=1)
effective_length = lengths[rownames(raw_counts), 1] - 75 + 1

# Min transcript length set to opt$readlen to avoid negative length
effective_length = ifelse(effective_length > 75, effective_length, 75)

# Compute TPM
tpm = apply(raw_counts, 2, function(v) v/effective_length / sum(v/effective_length) * 1e6)

# Output
write.table(tpm, sep="\t", quote=F, file="tpm_table.txt")
