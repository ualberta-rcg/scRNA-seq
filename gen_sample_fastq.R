# ???????????????
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ShortRead", quietly = TRUE)) BiocManager::install("ShortRead")
library(ShortRead)
library(Biostrings)

# Step 1: ?????????TP53??????cDNA?????? (????????????,???????????????)
tp53_cDNA <- DNAString(paste0(
  "AGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGAGGAGCGAAAG",
  "CCTAGCAGCCAGCTGGAAGGAGCCAGGGGGAGCAGAGACCGGCGGAGGGGCTGGGCGGAGGCTCCGAGAGAGCTGAATGAGGCCTT",
  "GGAACTCAAGGATGCCCAGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGC"
))

# Step 2: ????????????
n_reads <- 500   # ???????????????reads
output_prefix <- "tp53_test_S1_L001"  # ?????????????????????

# ????????????????????????
random_seq <- function(length) {
  paste(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
}

# ????????????Illumina???read ID
generate_id <- function() {
  sprintf("A00469:87:H5WY2DRX2:1:1101:%04d:%04d", sample(1000:9999, 1), sample(10000:19999, 1))
}

# Step 3: ??????reads
r1_seqs <- character(n_reads)
r2_seqs <- character(n_reads)
r1_quals <- character(n_reads)
r2_quals <- character(n_reads)
r1_ids <- character(n_reads)
r2_ids <- character(n_reads)

for (i in 1:n_reads) {
  start_pos <- sample(1:(length(tp53_cDNA) - 100), 1)
  r2_fragment <- subseq(tp53_cDNA, start=start_pos, width=sample(90:100, 1))
  
  barcode <- random_seq(16)
  umi <- random_seq(12)
  
  r1_seqs[i] <- paste0(barcode, umi)
  r2_seqs[i] <- as.character(r2_fragment)
  
  r1_quals[i] <- paste(rep("F", 28), collapse = "")
  r2_quals[i] <- paste(rep("F", nchar(r2_fragment)), collapse = "")
  
  id <- generate_id()
  r1_ids[i] <- paste(id, "1:N:0:ATCACG")
  r2_ids[i] <- paste(id, "2:N:0:ATCACG")
}

# Step 4: ??????fastq.gz
r1 <- ShortReadQ(sread = DNAStringSet(r1_seqs),
                 quality = BStringSet(r1_quals),
                 id = BStringSet(r1_ids))
r2 <- ShortReadQ(sread = DNAStringSet(r2_seqs),
                 quality = BStringSet(r2_quals),
                 id = BStringSet(r2_ids))

writeFastq(r1, paste0(output_prefix, "_R1_001.fastq.gz"), compress = TRUE)
writeFastq(r2, paste0(output_prefix, "_R2_001.fastq.gz"), compress = TRUE)

cat("??? Done! Generated files:\n",
    paste0(output_prefix, "_R1_001.fastq.gz"), "\n",
    paste0(output_prefix, "_R2_001.fastq.gz"), "\n")
