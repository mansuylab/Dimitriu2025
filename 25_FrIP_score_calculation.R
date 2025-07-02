library(Rsamtools)
library(chromVAR)
library(dplyr)

bamFolder <- "/mnt/groupMansuy/maria/NanoTag_07-2023/data/aligned/deduplicated"
bam_files <- list.files(bamFolder, pattern = "\\.bam$", full.names = TRUE)


# Select only TAF3 files
target <- c("TAF", "Nanog", "CTCF")
bam_files <- bam_files[grep(target[1], bam_files)]

# Use the TAF3 ChIP-seq peaks to calculate FrIP score
peaksFile <- "/mnt/groupMansuy/maria/dataToCompare/literature_ChIP/peaks/esc_h3k4me3_chip.bed"
peaks <- rtracklayer::import(peaksFile, format = "narrowPeak")

data <- data.frame(sample = character(length(bam_files)),
                   reads = character(length(bam_files)),
                   peaks = character(length(bam_files)),
                   number_of_reads = integer(length(bam_files)),
                   paired_or_not = logical(length(bam_files)),
                   number_of_reads_in_peaks = integer(length(bam_files)),
                   frip_score = numeric(length(bam_files)))

# Process BAM files and calculate statistics
for (i in 1:length(bam_files)) {
  data$sample[i] <- bam_files[i]
  data$reads[i] <- gsub(bamFolder, "", bam_files[i])
  reads <- bam_files[i]
  bam <- BamFile(reads)
  data$number_of_reads[i] <- countBam(bam)$records
  data$paired_or_not[i] <- !grepl("ENCODE", data$reads[i])
  fragment_counts <- getCounts(reads, peaks, paired = data$paired_or_not[i], by_rg = FALSE, format = "bam")
  inPeakN <- counts(fragment_counts)[, 1] %>% sum
  data$number_of_reads_in_peaks[i] <- inPeakN
}

# Calculate and round the FRiP score
for (i in 1:nrow(data)) {
  data$frip_score[i] <- data$number_of_reads_in_peaks[i] / data$number_of_reads[i]
}

data$frip_score <- round(data$frip_score, 2)
