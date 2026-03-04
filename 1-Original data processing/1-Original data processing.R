###########dada2_Sunjieyu############
setwd("/home/omics_2023060641/SS")
rm(list = ls())
library(dada2)
library(Rcpp)

# View files in the directory
cat("Files in the current directory：\n")
list.files(path = "/home/omics_2023060641/SS", pattern = "\\.fastq$", all.files = FALSE, full.names = FALSE, recursive = FALSE)

# Returns the full path of the sequencing file
fnFs <- sort(list.files(path = "/home/omics_2023060641/SS", pattern = "\\.fastq$", full.names = TRUE))

# Check the number of files
cat("\nFastq number of files found：", length(fnFs), "\n")
cat("The first few documents：\n")
print(head(fnFs))

# Extract the sample name
Sample.names <- sapply(strsplit(basename(fnFs), "\\.fastq$"), `[`, 1)

print(Sample.names)

# sequence file quality detection, drawing the quality diagram of the first two samples
plotQualityProfile(fnFs[1:2])

# Sequence filtering and clipping
# Create Filtered subdirectory
filtered_dir <- "/home/omics_2023060641/SS/Filtered"
if(!dir.exists(filtered_dir)) {
  dir.create(filtered_dir)
}

# Set the output file name
filtFs <- file.path(path = filtered_dir, paste0(Sample.names, "_filt.fastq"))

print(head(filtFs))

# Filter file output
filt.result <- filterAndTrim(fnFs, filtFs, 
                             truncQ = 2, 
                             truncLen = 0,  
                             trimLeft = 0, 
                             minLen = 100,  
                             maxN = 0,      
                             minQ = 0, 
                             maxEE = 2,     
                             rm.phix = TRUE, 
                             compress = FALSE,  
                             multithread = TRUE)

print(filt.result)

# Calculation of Filter Retention Rate
retention_rate <- filt.result[,2] / filt.result[,1] * 100
print(data.frame(Sample = Sample.names, 
                 Input = filt.result[,1], 
                 Output = filt.result[,2], 
                 Retention_rate = paste0(round(retention_rate, 1), "%")))


errFs <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errFs, nominalQ = TRUE)

# Removing repetitive sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)

names(derepFs) <- Sample.names

# Further quality control based on error model
dadaFs <- dada(derepFs, err = errFs, multithread = TRUE, verbose = TRUE)

print(dadaFs[[1]])


plotErrors(dadaFs[[1]])


ASV.table <- makeSequenceTable(dadaFs)

print(dim(ASV.table))

seqlen <- nchar(colnames(ASV.table))
print(table(seqlen))

# Remove chimera
ASV.table.nochim <- removeBimeraDenovo(ASV.table, method = "consensus", multithread = TRUE, verbose = TRUE)


cat("\nComparison before and after removal of chimeras：\n")
cat("Number of ASVs before removal：", ncol(ASV.table), "\n")
cat("The number of ASV after removal：", ncol(ASV.table.nochim), "\n")
cat("Chimera proportion：", round(1 - ncol(ASV.table.nochim)/ncol(ASV.table), 4)*100, "%\n")


cat("\nThe final ASV table dimension：\n")
print(dim(ASV.table.nochim))


# write.csv(t(ASV.table.nochim), file = "ASV_table.csv")


dim(ASV.table)

# Inspect distribution of sequence lengths
table(nchar(getSequences(ASV.table)))


ASV.table.nochim <- removeBimeraDenovo(ASV.table, method = "consensus", multithread = TRUE, verbose = TRUE)
table(nchar(getSequences(ASV.table.nochim)))

# chimera result statistics, the dim ( ) function outputs the number of rows of a matrix, followed by the number of columns of the output matrix.
dim(ASV.table.nochim)
sum(ASV.table.nochim)/sum(ASV.table)


GetN <- function(x) sum(getUniques(x))

# Merge the data volume of each sample step by step
track <- cbind(filt.result, sapply(dadaFs, GetN), sapply(derepFs, GetN), rowSums(ASV.table.nochim), round((rowSums(ASV.table.nochim)/sapply(derepFs, GetN))*100, 2), round((rowSums(ASV.table.nochim)/filt.result[,1])*100, 2))


colnames(track) <- c("Input", "Filtered", "DenoisedF",  "derepFs", "Nonchim", "Total reads lost", "ASVs")


rownames(track) <- Sample.names


head(track, n = 10)


# Sequence species annotation
Taxa <- assignTaxonomy(ASV.table.nochim, paste0(path = "/home/omics/database/dada2_db/silva_nr99_v138_train_set.fa.gz"), multithread = TRUE)


Taxa <- addSpecies(Taxa, paste0(path = "/home/omics/database/dada2_db/silva_species_assignment_v138.fa.gz"), allowMultiple = TRUE)

# Save the species annotation information, remove the sequence name, and only display the species information
Taxa.print <- Taxa
rownames(Taxa.print) <- NULL
head(Taxa.print)


setwd("/home/omics_2023060641/SS")

# # to obtain the ASV sequence.
Seqs <- getSequences(ASV.table.nochim) 

ASV_headers <- vector(dim(ASV.table.nochim)[2], mode = "character") 

for (i in 1:dim(ASV.table.nochim)[2]) 
{
  ASV_headers[i] <- paste("> ASV", i, sep = "_")
}

ASV_fasta <- c(rbind(ASV_headers, Seqs))

# # Save annotation information and ASV _ table
ASV_ID <- vector(dim(ASV.table.nochim)[2], mode = "character") 

for (i in 1:dim(ASV.table.nochim)[2]) 
{
  ASV_ID[i] <- paste("ASV", i, sep = "_")
}

Taxonomy <- cbind('#seq' = rownames(Taxa), ASV_ID, Taxa)

ASV_table <- cbind('#seq' = rownames(Taxa), ASV_ID, t(ASV.table.nochim))

ASV_table.taxonomy <- cbind('#seq' = rownames(Taxa), ASV_ID, t(ASV.table.nochim), Taxa)

Track <- cbind('#SampleID' = rownames(track), track)


# # Results derived
write.table(ASV_table, "ASV.counts.txt", sep = "\t", quote = F, row.names = F)

write.table(Taxonomy, "ASV.Taxonomy.txt", sep = "\t", quote = F, row.names = F)

write.table(ASV_table.taxonomy, "ASV.counts_taxon_species.txt", sep = "\t", quote = F, row.names = F)

write.table(Track, "ASV.track.txt", sep = "\t", quote = F, row.names = F)

write(ASV_fasta, "ASVs.fa")



