#####
#
#1) Code to install and read in packages
#
#BiocManager::install("XXXX")
####



#####
#
#2) Setting up input PATH
#Define input path
path <- "XXXXX/" # CHANGE THIS to the directory containing the fastq files after unzipping.
#
#make sure you see a set of 24 files
list.files(path)
#
#####

#####
#
#3) gathering file lists
#
# Forward and reverse fastq filenames have format: SRRXXXX_1.fastq and SRRXXX_1.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SRRXXXX_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
##

#####
#
#4) Quality profiles
#
#plotting the quality scores of the first two forward reads
plotQualityProfile(fnFs[1:2])
#plotting the quality scores of the first two forward reads
plotQualityProfile(fnRs[1:2])

#plotting the quality scores of the 3rd and 4th forward reads
plotQualityProfile(fnFs[XXX:XXX])

#####

#####
#
#5) filtering data
#
#
#First set up new filenames
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# set maxEE to 2 as it is a better measure than average quality score
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,XXXX),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=8, verbose = TRUE) # 


# calculate percentage passing filter and trim (out_perc)
out_perc = as.data.frame(out)
out_perc$perc_pass = out_perc$reads.out/out_perc$reads.in*100

#####
#
#6) Learning the error models
#
#####

# Make error models 
errF <- learnErrors(filtFs, multithread=8, verbose = TRUE, nbases = 1e9, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=8, verbose = TRUE, nbases = 1e9, MAX_CONSIST = 20)

# visualize the estimated error rates
#The error rates for each possible transition (A→C, A→G, …) are shown.
# Points are the observed error rates for each consensus quality score.
#The black line shows the estimated error rates after convergence of the machine-learning algorithm.
#The red line shows the error rates expected under the nominal definition of the Q-score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#####
#
#7) Run the denoising algorithm
#
#####

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#####
#
#8) merge two pairs of reads
#
#####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
hist((nchar(getSequences(seqtab))))

#####
#
#9) Remove chimeras/bimeras
#
#####

#identify chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#The fraction of chimeras 
sum(seqtab.nochim)/sum(seqtab)

#####
#
#10) statistic of the preprocessing
#
#####

#Fancy stuff to get a nice table
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

track

track/6250



#####
#
#11) Download classification database
#
#####

download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
              dest="silva_nr99_v138.1_train_set.fa.gz")

#####
#
#12) Classify reads
#
#####

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#taxa2 <- addSpecies(taxa, "silva_nr99_v138.1_train_set.fa.gz")

saveRDS(object = taxa, file = "taxa.Rdata")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)
