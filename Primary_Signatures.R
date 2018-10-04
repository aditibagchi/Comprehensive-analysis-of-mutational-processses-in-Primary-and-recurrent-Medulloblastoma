##Primary medulloblastoma from ICGC

Germline_Primary <- read.csv("~/Desktop/Data_Analysis/Data_Analysis_July_-December_2018/July_december_2018/Germline_Primary_137_final_Mut_sigs.csv")
View(Germline_Primary)
##deconstructSigs
##BSgenome.Hsapiens.UCSC.hg19
Germline_Primary_sigs <- mut.to.sigs.input(Germline_Primary, sample.id = "Sample", chr = "chr", pos = "X",
                  ref = "ref", alt = "alt", bsg = NULL)
View(Germline_Primary_sigs)

PrimaryMB02 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-02",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB02)
makePie(PrimaryMB02, sub = "PrimaryMB02", add.color = NULL)
PrimaryMB03 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-03",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB03)
makePie(PrimaryMB03, sub = "PrimaryMB03", add.color = NULL)
PrimaryMB04 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-04",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB04)
makePie(PrimaryMB04, sub = "PrimaryMB04", add.color = NULL)
makePie(PrimaryMB03, sub = "PrimaryMB03", add.color = NULL)


PrimaryMB05 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-05",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB05)
makePie(PrimaryMB05, sub = "PrimaryMB05", add.color = NULL)

PrimaryMB06 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-06",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB06)
makePie(PrimaryMB06, sub = "PrimaryMB06", add.color = NULL)


PrimaryMB07 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-07",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB07)
makePie(PrimaryMB07, sub = "PrimaryMB07", add.color = NULL)

PrimaryMB08 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-08",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB08)
makePie(PrimaryMB08, sub = "PrimaryMB08", add.color = NULL)
PrimaryMB09 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-09",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB09)
makePie(PrimaryMB09, sub = "PrimaryMB09", add.color = NULL)

PrimaryMB10 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-10",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB10)
makePie(PrimaryMB10, sub = "PrimaryMB10", add.color = NULL)


PrimaryMB11 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-11",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB11)
makePie(PrimaryMB11, sub = "PrimaryMB11", add.color = NULL)

PrimaryMB12 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-12",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB12)
makePie(PrimaryMB12, sub = "PrimaryMB12", add.color = NULL)


PrimaryMB13 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-13",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB13)
makePie(PrimaryMB13, sub = "PrimaryMB13", add.color = NULL)

PrimaryMB14 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-14",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB14)
makePie(PrimaryMB14, sub = "PrimaryMB14", add.color = NULL)

PrimaryMB15 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-15",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB15)
makePie(PrimaryMB15, sub = "PrimaryMB15", add.color = NULL)

PrimaryMB16 <- whichSignatures(tumor.ref = Germline_Primary_sigs, sample.id = "MB-REC-16",
                               signatures.ref = signatures.nature2013, associated = c(),
                               signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                               tri.counts.method = "default")
plotSignatures(PrimaryMB16)
makePie(PrimaryMB16, sub = "PrimaryMB16", add.color = NULL)

##########

Packages Used
##MutationalPatterns
##BSgenome.Hsapiens.UCSC.hg19

Germline_Primary_sigs_t  <- t(Germline_Primary_sigs)
View(Germline_Primary_sigs_t)
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(Germline_Primary_sigs_t), cancer_signatures$Somatic.Mutation.Type)
#Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
cos_sim(Germline_Primary_sigs_t[,1], cancer_signatures[,1])
cos_sim_samples_signatures = cos_sim_matrix(Germline_Primary_sigs_t, cancer_signatures)
# Plot heatmap with specified signature order
plot_cosine_heatmap(cos_sim_samples_signatures,cluster_rows = TRUE)
##Fit Mutation matrix to the cosmic mutational signature
fit_res <- fit_to_signatures(Germline_Primary_sigs_t, cancer_signatures)
# Select signatures with some contribution
select <- which(rowSums(fit_res$contribution) > 10)
# Plot contribution barplot
plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = FALSE, mode = "absolute")
