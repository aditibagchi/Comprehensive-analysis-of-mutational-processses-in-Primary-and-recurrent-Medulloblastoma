######
Work Flow working genrating a final Sigs_Input file for recurrent medulloblastoma cases.  
#Packages used
maftools::
deconstructSigs::  
BSgenome.Hsapiens.UCSC.hg19
BSgenome()

### 
#Cases : MB-REC-02, 03,41,06,09,10 from ICGC
#Variants called against germline. Recurrence at mestatic site. 
Germline_Met_Recur<- read.delim("~/Desktop/Data_Analysis/last MAF files/ICGC_germline_metastatic_reduced_maf_final.txt", header=TRUE)
View(Germline_Met_Recur)
Germline_Met_Recur_Maf <- read.maf(maf = "~/Desktop/Data_Analysis/last MAF files/ICGC_germline_metastatic_reduced_maf_final.txt")
mafSummary(Germline_Met_Recur_Maf)
###
##$variants.per.sample
##Tumor_Sample_Barcode Variants
1:            MB-REC-06   511043 
2:            MB-REC-09     6293
3:            MB-REC-10     4690
4:            MB-REC-02     3728
5:            MB-REC-41     2181
6:            MB-REC-03     1419
#Cases : MB-REC-16, 14,11,07,08,04,13,12,15 from ICGC
#Variants called against germline. Recurrence at Primary site. 
Germline_Recur<- read.delim("~/Desktop/Data_Analysis/last MAF files/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt", header=FALSE, comment.char="#")
View(Germline_Recur)
Germline_Recur_Maf <- read.maf(maf = "~/Desktop/Data_Analysis/last MAF files/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt" )
mafSummary(Germline_Recur_Maf)
###
##$variant.per.sampple
##Tumor_Sample_Barcode Variants
1:            MB-REC-16    28931
2:            MB-REC-14    11212
3:            MB-REC-11     9789
4:            MB-REC-07     7300
5:            MB-REC-08     6913
6:            MB-REC-04     6250
7:            MB-REC-13     5689
8:            MB-REC-12     5332
9:            MB-REC-15     4453
###
## To create a DF with rows with  number of cases (Recurrent) and 96 columns of of the types of single nucleotide subsitution.
##Subset_Met_GermLine_73.csv was created earlier after deleting the 73 variants from the original MAF file that didnot match the context.
Germline_Met_Recur_Sigs_Input_file <- read.csv("~/Desktop/Data_Analysis/last MAF files/SUbset_Mets_Germline_73.csv")
Germline_Met_Recur_Sigs <- mut.to.sigs.input(Germline_Met_Recur_Sigs_Input_file, sample.id = "Sample", chr = "chr", pos = "X",
                                           ref = "ref", alt = "alt", bsg = NULL)
View(Germline_Met_Recur_Sigs)

##To create a DF with rows with  number of cases (Recurrent) and 96 columns of of the types of single nucleotide subsitution
## This will be done on the 9 cases in Germline_Recur_Maf
#Create a CSV File with "Sample", "chr", 'pos", "ref", "alt" ; a subset of the maffile.
Germline_Recur_Maf_Subset <- read.csv("~/Desktop/Data_Analysis/last MAF files/Germline_Recur_Maf_Subset.csv")
Germline_Recur_Maf_Subset <- as.data.frame(Germline_Recur_Maf_Subset)
Germline_Met_Recur_Sigs <- mut.to.sigs.input(Germline_Recur_Maf_Subset, sample.id = "sample", chr = "chr", pos = "pos",ref = "ref", alt = "alt", bsg = NULL)
View(Germline_Met_Recur_Sigs)
##The DF
             A[C>A]A A[C>A]C A[C>A]G A[C>A]T C[C>A]A C[C>A]C C[C>A]G C[C>A]T G[C>A]A G[C>A]C G[C>A]G G[C>A]T T[C>A]A T[C>A]C T[C>A]G
MB-REC-04     115      73       5      36      90      51      19      81      83      43      10      30      90      57       4
MB-REC-07     134      59       3      42      78      59      16      64     122      89      17      88      88      43      12
MB-REC-08      59      32       8      40      51      45       5      53      43      39       6      23      39      34      10
MB-REC-11     192     154      27     123     115     136      21     155      70      82      13      60     141     143      10
MB-REC-12      95      56      11      50      92      60      12     105      38      38      10      37      76      81       5
MB-REC-13      99      74       8      62     122      81      11     146     112      64      15      73      80      97       8
MB-REC-14     195     188      33     138     236     143      20     246     120     127      19      86     145     153      16
MB-REC-15      69      64      11      41      61      56      20      71      49      61      17      28      70      55      21
MB-REC-16     240     137      33     198     235     184      56     268     173     117      29     133     236     129      22
             T[C>A]T A[C>G]A A[C>G]C A[C>G]G A[C>G]T C[C>G]A C[C>G]C C[C>G]G C[C>G]T G[C>G]A G[C>G]C G[C>G]G G[C>G]T T[C>G]A T[C>G]C
MB-REC-04     102      40      22       9      24      22      15      12      28      29      25       5      18      32      18
MB-REC-07     109      46      19      11      36      28      31      10      27      54      38      22      29      29      27
MB-REC-08      54      48      41       8      47      24      23      10      20      27      28       8      27      36      40
MB-REC-11     173      57      79      21      51      24      35      10      30      22      37       7      35      29      56
MB-REC-12      73      32      53       8      44      19      28       7      23      23      42       4      38      28      47
MB-REC-13     118      32      23      13      42      38      32       6      41      41      38      17      25      39      33
MB-REC-14     212      75      69      20      71      39      36      13      62      34      33      11      50      60      63
MB-REC-15      85      29      19      13      26       9      18       6      20      15      11       8      12      23      21
MB-REC-16     231     115      45      26     103      48      47      16      38      88      53      16      63     133      60
             T[C>G]G T[C>G]T A[C>T]A A[C>T]C A[C>T]G A[C>T]T C[C>T]A C[C>T]C C[C>T]G C[C>T]T G[C>T]A G[C>T]C G[C>T]G G[C>T]T T[C>T]A
MB-REC-04       7      43      96      45     270      74      79      73     132      81      92      73     206      64      68
MB-REC-07       2      37     123      59     120      86      81      80      67      94      90      66      91     103      77
MB-REC-08      11      51     108      55     105      90      86      68      90      68      60      55     103      59      65
MB-REC-11       7      79     239     134     494     151     162     205     312     207     127     152     357     133     165
MB-REC-12       5      84     104      77      93      86      68     127      68     119      61     104      56      92      65
MB-REC-13       9      46     105      83     101     124      62     117      85     124     113     128     122     126      79
MB-REC-14       8     120     229     153     437     175     168     219     238     208     107     149     252     127     133
MB-REC-15       8      21      81      59     334      66      70     156     210     111      65      73     257      54      48
MB-REC-16      12     132     219     119     130     263     152     112      80     234     163     101     111     182     201
              T[C>T]C T[C>T]G T[C>T]T A[T>A]A A[T>A]C A[T>A]G A[T>A]T C[T>A]A C[T>A]C C[T>A]G C[T>A]T G[T>A]A G[T>A]C G[T>A]G G[T>A]T
MB-REC-04      68     117      79      37      22      24      49      17      26      31      37      23       8      17      16
MB-REC-07      76      42     114      33      20      27      30      28      36      37      42      15      15      18      26
MB-REC-08      80      52      87      23      31      18      24      22      29      24      30       9      19      26      20
MB-REC-11     193     201     164      57      57      65      92      34      66      74      82      32      29      37      42
MB-REC-12     126      40     144      48      25      44      53      33      44      39      56      27      20      20      23
MB-REC-13     107      50     136      36      31      41      68      51      63      68      86      39      22      40      35
MB-REC-14     206     158     192      78      66      78     126     104     121     125     169      46      54      58      60
MB-REC-15      73     125      77      19      17      19      44      20      25      22      38      13       8      17      12
MB-REC-16     164      56     289     153      70     103     157     102      74     102     155      80      31      63      84
              T[T>A]A T[T>A]C T[T>A]G T[T>A]T A[T>C]A A[T>C]C A[T>C]G A[T>C]T C[T>C]A C[T>C]C C[T>C]G C[T>C]T G[T>C]A G[T>C]C G[T>C]G
MB-REC-04      64      34      14      47      88      51      70      75      36      55      53      49      44      40      53
MB-REC-07      39      18      22      55      89      47      53      57      35      59      50      49      59      47      50
MB-REC-08      23      13      12      25      89      36      69      66      41      40      64      46      43      42      44
MB-REC-11      62      56      48     115     105      63      90     124      34      68      56      64      59      44      66
MB-REC-12      41      24      23      68      88      38      43      77      30      36      41      58      43      34      28
MB-REC-13      65      36      26      83      83      44      39      84      54      51      33      53      72      53      44
MB-REC-14     110      69      66     153     160      63     115     159      77     106      91     108      82      64      92
MB-REC-15      46      22      15      37      49      23      39      58      26      21      22      38      25      35      21
MB-REC-16     156      87      84     207     123      57      52     119      80      62      50     104      70      50      41
              G[T>C]T T[T>C]A T[T>C]C T[T>C]G T[T>C]T A[T>G]A A[T>G]C A[T>G]G A[T>G]T C[T>G]A C[T>G]C C[T>G]G C[T>G]T G[T>G]A G[T>G]C
MB-REC-04      59      64      60      44      62      21      12      20      35      11      12      20      13     145      70
MB-REC-07      46      53      56      49      85      27      13      18      18      12       6      21      16      33      18
MB-REC-08      62      44      53      49      49      23      14      28      32       5      16      15      27      17      20
MB-REC-11      75      52      61      49      83      20      25      49      49      21      22      29      27      42      35
MB-REC-12      48      54      36      29      68      25      26      35      52      12      16      21      20      33      28
MB-REC-13      92      52      45      30     105      16       9      16       4      19      10      28      16      14      13
MB-REC-14     106      98      99      73     133      55      41      58      90      19      29      34      26      31      36
MB-REC-15      33      27      40      19      45      16       5      17      18       8       7       7      13       5      13
MB-REC-16      72     112      86      41     148      87      42      65      91      27      20      29      25      53      23
              G[T>G]G G[T>G]T T[T>G]A T[T>G]C T[T>G]G T[T>G]T
MB-REC-04     411     129      25      15      21      45
MB-REC-07      96      51      24      23      38      52
MB-REC-08      51      49      17      24      21      51
MB-REC-11     181      62      19      32      39      81
MB-REC-12     110      52      26      20      24      40
MB-REC-13      59      24      25       5      24      33
MB-REC-14      78      57      32      38      47      74
MB-REC-15      29      10      16      11      12      31
MB-REC-16      66      36      80      54      45     109
##To create a DF with rows with  number of cases (Recurrent) and 96 columns of of the types of single nucleotide subsitution
## This will be done on the  cases in Primary_Recurrent_ICGCDF_96.1.csv , these are recurrent cases in data base with no matched primary but matched germline and also recurrent cases with matched primary and no germline control. The variants for cases with no germline control were called against primary medullblastoma. 

Recurrent_Med<- read.csv("~/Desktop/Data_Analysis/Data_analysis_June _2018/Initial_plan_for_analysis/Primary_recurrent_ICGC_DF_96_1.csv")
######
##Final DF with all recurrent cases with rows as number of cases and 96 columns 
Final_recurrant_DF <- read.csv("~/Desktop/Data_Analysis/last MAF files/DF_Sigs_Input_Recurrent.csv")
Final_recurrant_DF <- as.data.frame(Final_recurrant_DF)
View(Final_recurrant_DF)

Packages Used
##MutationalPatterns
##BSgenome.Hsapiens.UCSC.hg19
Final_recurrent_DF1 <- read.csv("~/Desktop/Final_recurrent_DF.csv")
Final_recurrent_DF1 <- as.data.frame(Final_recurrent_DF1)
View(Final_recurrent_DF1)
Final_recurrent_DF1_t <- t(Final_recurrent_DF1)
View(Final_recurrent_DF1_t)
Germline_Met_Recur_Sigs_t <- t(Germline_Met_Recur_Sigs)

Packages Used
##MutationalPatterns
##BSgenome.Hsapiens.UCSC.hg19

sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(DF_sigs_input_recurrent_t), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])