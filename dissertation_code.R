# load in data 
load("ASD_CB_data.rda")
# load in phenotypes 
read.csv("CBpheno.csv", header=T) -> pheno
# slice data columns, remove unwanted columns 
pheno[ ,c('Sample.1', 'Sentix', 'Diagnosis', 'Sex', 'Age.Yrs')] -> phenosample
# slice data rows, remove unwanted rows 
phenosample[!phenosample$Diagnosis=="dup15q", ] -> phenosample
read.csv("phenosamplecopy.csv", header=T) -> phenosample
# see structure of CB 
str(CB)
# see list of sentix IDs to filter down data 
sampleNames(CB)
CB[ ,sampleNames(CB) %in% phenosample$Sentix]-> CB.sub 
# extract CB betas into 2 dimensional object 
betas(CB) -> betas.CB 
# slice CB data, remove the 8 15q samples
betas.CB[ ,colnames(betas.CB) %in% phenosample$Sentix] -> sliced.betas.CB
# check beta file matches pheno file, reorder to match 
phenosample[order(phenosample$Sentix),]-> phenosample 
sliced.betas.CB[, order(colnames(sliced.betas.CB))]->sliced.betas.CB
# save sliced betas output 
save(sliced.betas.CB, file='sliced.betas.CB.rda')
# check that first 5 sentix match between betas and pheno 
phenosample$Sentix[1:5]
colnames(sliced.betas.CB)[3,6]
sliced.betas.CB[1:3,1:6]
# run initial QC on beta file 
library('ChAMP')
champ.QC(beta=sliced.betas.CB,
         pheno=phenosample$Diagnosis,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=FALSE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages/")
# initial normalization  
na.omit(sliced.betas.CB)->sliced.betas.CB2 
out <- champ.norm(beta=sliced.betas.CB2,
           resultsDir="./CHAMP_Normalization/",
           method="BMIQ",
           plotBMIQ=TRUE,
           arraytype="450K",
           cores=3)
# remove outlier from data 
load(sliced.betas.CB)
sliced.betas.CB[ ,c(1:44, 46:60)]-> sliced.betas.CB
colnames(sliced.betas.CB)
# removed "3999979145_R02C02", remove same sample from pheno file 
phenosample[c(1:44,46:60),]-> phenosample 

# rerun QC and normalisation without outlier 
champ.QC(beta=sliced.betas.CB,
         pheno=phenosample$Diagnosis,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=FALSE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_no_outlier/")
# normalisation 
na.omit(sliced.betas.CB)->sliced.betas.CB2 
out <- champ.norm(beta=sliced.betas.CB2,
                  resultsDir="./CHAMP_Normalization_no_outlier/",
                  method="BMIQ",
                  plotBMIQ=TRUE,
                  arraytype="450K",
                  cores=3)

# rename random phenofiles to be background
X3999979055_R02C02
X3999979063_R01C02 

# check that diagnosis have been changed to background 
read.csv("phenosamplecopy_background.csv", header=T) -> phenosamplecopy_background
# remove outlier from phenosamplecopy_background
phenosamplecopy_background[c(1:46,48:60),]-> phenosamplecopy_background
# make intensity files 
CB[ ,sampleNames(CB) %in% phenosamplecopy_background$Sentix]-> CB.sub 
# generate CNVs 
cnv_results <- champ.CNA(
  intensity = total.intensity(CB.sub),
  pheno = phenosamplecopy_background$Diagnosis, # analyse by column
  control = TRUE,  # using own background controls
  controlGroup = "BACKGROUND", # background samples in the pheno column
  sampleCNA = TRUE,  # CNA by individual samples
  groupFreqPlots = TRUE,  # frequency by groups
  Rplot = FALSE,  # Do you want R to print plots as it runs?
  PDFplot = TRUE,  # Do you want R to save PDF plots for you?
  freqThreshold = 0.3,  # This is the detection level for a CNV. 
  resultsDir = "/mnt/data1/hazel/CHAMP.CNA", # save results
  arraytype = "450K"  # Which array are we using (450k)
)
save (cnv_results,file='cnv_results.rda')

# save into tables for statistics 
load("cnv_results.rda")
> write.csv(cnv_results$groupResult$`BACKGROUND `, file='cnv_results_BACKGROUND.csv')
> write.csv(cnv_results$groupResult$CTL, file='cnv_results_CTL.csv')
> write.csv(cnv_results$groupResult$ASD, file='cnv_results_ASD.csv')
  
