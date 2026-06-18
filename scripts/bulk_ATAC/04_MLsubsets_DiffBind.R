#########################################################################################
####  In this script, we use the DiffBind package (v3.2) to get differential peaks with DESeq2. The output of this
####  script can be input into Chipseeker to annotate the DARs, and it can be used to visualize the DARs in volcano plots.
####  A PCA plot of all sites labeled by batch will also be generated in this script.
#########################################################################################


#Before running this Script be sure to input your sample arguements and create Saample sheet
#Create sampleSheet with the following header columns (save as csv file)
#SampleID, Condition, Treatment, bamReads, Peaks, Peakcaller, Replicate
# below is the docker to use. Copy and paste everything after the #
#LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' bsub -Is -G compute-oncology -q oncology-interactive -M 200G -R 'select[mem>200G] rusage[mem=200G]' -a 'docker(naotokubota/diffbind:3.2)' /bin/bash
#Before proceeding, start R, by type R then press enter
library(DiffBind)
print("these are my variables for this run")
sampleSheet = "/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/MLsubsets_samplesheet.csv" 
# name your variables
sample = "MLsubsets"
resultsdir = "MLsubsets_DiffBind"
setwd(resultsdir)
date = gsub("-","",Sys.Date());
run = "MLsubsets"
print("starting diffbind analysis")
options(echo = TRUE)
##Create dba object
object <- dba(sampleSheet = sampleSheet)
object <- dba.blacklist(object, greylist = FALSE)

##Generate affinity binding matrix
object <- dba.count(object)

dba.save(object, file = sprintf("%s_%s_%s_dbacount", sample, date, run))

###20230816#####
##Load in object from dba analysis###
# if you want to save and load in the future this is an example: object <- dba.load(file = "MLsubsets_20230105_MLsubsets_dba_analysis")
object <- dba.normalize(object, method=DBA_DESEQ2, normalize=DBA_NORM_NATIVE, library = DBA_LIBSIZE_DEFAULT) # choose these parameters to be default DESeq2
object <- dba.contrast(object, design = "~Replicate + Treatment", minMembers=2)
print("visualizing contrast groups")
object <- dba.analyze(object, method=DBA_DESEQ2)
object
dba.show(object, bContrasts = TRUE) # visualize number of DEGs in each contrast
analysis.DB <- dba.report(object) # this is your output, which you will need to convert to gene names later, defaults to the first comparison
analysis.DB

##To get the sample size factor to generate appropriate bigwigs

#Copy factors below for bigwig reNorm
object$DESeq2$DEdata$sizeFactor
X146724 X146724.1 X146724.2 X146724.3   X147429 X147429.1 X147429.2 X147429.3 
1.7459019 1.1043424 1.2762536 0.9749536 2.0340275 0.9896347 0.8452760 1.9037192 
X153147 X153147.1 X153147.2 X153147.3   X172784 X172784.1 X172784.2 X172784.3 
0.8122087 0.6809684 1.2409310 0.3437656 0.7393783 0.6566474 1.0432068 1.3591269 

dba.save(object, sprintf("%s_%s_%s_dba_analysis", sample, date, run))


#PCA plot of all comparisons
pdf(sprintf("%s_%s_%spca_ofallsites_labeled_batch.pdf", sample, date, run), width = 10, height = 10)
dba.plotPCA(object, attributes = DBA_TREATMENT, score = DBA_SCORE_RPKM, label = DBA_FACTOR, vColors=c("gray45", "purple2", "mediumblue", "green4")) # pca of all sites
dev.off()


##CONTRAST 7
contrast = 7 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))

##CONTRAST 8
contrast = 8 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))


##CONTRAST 9
contrast = 9 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))


##CONTRAST 10
contrast = 10 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))

##CONTRAST 11
contrast = 11 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))

##CONTRAST 12
contrast = 12 # insert number here 
analysis.DB <- dba.report(object, contrast=contrast) # here you specify which data you want to see, here we wanted #3 from contrast so contrast=3
save(analysis.DB, file = sprintf("%s_%s_%s_contrast%s_granges.Rdata", sample, date, run, contrast)) # this is the input into chippeak and your differential peakresults

#Obtain file for volcano data
plot <- dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod, th=1, fold=0)
write.csv(plot, file = sprintf("%s_%s_%s_contrast%svolcanodata.csv", sample, date, run, contrast))

#Plot volcano from DiffBind as a reference for x and y axis 

contrast = 7 
pdf("20230815_contrast7_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

contrast = 8 
pdf("20230815_contrast8_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

contrast = 9
pdf("20230815_contrast9_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

contrast = 10
pdf("20230815_contrast10_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

contrast = 11
pdf("20230815_contrast11_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

contrast = 12
pdf("20230815_contrast12_volcano.pdf", width = 10, height = 10)
dba.plotVolcano(object, contrast =contrast, method = object$config$AnalysisMethod)
dev.off()

#PCA plot of all comparisons
pdf(sprintf("%s_%s_%spca_ofallsites_labeled_batch.pdf", sample, date, run), width = 8, height = 8)
dba.plotPCA(object, attributes = DBA_TREATMENT, dotSize = 2, score = DBA_SCORE_RPKM, label = DBA_FACTOR, vColors=c("cyan4", "slateblue", "grey20", "turquoise")) # pca of all sites
dev.off()





# if interactive: 
# q()
# Save workspace image? [y/n/c]: n

