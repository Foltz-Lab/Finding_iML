#########################################################################################
####  In this script, we generate bigwig files from our final bam files. After generating the bigwigs, 
####  we run DiffBind wth DESeq2 to get the size factors for re-normalization. Then we can re-normalize the
####  bigwig files with the 'docker(jafoltz/biggywiggy:1.0)' and the Rscript bwReNorm.R
####  The re-normalized bigwig files can be used for data visualization with tools such as IGV and deepTools
#########################################################################################

##Generate bam to bigwig files
cd /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051

######Donor 146724#######
bsub -oo bamtobw1.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.146724_IL121518_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/146724_IL121518_PanML.bw -p 12

bsub -oo bamtobw1.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.146724_IL121518_ML3.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/146724_IL121518_ML3.bw -p 12

bsub -oo bamtobw3.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.146724_IL15_cNK.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/146724_IL15_cNK.bw -p 12

bsub -oo bamtobw4.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.146724_IL15_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/146724_IL15_PanML.bw -p 12

######Donor 147429#######
bsub -oo bamtobw5.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.147429_IL121518_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/147429_IL121518_PanML.bw -p 12

bsub -oo bamtobw6.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.147429_IL121518_ML3.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/147429_IL121518_ML3.bw -p 12

bsub -oo bamtobw7.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.147429_IL15_cNK.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/147429_IL15_cNK.bw -p 12

bsub -oo bamtobw8.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.147429_IL15_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/147429_IL15_PanML.bw -p 12

######Donor 153147#######
bsub -oo bamtobw9.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.153147_IL121518_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/153147_IL121518_PanML.bw -p 12

bsub -oo bamtobw10.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.153147_IL121518_ML3.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/153147_IL121518_ML3.bw -p 12

bsub -oo bamtobw11.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.153147_IL15_cNK.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/153147_IL15_cNK.bw -p 12

bsub -oo bamtobw12.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.153147_IL15_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/153147_IL15_PanML.bw -p 12

######Donor 172784#######
bsub -oo bamtobw13.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.172784_IL121518_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/172784_IL121518_PanML.bw -p 12

bsub -oo bamtobw14.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.172784_IL121518_ML3.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/172784_IL121518_ML3.bw -p 12

bsub -oo bamtobw15.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.172784_IL15_cNK.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/172784_IL15_cNK.bw -p 12

bsub -oo bamtobw16.out -G compute-oncology -q oncology  -R 'rusage[mem=100G]' \
-a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
bamCoverage -b /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/finalbam/final.172784_IL15_PanML.bam \
--normalizeUsing None -o ./bigwigs_MLsubsets/172784_IL15_PanML.bw -p 12

# Get these size factors from the DiffBind object after running MLsubsets_Diffbind.R at object$DESeq2$DEdata$sizeFactor
#####SIZE FACTORS#############
#Copy factors below for bigwig reNorm
#X146724 X146724.1 X146724.2 X146724.3   X147429 X147429.1 X147429.2 X147429.3 
#1.7459019 1.1043424 1.2762536 0.9749536 2.0340275 0.9896347 0.8452760 1.9037192 
#X153147 X153147.1 X153147.2 X153147.3   X172784 X172784.1 X172784.2 X172784.3 
#0.8122087 0.6809684 1.2409310 0.3437656 0.7393783 0.6566474 1.0432068 1.3591269 

# copy below command for running bwReNorm.R
LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' bsub -oo getthemdone%J.out -G compute-jennifer.a.foltz -q general -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' "Rscript bwReNorm.R -i pathtosamplesbigwigfilehere.bw -f sizefactornumberhere"

#be sure to input -i for each sample you run
#
#and -f for each sample you run
#
######Renorm Donor 146724#######
LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/146724_IL121518_PanML.bw -f 1.7459019"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/146724_IL121518_ML3.bw -f 1.1043424"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/146724_IL15_cNK.bw -f 1.2762536"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/146724_IL15_PanML.bw -f 0.9749536"

######Renorm Donor 147429#######
LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/147429_IL121518_PanML.bw -f 2.0340275"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/147429_IL121518_ML3.bw -f 0.9896347"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/147429_IL15_cNK.bw -f 0.8452760"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/147429_IL15_PanML.bw -f 1.9037192"

######Renorm Donor 153147#######
LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/153147_IL121518_PanML.bw -f 0.8122087"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/153147_IL121518_ML3.bw -f 0.6809684"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/153147_IL15_cNK.bw -f 1.2409310"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/153147_IL15_PanML.bw -f 0.3437656"

######Renorm Donor 172784#######
LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/172784_IL121518_PanML.bw -f 0.7393783"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/172784_IL121518_ML3.bw -f 0.6566474"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs/172784_IL15_cNK.bw -f 1.0432068"

LSF_DOCKER_VOLUMES='/storage1/fs1/jennifer.a.foltz/Active:/storage1/fs1/jennifer.a.foltz/Active /scratch1/fs1/jennifer.a.foltz:/scratch1/fs1/jennifer.a.foltz' \
bsub -oo getthemdone%J.out -G compute-oncology -q oncology -M 100G -R 'select[mem>100G] rusage[mem=100G]' -a 'docker(jafoltz/biggywiggy:1.0)' \
"Rscript bwReNorm.R -i /storage1/fs1/tfehnige/Active/matthewmueller/JT048_JT051/bigwigs/172784_IL15_PanML.bw -f 1.3591269"



