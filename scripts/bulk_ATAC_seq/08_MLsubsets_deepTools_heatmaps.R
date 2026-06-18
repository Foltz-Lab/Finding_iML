#########################################################################################
####  In this script, we will use deepTools' computeMatrix to generate the data needed for plotHeatmap. The output from this
####  script will be a matrix holding significant DARs for visualization, and a heatmap for each NK type that highlights
####  the DARs in that NK type compared to the other two (eg. DARs up in iML compared to eML and cNK)
#########################################################################################


# Run computematrix on each sample

bsub -oo eML_mtx.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
computeMatrix reference-point -R /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/eML_contrast8_diffpeaks_significant.bed \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/eML_contrast9_diffpeaks_significant.bed \
-S /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/eML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/iML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/cNK.bw \
-o eML_mtx.gz -b 2000 -a 2000 -bs 50 -p 12

bsub -oo iML_mtx.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
computeMatrix reference-point -R /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/iML_contrast9_diffpeaks_significant.bed \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/iML_contrast12_diffpeaks_significant.bed \
-S /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/eML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/iML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/cNK.bw \
-o iML_mtx.gz -b 2000 -a 2000 -bs 50 -p 12

bsub -oo cNK_mtx.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
computeMatrix reference-point -R /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/cNK_contrast8_diffpeaks_significant.bed \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bed_files/cNK_contrast12_diffpeaks_significant.bed \
-S /storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/eML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/iML.bw \
/storage1/fs1/jennifer.a.foltz/Active/Matt/JT048_JT051/bigwigs_MLsubsets/cNK.bw \
-o cNK_mtx.gz -b 2000 -a 2000 -bs 50 -p 12



# Use matrix to plot data

bsub -oo eML_heatmap.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
plotHeatmap -m eML_mtx.gz --colorMap Reds --missingDataColor 0.85 --regionsLabel eML-vs-cNK eML-vs-iML --samplesLabel eML.bw iML.bw cNK.bw -o eML_heatmap.pdf

bsub -oo iML_heatmap.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
plotHeatmap -m iML_mtx.gz --colorMap Reds --missingDataColor 0.85 --regionsLabel iML-vs-eML iML-vs-cNK --samplesLabel eML.bw iML.bw cNK.bw -o iML_heatmap.pdf

bsub -oo cNK_heatmap.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
plotHeatmap -m cNK_mtx.gz --colorMap Reds --missingDataColor 0.85 --regionsLabel cNK-vs-eML cNK-vs-iML --samplesLabel eML.bw iML.bw cNK.bw -o cNK_heatmap.pdf


