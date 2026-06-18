#########################################################################################
####  In this script, we will check the results of the bigwig renormalization to make sure normalized values fall 
####  within the expected values. We can visualize the normalization quality using computeMatrix and plotHeatmap from deepTools
#########################################################################################


bsub -oo normCheck_mtx.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
computeMatrix reference-point -o normCheck_mtx.gz -p 12 \
-S /storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/146724_IL121518_PanML_reNorm_1.75.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/146724_IL15_cNK_reNorm_1.28.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/146724_IL15_PanML_reNorm_0.97.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/153147_IL121518_PanML_reNorm_0.81.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/153147_IL15_cNK_reNorm_1.24.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/153147_IL15_PanML_reNorm_0.34.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/147429_IL121518_PanML_reNorm_2.03.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/147429_IL15_cNK_reNorm_0.85.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/147429_IL15_PanML_reNorm_1.9.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/172784_IL121518_PanML_reNorm_0.74.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/172784_IL15_cNK_reNorm_1.04.bw \
/storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bigwigs_MLsubsets/172784_IL15_PanML_reNorm_1.36.bw \
-R /storage1/fs1/tfehnige/Active/JTran/JT048_JT051/bed_files/hg38.HouseKeepingGenes.bed -a 2000 -b 2000 --missingDataAsZero


bsub -oo normCheck_heatmap.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/deeptools:3.5.1--py_0)" \
plotHeatmap -m normCheck_mtx.gz -out normCheck_heatmap.pdf --colorMap Reds \
--samplesLabel 146724_eML_reNorm 146724_cNK_reNorm 146724_iML_reNorm \
153147_eML_reNorm 153147_cNK_reNorm 153147_iML_reNorm 147429_eML_reNorm \
147429_cNK_reNorm 147429_iML_reNorm 172784_eML_reNorm 172784_cNK_reNorm 172784_iML_reNorm