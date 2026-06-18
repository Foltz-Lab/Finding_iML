#########################################################################################
####  In this script, we are merging the bigwigs so we have a representative iML, eML, cNK, and effcNK bigwig for
####  visualization in deepTools. First, we merge the bigwigs to get a bedgraph file, then we sort the bedgraph before
####  using bedGraphToBigWig to get the final merged files.
#########################################################################################

# Merge the bigwigs for each sample group and get merged bedgraph
bsub -oo mergebigwigs_iML.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/ucsc-bigwigmerge:377--h446ed27_1)" \
bigWigMerge /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/147429_IL15_PanML_reNorm_1.9.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/153147_IL15_PanML_reNorm_0.34.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/146724_IL15_PanML_reNorm_0.97.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/172784_IL15_PanML_reNorm_1.36.bw iML.bedgraph

bsub -oo mergebigwigs_eML.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/ucsc-bigwigmerge:377--h446ed27_1)" \
bigWigMerge /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/147429_IL121518_PanML_reNorm_2.03.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/153147_IL121518_PanML_reNorm_0.81.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/146724_IL121518_PanML_reNorm_1.75.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/172784_IL121518_PanML_reNorm_0.74.bw eML.bedgraph

bsub -oo mergebigwigs_cNK.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/ucsc-bigwigmerge:377--h446ed27_1)" \
bigWigMerge /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/147429_IL15_cNK_reNorm_0.85.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/153147_IL15_cNK_reNorm_1.24.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/146724_IL15_cNK_reNorm_1.28.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/172784_IL15_cNK_reNorm_1.04.bw cNK.bedgraph

bsub -oo mergebigwigs_effcNK.out -G compute-jennifer.a.foltz -q general -R 'rusage[mem=100G]' -a "docker(quay.io/biocontainers/ucsc-bigwigmerge:377--h446ed27_1)" \
bigWigMerge /storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/147429_IL121518_ML3_reNorm_0.99.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/153147_IL121518_ML3_reNorm_0.68.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/146724_IL121518_ML3_reNorm_1.1.bw \
/storage1/fs1/jennifer.a.foltz/Active/matthewmueller/JT048_JT051/bigwigs_MLsubsets/172784_IL121518_ML3_reNorm_0.66.bw effcNK.bedgraph


# sort bedgraphs

LC_COLLATE=C sort -k1,1 -k2,2n iML.bedgraph > sorted.iML.bedgraph


# Convert bedgraph back to bigwig
bsub -oo bigwig_iML.out -G compute-jennifer.a.foltz -q general  -R 'rusage[mem=100G]'  -a "docker(quay.io/biocontainers/ucsc-bedgraphtobigwig:357--h446ed27_4)" bedGraphToBigWig sorted.iML.bedgraph hg38.chrom.sizes iML.bw
bsub -oo bigwig_iML.out -G compute-jennifer.a.foltz -q general  -R 'rusage[mem=100G]'  -a "docker(quay.io/biocontainers/ucsc-bedgraphtobigwig:357--h446ed27_4)" bedGraphToBigWig sorted.eML.bedgraph hg38.chrom.sizes eML.bw
bsub -oo bigwig_iML.out -G compute-jennifer.a.foltz -q general  -R 'rusage[mem=100G]'  -a "docker(quay.io/biocontainers/ucsc-bedgraphtobigwig:357--h446ed27_4)" bedGraphToBigWig sorted.cNK.bedgraph hg38.chrom.sizes cNK.bw
bsub -oo bigwig_iML.out -G compute-jennifer.a.foltz -q general  -R 'rusage[mem=100G]'  -a "docker(quay.io/biocontainers/ucsc-bedgraphtobigwig:357--h446ed27_4)" bedGraphToBigWig sorted.effcNK.bedgraph hg38.chrom.sizes effcNK.bw
