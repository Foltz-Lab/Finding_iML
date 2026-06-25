# Finding_iML
![Python](https://img.shields.io/badge/Python-3.9%20|%203.10%20|%203.12-blue)
![scRNA-seq](https://img.shields.io/badge/analysis-scRNA--seq-purple)
![NK cells](https://img.shields.io/badge/cell%20type-NK%20cells-teal)
![pan-cancer](https://img.shields.io/badge/scope-pan--cancer-orange)
![Machine Learning](https://img.shields.io/badge/method-machine%20learning-yellow)

Code repository for the paper:
> **Memory-like NK cells naturally occur in healthy individuals and are enriched within solid tumors**

Package documentation, classifier scripts and cohort analysis pipelines for identifying and analyzing memory-like (ML) NK cells using machine learning.

---

## 🗂 Repository structure

```
classifier_scripts/
    01_classifier_buildandclassifyquery_batch_dataset_v2_Netskar_et_al.py
    02_classifier_buildandclassifyquery_batch_datasource_chemistry_Ding_et_al.py
    03_classifier_buildandclassifyquery_batch_none_Xing_et_al.py
scripts/
    bulk_ATAC_seq/
    Cohort1_Netskar_et_al/
    Cohort2_Ding_et_al/
    Cohort3_SI_CIML/
    Cohort4_Xing_et_al/
    PBMC_combined/
    Pan_cancer_combined/
Finding_eML.zip       ← full package documentation archive
README.md
```

---

## 📖 Package Documentation - iML/eML
Full documentation is hosted online:

📄 **[https://finding-eml.readthedocs.io](https://finding-eml.readthedocs.io)**

---

## ⚙️ Requirements for scripts

**🐳 Docker images**

| Used for | Image |
|----------|-------|
| Finding eML / iML package | `veda504/finding_eml:v1.1` |
| Classifier scripts | `evelyns2000/foltz_tools:classifier2024` |
| R scripts (scRNA / multiome) | `kaushalmadhurima/scrna_multiome` |

**🐍 Python — Jupyter notebook scripts**

| Python version | Used for |
|----------------|----------|
| 3.9.6 | Most scripts (default) |
| 3.10.14 | `Marsilea`, `PyDESeq2` |
| 3.12.12 | `huCIRA` |

**📊 R scripts** — local version `4.5.2` / Docker: `kaushalmadhurima/scrna_multiome`

**🧬 Bulk ATAC-seq scripts** 
 
| Used for | Image |
|----------|-------|
| Seurat v4.1 with DSB | `jafoltz/seuratv4.1withdsb` |
| BigWig processing | `jafoltz/biggywiggy:1.0` |
| deepTools | `quay.io/biocontainers/deeptools:3.5.1--py_0` |
| DiffBind | `naotokuboda/diffbind:3.2` |
| UCSC BigWigMerge | `quay.io/biocontainers/ucsc-bigwigmerge:377--h446ed27_1` |

---

## 🔬 Classifier scripts

The `classifier_scripts/` directory contains the core ML classifier scripts, each tailored to a specific cohort's batch correction strategy.

**01 — Netskar et al.**  
`01_classifier_buildandclassifyquery_batch_dataset_v2_Netskar_et_al.py`  
Builds and applies the ML NK cell classifier with dataset-level batch correction.

**02 — Ding et al.**  
`02_classifier_buildandclassifyquery_batch_datasource_chemistry_Ding_et_al.py`  
Classifier with batch correction accounting for data source and sequencing chemistry.

**03 — Xing et al.**  
`03_classifier_buildandclassifyquery_batch_none_Xing_et_al.py`  
Classifier applied without batch correction for the Xing et al. cohort.

---

## 🔍 Cohort analysis scripts

The `scripts/` directory contains per-cohort downstream analysis pipelines:

| Folder | Description |
|--------|-------------|
| `bulk_ATAC_seq/` | Chromatin accessibility analysis |
| `Cohort1_Netskar_et_al/` | Cohort 1 analysis — Netskar et al. |
| `Cohort2_Ding_et_al/` | Cohort 2 analysis — Ding et al. |
| `Cohort3_SI_CIML/` | Cohort 3 analysis — SI CIML |
| `Cohort4_Xing_et_al/` | Cohort 4 analysis — Xing et al. |
| `PBMC_combined/` | Combined PBMC cohort analysis |
| `Pan_cancer_combined/` | Pan-cancer integrated cohort analysis |

---

💡 **Note on naming:** Within scripts, dataset names and author names are used interchangeably — Netskar et al. data may appear as `malmberg` or `malm`; Xing et al. data may appear as `fan`. Additionally, iML are referred to as `eML` in some scripts.

## 📚 Citations

If you use `Finding eML` package in your research, please cite:

**Finding eML:**
Foltz JA, Tran J, Wong P, et al. Cytokines drive the formation of memory-like NK cell subsets via epigenetic rewiring and transcriptional regulation. *Science Immunology.* 2024;9(96):eadk4893. https://doi.org/10.1126/sciimmunol.adk4893

**TotalVI:**
Gayoso A, Steier Z, Lopez R, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. *Nat Methods.* 2021;18(3):272-282. https://doi.org/10.1038/s41592-020-01050-x

**scVI (scvi-tools):**
Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol.* 2022;40(2):163-166. https://doi.org/10.1038/s41587-021-01206-w

**SCANPY:**
Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 2018;19(1):15. https://doi.org/10.1186/s13059-017-1382-0

**Imbalanced-learn:**
Lemaitre G, Nogueira F, Aridas CK. Imbalanced-learn: A Python Toolbox to Tackle the Curse of Imbalanced Datasets in Machine Learning. *arXiv.* 2016. https://doi.org/10.48550/ARXIV.1609.06570

