# Finding_iML

![Python](https://img.shields.io/badge/Python-3.9%20|%203.10%20|%203.12-blue)
![scRNA-seq](https://img.shields.io/badge/analysis-scRNA--seq-purple)
![NK cells](https://img.shields.io/badge/cell%20type-NK%20cells-teal)
![pan-cancer](https://img.shields.io/badge/scope-pan--cancer-orange)
![Machine Learning](https://img.shields.io/badge/method-machine%20learning-yellow)

Code repository for the paper:
> **Memory-like NK cells naturally occur in healthy individuals and are enriched within solid tumors**

Classifier scripts and cohort analysis pipelines for identifying and analyzing memory-like (ML) NK cells using machine learning.

---

## 📖 Documentation

Full documentation is included as `Finding_eML.zip`. Decompress the archive and open the local HTML docs:

```bash
# macOS
open Finding_eML/docs/_build/html/index.html

# Linux
xdg-open Finding_eML/docs/_build/html/index.html

# Windows
start Finding_eML/docs/_build/html/index.html
```

---

## ⚙️ Requirements

**🐳 Docker images**

| Used for | Image |
|----------|-------|
| Classifier scripts | `evelyns2000/foltz_tools:classifier2024` |
| Finding eML / iML package | `veda504/finding_eml:v1.1` |
| R scripts (scRNA / multiome) | `kaushalmadhurima/scrna_multiome` |

**🐍 Python — Jupyter notebook scripts**

| Python version | Used for |
|----------------|----------|
| 3.9.6 | Most scripts (default) |
| 3.10.14 | `marsilea`, `pydeseq2` |
| 3.12.12 | `hucira` |

**📊 R scripts** — version `4.5.2` · Docker: `kaushalmadhurima/scrna_multiome`

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
Finding_eML.zip       ← full documentation archive
README.md
```

---

## 🔬 Classifier scripts

The `classifier_scripts/` directory contains the core ML classifiers, each tailored to a specific cohort's batch correction strategy.

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

## 📂 Cohort analysis scripts

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

**Note on naming:** Within scripts, dataset names and author names are used interchangeably — Netskar et al. data may appear as `malmberg` or `malm`; Xing et al. data may appear as `fan`. Additionally, iML are referred to as `eML` in some scripts.
