# scds

**s**ingle **c**ell **d**oublet **s**coring: In-silico doublet annotation for single cell RNA sequencing data

```scds``` is an ```R``` package for computational doublet annotation of single cell RNA sequencing data. It interfaces with the S4 ```SingleCellExperiment``` class [(see here)](https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment), so it should easily integrate into many ```R```/Bioconductor scRNA-seq analysis workflows. 

#### Installation

```
devtools::install_github('kostkalab/scds',ref="nBC")
```

#### Quick Start

In the following ```sce``` is a ```SingleCellExperiment``` holding at least raw counts in an assay called ```counts```.

```
#- Annotate doublet using co-expression based doublet scoring:
sce = cxds(sce)

#- Annotate doublet using binary classification based doublet scoring:
sce = bcds(sce)

#- Combine both annotations into a hybrid annotation
sce = cxds_bcds_hybrid(sce)

#- Doublet scores are now available via colData:
CD  = colData(sce)
head(cbind(CD$cxds_score,CD$bcds_score, CD$hybrid_score))
rm(CD)

```

#### Other doublet detection tools:
* [DoubletCells](https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-6-doublet.html) as part of [SimpleSingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html)
* [DoubletDecon](https://github.com/EDePasquale/DoubletDecon)  
* [DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
* [Scrublet](https://github.com/AllonKleinLab/scrublet)



