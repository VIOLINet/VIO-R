# VIO-R

## limma_one_series.R

Run limma for RNA microarray data and compare two or more groups of GEO Samples in a GEO Series.

### Dependencies
R version >=3.6.1

Biobase version >=2.46.0

GEOquery version >=2.54.1

GEOmetadb version >=1.48.0

limma version >=3.42.0

### Example

Rscript limma_one_series.R GSM733843,GSM733844 GSM733852,GSM733853

Rscript limma_one_gse.R GSM733843,GSM733852 GSM733844,GSM733853 GSM733845,GSM733854
