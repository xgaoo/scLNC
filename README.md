
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scLNC

<!-- badges: start -->
<!-- badges: end -->

The goal of scLNC is to investigate lncRNA function in single-cell
RNA-seq.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("xgaoo/scLNC")
```

\##Examples

1.Create a scLNC object for analysis.

``` r
LNCobject=scLNC_1_CreateLNCobject (count_table1= HCC_rawcount, cell_info1= HCC_cellinfo)
```

2.Perform calculation, filtration and annotation of lncRNA-mRNA
co-expression relationship pairs.

``` r
LNCobject=scLNC_2_Relation(objectInput=LNCobject,relation.fileName='test',pairs.filter=TRUE,corcut=0.6,anno.pairs=TRUE)
```

3.Calculate lncRNA units activity score.

``` r
LNCobject=scLNC_3_Unit (objectInput=LNCobject,AUC=TRUE,displayLncRNA=NULL,DEAUC=TRUE,item.add=NULL,DEitem='majorCluster')
```

4.Compare units from different items.

``` r
scLNC_4_Compare(DisplayUnit=TRUE,ob.ls,lncRNA,DeGO=TRUE,CopareGOfile,CopareGeneList)
```
