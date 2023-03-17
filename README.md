
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

## Example data

An example data can be downloaded from
[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638) with
the commands:

``` r
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.count.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/miniml/GSE98638_family.xml.tgz
```

## Examples

1.  Create a scLNC object first.

``` r
LNCobject=scLNC_1_CreateLNCobject (count_table1= HCC_rawcount, cell_info1= HCC_cellinfo)
```

2.  Build lncRNA units and annotate the lncRNA-mRNA co-expression pairs.

``` r
LNCobject=scLNC_2_Relation(objectInput=LNCobject,relation.fileName='HCC',pairs.filter=TRUE,corcut=0.6,anno.pairs=TRUE)
```

3.  Calculate the activity score of lncRNA units.

``` r
LNCobject=scLNC_3_Unit (objectInput=LNCobject,AUC=TRUE,displayLncRNA=NULL,DEAUC=TRUE,item.add=NULL,DEitem='majorCluster')
```

4.  Compare units from different experimental conditions.

``` r
scLNC_4_Compare(DisplayUnit=TRUE,ob.ls,lncRNA,DeGO=TRUE,CopareGOfile,CopareGeneList)
```
