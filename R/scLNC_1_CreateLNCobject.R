
#' Create a scLNC object for analysis.
#'
#' @param count_table1 Gene count matrix with cells as columns and features as rows.
#' @param cell_info1 A list of cell information. A data.frame where the rows are cell names and the columns are cell informations. Row names in the cell_info1 need to match the column names of the count_table1.
#' @param geneid Boolean values determining if gene names convert to Ensemble IDs. Default is TRUE.
#' @param create Boolean values determine whether the object is to be created. The default value is TRUE.
#' @param splitobject Boolean values determine whether objects should be split into items. Default is TRUE.
#' @param split.item1 Attribute for splitting.
#'
#'
#'
#' @return Return a scLNC object or a list of scLNC object by different items.
#' @export
#'
#' @examples
#'  \dontrun{
#' data(HCC_rawcount)
#' data(HCC_cellinfo)
#'
#' object=scLNC_1_CreateLNCobject (count_table1= HCC_rawcount,
#' cell_info1= HCC_cellinfo, geneid=TRUE, splitobject=FALSE, split.item1=NA)
#' object.list=scLNC_1_CreateLNCobject (count_table1= HCC_rawcount,
#' cell_info1= HCC_cellinfo, geneid=TRUE, splitobject=TRUE,split.item1="Tissue")
#'}
#'
scLNC_1_CreateLNCobject <- function(count_table1, cell_info1, geneid=TRUE,create=TRUE, splitobject=TRUE, split.item1){
if(geneid){
rownames(count_table1)=id_convert_fromgtf(genelist=rownames(count_table1), gtf.info = scLNCgencode)#name转id
}
if(create){
LNCobject <- createLRAT(count_table=count_table1, cell_info=cell_info1, min.cells = 15, min.genes = 200)#创建object
}
if(splitobject){
LNCobject.list <- splitLRAT(LNCobject, split.item=split.item1)#拆分object
return(LNCobject.list)
}else{
return(LNCobject)
}
}


