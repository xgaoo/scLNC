
#' Perform calculation, filtration and annotation of lncRNA-mRNA co-expression relationship pairs.
#'
#' @param objectInput A scLNC object including expression matrix slot.
#' @param relation Boolean values determine whether to calculate the co-expression of lncRNA-mRNA. The default value is TRUE.
#' @param relation.fileName The name of the output file.
#' @param pairs.filter The Boolean value determines whether the calculated relationship pairs should be filtered. The default value is TRUE.
#' @param corcut Choose the relationship pairs with higher correlation coefficient.
#' @param anno.pairs Boolean values determine whether the relationship pairs should be annotated, including TF, cytokine, known lncRNA-mRNA databases, promoter and enhancer information. The default value is TRUE.
#'
#'
#'
#' @return A scLNC object including lncRNA-mRNA co-expression relationship pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject)
#'LNCobject=scLNC_2_Relation(objectInput=LNCobject,relation.fileName='test',
#'pairs.filter=TRUE,corcut=0.6,anno.pairs=TRUE)
#'}
#'
scLNC_2_Relation <- function(objectInput,relation=TRUE,relation.fileName,pairs.filter=TRUE,corcut,anno.pairs=TRUE){
if(relation){
rdata=CaculateCorelation(object=objectInput,fileName=relation.fileName)
}
if(pairs.filter){
objectInput=FilterPairs(data=rdata, object=objectInput, corcut=corcut)
}
if(anno.pairs){
objectInput=PairsAnnotation(object=objectInput,TF=TRUE,cytokine=TRUE,LMpairs=TRUE)
}
return(objectInput)
 }








