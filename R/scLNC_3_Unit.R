
#' Calculate lncRNA units activity score.
#'
#' @param objectInput A scLNC object including expression matrix slot.
#' @param AUC The Boolean value determines whether to calculate the activity of the unit. The default value is TRUE.
#' @param displayLncRNA lncRNA list. Only keep a subset of lncRNA units, defaults to all lncRNA units.
#' @param DEAUC Boolean values determine whether to compare differential activity of units between two groups. The default value is TRUE.
#' @param item.add A new list of cell information, including new cell item column.
#' @param DEitem ttribute for comparison.
#'
#'
#' @return Return a scLNC object with AUC slot.
#' @export
#'
#' @examples
#' \dontrun{
#'data(LNCobject)
#'LNCobject=scLNC_3_Unit (objectInput=LNCobject,AUC=TRUE,displayLncRNA=NULL,DEAUC=TRUE,item.add=NULL,DEitem='majorCluster')
#'}
#'
scLNC_3_Unit <- function(objectInput,AUC=TRUE,displayLncRNA=NULL,DEAUC=TRUE,item.add=NULL,DEitem='majorCluster'){
  if(AUC){
    objectInput <- AUCell_score(object=objectInput,lnclist = displayLncRNA)
  }


  if(DEAUC){
    if(!(is.null(item.add))){
      LRAT.object_T@ cell.info=item.add

    }

    objectInput <- DeActivity(object=objectInput,item=DEitem,FC=0.1,pvalue=0.05,min.pct=0.3,padj=1)
  }

  return(objectInput)
}

