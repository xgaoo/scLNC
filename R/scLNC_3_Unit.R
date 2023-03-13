
#' Calculate lncRNA units activity score.
#'
#' @param objectInput A scLNC object including expression matrix slot.
#' @param AUC The Boolean value determines whether to calculate the activity of the unit. The default value is TRUE.
#' @param displayLncRNA lncRNA list. Only keep a subset of lncRNA units, defaults to all lncRNA units.
#' @param Go The Boolean value determines whether to visualize the functional enrichment results of several units. The default value is FALSE.
#' @param MetascapeGofile Metascape results.
#'
#'
#' @return Return a scLNC object with AUC slot.
#' @export
#'
#' @examples
#' \dontrun{
#'data(LNCobject)
#'LNCobject=scLNC_3_Unit (objectInput=LNCobject,AUC=TRUE,displayLncRNA=NULL)
#'}
#'
scLNC_3_Unit <- function(objectInput,AUC=TRUE,displayLncRNA=NULL,Go=FALSE,
                         MetascapeGofile){
if(AUC){
objectInput <- AUCell_score(object=objectInput,lnclist = displayLncRNA)
return(objectInput)
 }
if(Go){
CompareGo(gores=MetascapeGofile,cutoff_p = 0.05,cutoff_top = 10)
}
}

