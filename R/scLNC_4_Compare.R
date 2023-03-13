
#' Compare units from different items.
#'
#' @param Activity Boolean values determine whether to compare differential activity of units between two groups. The default value is TRUE.
#' @param objectInput A scLNC object with AUC slot.
#' @param DEitem ttribute for comparison.
#' @param DisplayUnit Boolean values determine whether to compare composition of units between two groups. The default value is TRUE.
#' @param objectCtrl A scLNC object of group1.
#' @param objectCondi A scLNC object of group2.
#' @param lncRNA Interesting lncRNA unit name.
#' @param corcutCtrl Top pairs with high correlation coefficient in group1.
#' @param corcutCondi Top pairs with high correlation coefficient in group2.
#' @param DeGO Boolean values determine whether to compare differences in functional enrichment between the two groups of units. The default value is TRUE.
#' @param CopareGOfile GO_file Metascape results.
#' @param CopareGeneList genelist Two mRNA lists from two units.
#'
#'
#' @return a scLNC object with Differential activity of units slot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject)
#' LNCobject=scLNC_3_Unit (objectInput=LNCobject,AUC=TRUE,displayLncRNA=NULL)
#' LNCobject=scLNC_4_Compare(Activity=TRUE,objectInput=LNCobject,
#' DEitem='majorCluster',DisplayUnit=FALSE,DeGO=FALSE)
#'}
#'
scLNC_4_Compare <- function(Activity=TRUE,objectInput,DEitem='majorCluster',
					DisplayUnit=TRUE,objectCtrl,objectCondi,lncRNA,corcutCtrl,corcutCondi,
					DeGO=TRUE,CopareGOfile,CopareGeneList){

if(Activity){
objectInput <- DeActivity(object=objectInput,item=DEitem,FC=0.1,pvalue=0.05,min.pct=0.3,padj=1)
return(objectInput)
}
if(DisplayUnit){
display_unit(object1=objectCtrl,object2=objectCondi,myunit=lncRNA,corcut1=corcutCtrl,corcut2=corcutCondi)

}
if(DeGO){
lnc_network(GO_file=CopareGOfile,genelist=CopareGeneList)
}

}


