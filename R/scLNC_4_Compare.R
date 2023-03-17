
#' Compare units from different items.
#' @import plyr
#'
#' @param DEtarget Boolean values determine whether to compare composition of units between two groups. The default value is TRUE.
#' @param ob.ls A list of objects including two scLNC objects that need to be compared.
#' @param lncRNA Interesting lncRNA unit name.
#' @param DeGO Boolean values determine whether to compare differences in functional enrichment between the two groups of units. The default value is TRUE.
#' @param CopareGOfile A GO_file Metascape results of two groups.
#' @param CopareGeneList A data frame including two mRNA lists from two groups.
#'
#' @return Corresponding plot.
#' @export
#'
#' @examples
scLNC_4_Compare <- function(DEtarget=TRUE,ob.ls,lncRNA,
                            DeGO=TRUE,CopareGOfile,CopareGeneList){


  if(DEtarget){

    ob.ls=lapply(ob.ls,function(i){

      q1=setdiff(c('withEnhancer','withPromotor','withSeq','ndG','withcyto','withTF','withPairs','score'),colnames(i@ link.data$pairs))
      df=data.frame(matrix(ncol=length(q1),nrow=nrow(i@ link.data$pairs)))
      colnames(df)=q1
      i@ link.data$pairs=cbind(i@ link.data$pairs,df)
      return(i)
    })

    display_unit(object.list=ob.ls,myunit=lncRNA)
  }

  if(DeGO){

    lnc_network(object.list=ob.ls,lncRNA,GO_file=CopareGOfile,genelist=CopareGeneList)
  }


}

