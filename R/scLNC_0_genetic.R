#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#data(scLNCgencode)
setClass("LRAT", list(counts = "list", cell.info = "data.frame",
                      proj.name = "character", gtf.info = "list",
                      lnc.info = "data.frame", link.data = "list",
                      gene.list = "list", unit = "list"))


#' Normalize Data
#' @description Normalize the count data in the object.
#' @param count Count matrix.
#' @param scale Scale the data. Default is 10000.
#'
#' @return Return a matrix with the normalize and log transformed data.
#' @export
#'
#' @examples
#' \dontrun{
#' count_filter=matrix(data = rbinom(n = 50, size = 5, prob = 0.2), nrow = 10)
#' count_norm <- lognorm(count_filter, scale = 10000)
#' count_norm
#' }
#'
lognorm <- function(count, scale = 10000) {
  norm <- apply(count, 2, function(x) {
    log1p(scale * x/sum(x))
  })
  return(norm)
}


filter_matrix <- function(countsMat, cell_info, thresh.min = 0, min.cells = 0,
                          min.pct = 0, min.count = 0,
                          min.genes = 0, item = NULL) {
  if (!is.null(item) && length(unique(cell_info[, item])) > 1) {
    df <- data.frame(item = cell_info[, item], row.names = rownames(cell_info))
    df.list <- split(df, df, drop = TRUE)

    Per.nCellsPerGene.list <- lapply(df.list, function(x) {
      rowSums(x = countsMat[, rownames(x),
                      drop = FALSE] > thresh.min)/ncol(countsMat[, rownames(x)])
    })
    Per.nCellsPerGene.df <- do.call(cbind, Per.nCellsPerGene.list)
    alpha.min <- apply(Per.nCellsPerGene.df, 1, max)
    features <- names(alpha.min[alpha.min >= min.pct])

    countsMat_filtered <- countsMat[features, ]

  } else {
    nCellsPerGene <- rowSums(countsMat > thresh.min)
    countsMat_filtered <- countsMat[which(nCellsPerGene >= min.cells), ]
    Per.nCellsPerGene <- nCellsPerGene[rownames(countsMat_filtered)]/
      ncol(countsMat_filtered)

    countsMat_filtered <- countsMat_filtered[which
                                             (Per.nCellsPerGene >= min.pct), ]
  }
  nGenesPerCell <- colSums(countsMat_filtered >= min.count)
  countsMat_filtered <- countsMat_filtered[, which(nGenesPerCell >= min.genes)]
  return(countsMat_filtered)
}


#' Create a scLNC object
#' @description Create a scLNC object from a gene expression matrix and a list of cell information.
#' The expected format of the input matrix is genes x cells and cells x item, respectively.
#' @importFrom methods new
#' @param count_table Gene count matrix.
#' @param cell_info A list of cell information.
#' @param min.cells Include features detected in at least this many cells. Default is 0.
#' @param min.genes Include cells where at least this many genes are detected. Default is 0.
#' @param proj.name Project name.
#' @param gtf Gene annotation from Gencode. Default is Gencode.v32.

#'
#' @return Return a scLNC object with genes and cells information.
#' @export
#'
#' @examples
#' \dontrun{
#' data(HCC_rawcount)
#' data(HCC_cellinfo)
#' data(scLNCgencode)
#' LNCobject <- createLRAT(count_table = HCC_rawcount, cell_info = HCC_cellinfo,
#'  min.cells = 15, min.genes = 200, proj.name = "HCC")
#'}
#'
createLRAT <- function(count_table, cell_info, min.cells = 0,min.genes = 0,
                       proj.name = "test", gtf = scLNCgencode) {
  count_filter <- filter_matrix(countsMat = count_table, cell_info = cell_info,
                                min.cells = min.cells, min.genes = min.genes)
  count_norm <- lognorm(count_filter, scale = 10000)
  count.list <- list(filtered = count_filter, normalized = count_norm)
  cell_info <- cell_info[colnames(count_filter), ]
  gtf.simple <- unique(gtf[, c("gene_ID", "gene_type", "simple_type")])
  rownames(gtf.simple) <- gtf.simple$gene_ID
  geneid <- rownames(count_filter)
  lncRNA_all_geneid <- geneid[geneid %in%
                                as.character(gtf[gtf$simple_type ==
                                                "long_non_coding", "gene_ID"])]
  mRNA_geneid <- geneid[geneid %in%
                          as.character(gtf[gtf$simple_type
                                          == "protein_coding", "gene_ID"])]

  #data(cisfile)
  q1=subset(cisfile,locus_distance==0)
  q1$lncRNA=gsub('_.*','',q1$relation)
  lncRNA_0=unique(q1$lncRNA)
  lncRNA_geneid <- setdiff(lncRNA_all_geneid,lncRNA_0)
  LRAT.obj <- new("LRAT", counts = count.list, cell.info = cell_info,
                  proj.name = proj.name, gtf.info= list(info=gtf),
                  gene.list = list(mRNA = mRNA_geneid,
                      lncRNA = lncRNA_geneid,lncRNA_all = lncRNA_all_geneid))
  return(LRAT.obj)
}



#' Splits object into a list of subsetted objects.
#'
#' @param object A scLNC object.
#' @param split.item Attribute for splitting.
#' @param mincell.peritem Include attribute where at least this many cells.
#' @param min.cells.pergene Include features detected in at least this many cells. Default is 0.
#' @param min.gene.percell Include cells where at least this many genes are detected. Default is 0.
#'
#' @return Return a list of scLNC object, each containing a subset of cells from the original object.
#' @export
#'
#' @examples
#' \dontrun{
#' data(HCC_rawcount)
#' data(HCC_cellinfo)
#' LNCobject <- createLRAT(count_table = HCC_rawcount, cell_info = HCC_cellinfo,
#'  min.cells = 15, min.genes = 200, proj.name = "HCC")
#' LNCobject.list <- splitLRAT(LNCobject, split.item="Tissue")
#'}
#'
splitLRAT<-function(object,split.item,mincell.peritem=15,min.cells.pergene=0,
                    min.gene.percell=0){
  info.list<-split(object@cell.info,as.character(object@cell.info[,split.item]))

  item.level<-unique(object@cell.info[,split.item])
  obj.list<-sapply(names(info.list),function(x){
    cells<-rownames(info.list[[x]])
    if(length(cells)>mincell.peritem){
      count<-object@counts$filtered[,cells]
      cell_info<-data.frame(info.list[[x]])
      cell_info[,split.item]<-as.character(cell_info[,split.item])
      cell_info[,split.item]<-factor(as.character(cell_info[,split.item]),
                                     levels = item.level[item.level %in%
                                        as.character(cell_info[,split.item])],
                                     ordered = TRUE)

      to<-createLRAT(min.cells = min.cells.pergene,min.genes = min.gene.percell,
                     count_table = count,cell_info = cell_info,
                     proj.name = paste(object@proj.name,x,sep="_"),
                     gtf = object@gtf.info$info)
      return(to)
    }
  })
  item.list<-as.character(unlist(sapply(names(obj.list),
                        function(x){if(!is.null(obj.list[[x]])){return(x)}})))

  return(obj.list[item.list])

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# IdConversion
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#data(scLNCgencode)
#' Gene names converts to gene ensemble ids.
#'
#' @param genelist A gene names list.
#' @param gtf.info Gene annotation from Gencode.
#'
#' @return Return a gene ensemble ids list.
#' @export
#'
#' @examples
#' \dontrun{
#' data(scLNCgencode)
#' id_convert_fromgtf(genelist=c('LINC00589','LINC00963','LINC01281','CYTOR'),
#'  gtf.info = scLNCgencode)
#'}
#'
id_convert_fromgtf <- function(genelist, gtf.info = scLNCgencode) {
  id_info <- gtf.info[, c(5, 6, 7)]
  typelv <- unique(id_info$gene_type)
  typelv <- as.character(typelv[!typelv %in% c("lncRNA",
                                               "TEC", "protein_coding")])
  id_info$gene_type <- factor(id_info$gene_type,
          levels = c("lncRNA", "TEC", "protein_coding", typelv), ordered = TRUE)
  id_info <- id_info[with(id_info, order(gene_type)), ]
  id_info <- id_info[!duplicated(id_info$gene_name), ]
  rownames(id_info) <- id_info$gene_name
  id <- data.frame(genetoconv = genelist, row.names = genelist,
                   stringsAsFactors = FALSE)
  id$gene_ID <- as.character(id$genetoconv)
  g <- genelist[genelist %in% rownames(id_info)]
  id[g, c("gene_ID")] <- as.character(id_info[g, c("gene_ID")])

  return(id$gene_ID)
}

#' Gene ensemble ids converts to gene names.
#'
#' @param id A gene ensemble ids list.
#' @param gtf.info Gene annotation from Gencode.
#'
#' @return Return a gene names list.
#' @export
#'
#' @examples
#' \dontrun{
#' data(scLNCgencode)
#' gtfid2genename(id=c('ENSG00000251191','ENSG00000204054','ENSG00000235304',
#' 'ENSG00000222041'),gtf.info = scLNCgencode)
#'}
#'
gtfid2genename<-function(id,gtf.info = scLNCgencode){
  gene.df<-unique(gtf.info[,c("gene_ID","gene_name")])
  rownames(gene.df)<-gene.df$gene_ID
  return(as.character(gene.df[id,"gene_name"]))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# lncRNA_annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' LncRNA annotation information from lncRNA-related databases.
#' @description lncRNA-related databases included disease, drug, and function related information from Lncbook v1.0, FANTOM5, LncSEA and Lnc2Cancer v3.0 .
#' @import ggplot2
#' @import RColorBrewer
#' @param features Interesting lncRNA ensemble ids, defaults to all lncRNAs.
#'
#' @return Return a statistical graph of lncRNAs about information related to disease, drug, and function.
#' @export
#'
#' @examples
#' \dontrun{
#'lncRNADatabase(features=c('ENSG00000251191','ENSG00000204054',
#''ENSG00000235304','ENSG00000222041'))
#'}
#'
lncRNADatabase=function(features=NULL){
 # data(lncRNA_anno_by_dbs)
  q4=lncRNA_anno_by_dbs
  if(!is.null(features)){
    q4=subset(q4,Ensembl_gene_id%in%features)
  }
  RE=NULL
  for(i in 1:ncol(q4)){

    a2=q4[!is.na(q4[,i]),c(1,i)]
    a2=a2[which(a2[,2]!=''),]
    a2=a2[which(a2[,2]!='NA'),]
    temp=length(unique(a2$Ensembl_gene_id))
    RE=rbind(RE,temp)
  }
  rownames(RE)=colnames(q4)
  RE=as.data.frame(RE)
  RE$annotation=rownames(RE)
  RE=RE[c('Lncbook_Disease','Lncbook_Mesh_Ontology',
          'FANTOM_associated_sample_ontology','FANTOM_Associated_trait',
          'LncSEA_drug','LncSEA_disease','Lnc2Cancer_cancer.type'),]
  RE$annotation=factor(RE$annotation,levels=c('Lncbook_Disease',
          'Lncbook_Mesh_Ontology','FANTOM_associated_sample_ontology',
          'FANTOM_Associated_trait','LncSEA_drug','LncSEA_disease',
          'Lnc2Cancer_cancer.type'))

  p=ggplot(RE,aes(y=V1,x=annotation,fill=annotation))+geom_bar( stat="identity", alpha=0.7,width = 0.8)+#position="stack",
    labs(x = "", y = "", title = "")+
    theme_bw()+
    theme(panel.grid =element_blank(),axis.title.y = element_text(hjust = 0.8),
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 0,hjust = 1))+
    geom_text(aes(label = V1), position = position_stack(vjust = 0.5),size=3)+
    scale_fill_manual(values = brewer.pal(7, "Set1"))+guides(fill=FALSE)
  #scale_color_manual(values = brewer.pal(7, "Set1"))
  return(p)
}



#' LncRNA annotation information from lncRNA-mRNA-related databases.
#' @description The lncRNA-mRNA information are from published databases, including ENCORI, LncRNA2Target v2.0, LncTarD v1.0, MONOCLdb, NPInter v4.0, RISE, RAIN v1.0, RNAInter and lncReg.
#' @import ggplot2
#' @import RColorBrewer
#' @param features Interesting lncRNA ensemble ids, defaults to all lncRNAs.
#'
#' @return Return a statistical graph of lncRNAs from lncRNA-mRNA databases.
#' @export
#'
#' @examples
#' \dontrun{
#' data(scLNCgencode)
#' lncRNA_mRNADatabase(features=c('ENSG00000251191','ENSG00000204054',
#' 'ENSG00000235304','ENSG00000222041'))
#'}
#'
lncRNA_mRNADatabase=function(features=NULL){
  #data(lm_anno_by_dbs)
  qa=lm_anno_by_dbs
  data2=qa
  if(!is.null(features)){
    data2=subset(qa,lncRNA%in%features)
  }
  data2$lncRNAnames=gtfid2genename(data2$lncRNA,gtf.info = scLNCgencode)
  data2$mRNAnames=gtfid2genename(data2$mRNA,gtf.info = scLNCgencode)

  d.ls=split(data2,data2$database)
  dn=as.data.frame(unlist(lapply(d.ls,function(x){length(unique(x$lncRNA))})))
  colnames(dn)='num'
  dn$database=rownames(dn)



  p=ggplot(dn, aes(x = database, y = num,color =database)) +
    geom_segment( aes(x = database, xend = database, y = 1, yend = num))+
    geom_point(size = 14, pch = 16, bg = 5) +
    geom_text(aes(label = num),
              position = position_stack(vjust = 1),size=2,color='black')+
    theme(plot.margin=unit(c(0.5,0.5, 1, 0),'cm'))
  return(p)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# relation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setClass("LRAT", list(counts = "list", cell.info = "data.frame",
                      proj.name = "character", gtf.info = "list",
                      lnc.info = "data.frame", link.data = "list",
                      gene.list = "list", unit = "list"))


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame( row = rownames(cormat)[row(cormat)[ut]],
              column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut] )
}

#' Calculate the co-expression correlation.
#' @importFrom bigSCale compute.network
#' @import dplyr
#' @param object A scLNC object including expression matrix slot.
#' @param fileName The name of the output file.
#'
#' @return Return a co-expression correlation matrix and a co-expression correlation list including genes, co-expression relationship pairs, gene types and locus biotypes.
#' @export
#'
#' @examples
#' \dontrun{
#' library('bigSCale')
#' data(LNCobject)
#' data(scLNCgencode)
#' rdata=CaculateCorelation(object=LNCobject,fileName='HCC')
#' }
#'#'
CaculateCorelation=function(object,fileName){
  T_data=object@ counts$ normalized
  results.ctl=compute.network(expr.data = T_data, gene.names = rownames(T_data))
  saveRDS(results.ctl,paste0(fileName,'_coexpression.rds'))

  T_coexpression_dataframe=flattenCorrMatrix(as.numeric(results.ctl$correlations))
  saveRDS(T_coexpression_dataframe,paste0(fileName,'_coexpression_dataframe.rds'))

  #lncRNA+mRNA(nooverlap)
  P=subset(T_coexpression_dataframe,row%in%c(object@ gene.list$mRNA,object@ gene.list$lncRNA)&column%in%c(object@ gene.list$mRNA,object@ gene.list$lncRNA))
  P$genetype_pair='test'
  P$genetype_pair[P$row%in%c(object@ gene.list$lncRNA)&P$column%in%c(object@ gene.list$lncRNA)]='lncRNA_lncRNA'
  P$genetype_pair[P$row%in%c(object@ gene.list$mRNA)&P$column%in%c(object@ gene.list$mRNA)]='mRNA_mRNA'
  P$genetype_pair[P$row%in%c(object@ gene.list$lncRNA)&P$column%in%c(object@ gene.list$mRNA)]='lncRNA_mRNA'
  P$genetype_pair[P$row%in%c(object@ gene.list$mRNA)&P$column%in%c(object@ gene.list$lncRNA)]='mRNA_lncRNA'
  P[which(P$genetype_pair=='mRNA_lncRNA'),]=P[which(P$genetype_pair=='mRNA_lncRNA'),] %>% select(column, row, cor, genetype_pair, everything())
  P$genetype_pair=gsub('mRNA_lncRNA','lncRNA_mRNA',P$genetype_pair)
  saveRDS(P,paste0("new_",fileName,"_3type.rds"))


  P_lm=subset(P,genetype_pair=='lncRNA_mRNA')
  #data(cisfile)
  P_lm$relation=paste(P_lm$row,P_lm$column,sep="_")
  lo_P_lm_lncm=left_join(P_lm,cisfile,c('relation'='relation'))
  lo_P_lm_lncm$locus_biotypes_rough[is.na(lo_P_lm_lncm$locus_biotypes_rough)]="trans"
  lo_P_lm_lncm$locus_biotypes_detailed[is.na(lo_P_lm_lncm$locus_biotypes_detailed)]="trans"
  lo_P_lm_lncm$locus_distance[is.na(lo_P_lm_lncm$locus_distance)]=250000
  lo_P_lm_lncm$locus_distance_range[is.na(lo_P_lm_lncm$locus_distance_range)]=">=250kb"
  lo_P_lm_lncm$locus2[is.na(lo_P_lm_lncm$locus2)]="trans"
  lo_P_lm_lncm$locus[is.na(lo_P_lm_lncm$locus)]="trans"
  lo_P_lm_lncm$rowname=gtfid2genename(id =as.character(lo_P_lm_lncm$row),gtf.info = scLNCgencode)
  lo_P_lm_lncm$columnname=gtfid2genename(id =as.character(lo_P_lm_lncm$column),gtf.info = scLNCgencode)
  saveRDS(lo_P_lm_lncm,paste0("new_",fileName,"_lm.rds"))
  return(lo_P_lm_lncm)
}





#' Filtering of co-expressed relationship pairs.
#'
#' @param data A co-expressed correlation list including genes, co-expressed relationship pairs and gene types.
#' @param corcut Choose the relationship pairs with higher correlation coefficient.
#' @param RPSL Boolean values determining if co-expressed relationship pairs should be removed the ribosome gene. Default is FALSE.
#' @param targetcut Selected lncRNAs with higher number of co-expressed mRNA. Default is 0.
#' @param object A scLNC object.
#'
#' @return Return a scLNC object with co-expressed relationship pairs slot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList)
#' data(LNCobject)
#' LNCobject=FilterPairs(data=CoExpreList, object=LNCobject, corcut=0.6,
#' RPSL=FALSE, targetcut=0)
#'}
#'
FilterPairs=function(data,corcut=0.6,RPSL=FALSE,targetcut=0,object){
  lo_P_new=subset(data,cor>corcut)
  if(RPSL==TRUE){
    lo_P_new=lo_P_new[-grep('^RP[SL].*',lo_P_new$columnname),]
  }
  targetfilter=names(table(lo_P_new$row))[table(lo_P_new$row)>=targetcut]
  lo_P_new=subset(lo_P_new,row%in%targetfilter)
  object@ link.data$pairs=lo_P_new
  object@unit$unitlist=split(object@ link.data$pairs$column,object@ link.data$pairs$row)

  return(object)
}


#' Check LncTar Tool
#'@import reticulate

checkLncTar = function(){
  stat = system("which LncTar 2> /dev/null")
  if(stat > 0){
    system("mkdir -p ~/.local/LncTar")
    system("wget -O ~/.local/LncTar/LncTar.zip http://www.cuilab.cn/lnctarapp/download")
    system("unzip -f  ~/.local/LncTar/LncTar.zip -d ~/.local/")
    system("echo '#/bin/bash\nPERL5LIB=~/.local/LncTar/ perl ~/.local/LncTar/LncTar.pl $@' > ~/.local/LncTar/LncTar && chmod 755 ~/.local/LncTar/LncTar")
    Sys.setenv('PATH'= paste('~/.local/LncTar', Sys.getenv('PATH'), sep=':'))
  }
  system("which LncTar 2> /dev/null", intern = TRUE)
}

#' Check TDF Tool
#' @import reticulate

checktdf = function(){
  stat = system("rgt-TDF --version 2> /dev/null")
  if(stat > 0){
    ### prepare conda

    if(!file.exists(miniconda_path())){
      install_miniconda()
    }
    use_miniconda()
    aaa=conda_list()
    ### tdf
    if(! 'TDF' %in% aaa$name){
      conda_install('TDF', 'python==3.9.2')
      py_install('setuptools==57.5.0', 'TDF', pip=TRUE)

      py_install('rgt==0.13.2', 'TDF', pip=TRUE)
      Sys.setenv('PATH'= paste('~/.local/share/r-miniconda/envs/TDF/bin/', Sys.getenv('PATH'), sep=':'))
      system("wget -O RGT-0.13.2.tar.gz https://github.com/CostaLab/reg-gen/archive/refs/tags/RGT-0.13.2.tar.gz && tar zxf RGT-0.13.2.tar.gz\ncd reg-gen-RGT-0.13.2 && python setup.py install --rgt-tool=TDF")
      system('perl -pi -e "s#/usr/bin/python#$(which python)#" $(which rgt-TDF)')
      system("rm -rf reg-gen-RGT-0.13.2 reg-gen-RGT-0.13.2.tar.gz")
    }else{
      Sys.setenv('PATH'= paste('~/.local/share/r-miniconda/envs/TDF/bin/', Sys.getenv('PATH'), sep=':'))
    }
  }
  if(! file.exists('~/rgtdata/hg38/genome_hg38.fa')){
    system("python ~/rgtdata/setupGenomicData.py --hg38 && rm -rf reg-gen", intern = TRUE)
  }
  system("which rgt-TDF 2> /dev/null", intern = TRUE)
  # system("rgt-TDF --version", intern = TRUE)
}

checkparallel = function(){
  stat = system("which parallel 2> /dev/null")
  if(stat > 0){
    system("mkdir -p ~/.local/parallel")
    system("wget -O ~/.local/parallel/parallel-20230122.tar.bz2 https://ftp.gnu.org/gnu/parallel/parallel-20230122.tar.bz2")
    system("cd ~/.local/parallel ; tar xjf ~/.local/parallel/parallel-20230122.tar.bz2")
    system("cd ~/.local/parallel/parallel-20230122; ./configure --prefix=$(echo ~/.local/parallel) && make && make install")
    Sys.setenv('PATH'= paste('~/.local/parallel/bin', Sys.getenv('PATH'), sep=':'))
  }
  system("which parallel 2> /dev/null", intern = TRUE)
}

#' Title Run LncTar for lncRNA-mRNA interactions.
#' @import utils
#' @importFrom biomaRt listEnsemblArchives getBM useMart
#' @param pairList The lncRNA-mRNA pairs,first column lncRNA, second column mRNA, and the form is ensembl id.
#' @param ensembl_version The version of ensemble id. Default is 105.
#'
#' @return A file of lncTar result, including lncRNA, mRNA and the normalized binding free energy (ndG).
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList)
#' path_res_lnctar=run_LncTar(pairList=CoExpreList[,1:2], ensembl_version = 105)
#' }
#'
run_LncTar = function(pairList, ensembl_version = 105){
  checkparallel()
  checkLncTar()

  randStr = paste(sample(letters, 6, replace=TRUE), collapse="")
  colnames(pairList) <- c('V1', 'V2')
  ids = unique(c(pairList$V1, pairList$V2))

  host = (listEnsemblArchives() %>% filter(version == ensembl_version))$url
  mart=useMart("ensembl",dataset="hsapiens_gene_ensembl", host=host)

  seqs = getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
               filters    = "ensembl_gene_id",
               values     = c(ids),
               mart       = mart)

  seq_longest = do.call('rbind', lapply(split(seqs, seqs$ensembl_gene_id), function(m){
    head(arrange(m, -m$transcript_length),1)
  })) %>% remove_rownames

  seqs_cdna = getBM(attributes = c("ensembl_transcript_id", "cdna"),
                    filters    = "ensembl_transcript_id",
                    values     = unique(seq_longest$ensembl_transcript_id),
                    mart       = mart)

  seq_longest = left_join(seq_longest, seqs_cdna, by=c('ensembl_transcript_id'))  %>% column_to_rownames('ensembl_gene_id')

  pairList$V1seq = seq_longest[pairList$V1, 'cdna']
  pairList$V2seq = seq_longest[pairList$V2, 'cdna']
  pairList_seq = pairList %>% dplyr::select(c('V1', 'V1seq', 'V2', 'V2seq'))
  outdir = paste0(getwd(), '/', '_LncTar_', randStr)
  dir.create(outdir)
  write.table(pairList_seq, paste0(outdir, '/lncTar_input.txt'), quote = F, sep='\t',col.names = F, row.names = F)
  system(paste0('cat ', outdir, '/lncTar_input.txt | parallel -j 50 -n 10  --pipe --cat LncTar  -p 2 -f {} -d -0.08 -s F -o ', outdir, '/res_{#}'))
  system(paste0("find ", outdir, "/ -name 'res_*' | xargs cat | grep -v ^Query > ", outdir, "/lncTar_result.txt"))
  return(paste0(outdir, "/lncTar_result.txt"))
}

#' Title Run rgt-TDF to identify the target regions related to lncRNA.
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom biomaRt listEnsemblArchives
#'
#' @param pairList The lncRNA-mRNA pairs,first column lncRNA, second column mRNA, and the form is ensemble id.
#' @param tdf_type The target regions include 'enhancer' or 'promoter'.
#' @param ensembl_version The version of ensemble id. Default is 105.
#'
#' @return A file of rgt-TDF regiontest result, including DBD binds significantly to the target regions.
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList)
#' path_res_enh=run_TDF(pairList=CoExpreList[,1:2], tdf_type = 'enhancer',
#'  ensembl_version = 105)
#' }
#'
run_TDF = function(pairList, tdf_type = 'enhancer', ensembl_version = 105){
  checkparallel()
  checktdf()


  colnames(pairList) <- c('V1', 'V2')
  ids = unique(pairList$V1)
  host = (listEnsemblArchives() %>% filter(version == ensembl_version))$url
  mart=useMart("ensembl",dataset="hsapiens_gene_ensembl", host=host)

  seqs = getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
               filters    = "ensembl_gene_id",
               values     = c(ids),
               mart       = mart)

  seq_longest = do.call('rbind', lapply(split(seqs, seqs$ensembl_gene_id), function(m){
    head(arrange(m, -m$transcript_length),1)
  })) %>% remove_rownames

  seqs_cdna = getBM(attributes = c("ensembl_transcript_id", "cdna"),
                    filters    = "ensembl_transcript_id",
                    values     = unique(seq_longest$ensembl_transcript_id),
                    mart       = mart)

  seq_longest = left_join(seq_longest, seqs_cdna, by=c('ensembl_transcript_id'))  %>% column_to_rownames('ensembl_gene_id')

  # set.seed(0)
  randStr = paste(sample(letters, 6, replace=TRUE), collapse="")

  #data(enhOut)
  enhOut = unique(enhOut[c('seqnames', 'start', 'end', 'GeneID')])
  tdfdir = paste0(getwd(), '/', '_tdf_', tdf_type, '_', randStr)
  dir.create(tdfdir,showWarnings = F)
  tmp <- lapply(split(pairList, pairList$V1), function(m){
    lnc = unique(m$V1)
    outdir = paste0(tdfdir, '/', lnc)
    dir.create(outdir, showWarnings = F)
    lncPath = paste0(outdir , '/', 'lncRNA.fa')
    ## lnc seq
    write.table(
      row.names = F,
      col.names = F,
      x = paste0('>', lnc, '\n', seq_longest[lnc, 'cdna']),
      file = lncPath, quote = F, sep = '\t'
    )
    if(tdf_type == 'promoter'){
      listfile = paste0(outdir , '/', 'targetGene.list')
      ## list
      write.table(
        row.names = F,
        col.names = F,
        x = unique(m$V2),
        file = listfile, quote = F, sep = '\t'
      )
      cmd = paste('cd', outdir, '&& rgt-TDF promotertest -r', lncPath, '-organism hg38 -de', listfile, '-rn result -o', outdir, '-pl 4000')
    }else if(tdf_type == 'enhancer'){
      listfile = paste0(outdir , '/', 'targetGene.bed')
      ## bed
      x = enhOut %>% filter(GeneID %in% m$V2)
      x = unique(x[c('seqnames', 'start', 'end')])
      # x = arrange(x[,1], x[,2], x[,3])
      write.table(
        row.names = F,
        col.names = F,
        x = x,
        file = listfile, quote = F, sep = '\t'
      )
      cmd = paste('cd', outdir, '&& rgt-TDF regiontest -r', lncPath ,'-organism hg38 -bed', listfile, '-rn result -o', outdir, '-obed -n 1000 -mp 1')
    }
  })
  ## run
  write.table(
    row.names = F,
    col.names = F,
    x = paste(tmp, collapse = '\n'),
    file = paste0(tdfdir, '/run.cmd'), quote = F, sep = '\t'
  )
  system(paste0('cat ', tdfdir, '/run.cmd | parallel -j 50 2>&1 > ', tdfdir, '/run.log'), intern = T)
  if(tdf_type == 'promoter'){
    resulthtml = 'spromoters.html'
    resultbed = 'spromoters.bed'
  }else if(tdf_type == 'enhancer'){
    resulthtml = 'starget_regions.html'
    resultbed = 'starget_regions.bed'
  }
  system(paste0("cd ",tdfdir,"; for x in */result/",resulthtml,"; do perl -e 'undef $/; $_=<>; if(m#(<tbody>.*?</tbody>)#s){$_=$1; s#[\\n\\t]##g; s#<td>#\\t#g; s#</tr>#\\n#g;s#<.*?>##g ; s#^\\t##; s#\\n\\t#\\n#g; print}' $x > ${x/.html/.bed} ; done\n"))
  system(paste0("cd ",tdfdir,"; grep -H ':' */result/", resultbed,"  | perl -pe 's#/result/", resultbed, ":#\\t#' > ",tdf_type,"_result.txt\n"))
  return(paste0(tdfdir, "/", tdf_type,"_result.txt"))
}

#' Title The results of TDF enhancer were collated and the regions were mapped to genes.
#' @import dplyr
#' @import utils
#' @param PATH Input the pathway of rgt-TDF enhancer result.
#'
#' @return The lnRNA-mRNA pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList)
#' path_res_enh=run_TDF(pairList=CoExpreList[,1:2], tdf_type = 'enhancer',
#' ensembl_version = 105)
#' load_tdf_enhancer(path_res_enh)
#' }
#'

#'
load_tdf_enhancer = function(PATH='~/test/packages/_tdf_enhancer_nydgab/enhancer_result.txt'){
  ### shi yong data tihuan in packages
  #data(enhOut)
  enhOut = unique(enhOut[c('seqnames', 'start', 'end', 'GeneID')])
  enhOut_tmp = enhOut
  enhOut_tmp$region = paste0(enhOut_tmp$seqnames, ':', enhOut_tmp$start, '-', enhOut_tmp$end)
  res = read.table(PATH) %>%
    rename(c(
      'lncRNA'='V1','region'='V3', 'TTS_count'='V5', 'TTS_norm'='V6',
      'TTS_cov'='V7', 'RankSum' = 'V8'
    ))
  res = left_join(
    res,
    enhOut_tmp[c('GeneID', 'region')] %>% rename(c('targetGene'='GeneID')),
    by=c('region')
  ) %>% select(-c('V2', 'V4'))
  return(res)
}

#' Title The results of TDF promoter were collated and the regions were mapped to genes.
#' @import dplyr
#' @import utils
#' @param PATH Input the pathway of rgt-TDF promoter result.
#'
#' @return The lnRNA-mRNA pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList)
#' path_res_pro=run_TDF(pairList=CoExpreList[,1:2], tdf_type = 'promoter',
#' ensembl_version = 105)
#' load_tdf_promoter(path_res_pro)
#' }
#'
#'
load_tdf_promoter = function(PATH='~/test/packages/_tdf_promoter_nydgab/promoter_result.txt'){
  res = read.table(PATH) %>%
    rename(c(
      'lncRNA'='V1','region'='V3', 'TTS_count'='V5', 'TTS_cov'='V6',
      'RankSum' = 'V7'
    ))
  res = left_join(
    res,
    read.table('~/rgtdata/hg38/alias_human.txt', sep='\t', quote = NULL) %>%
   select(c('V1', 'V2')) %>% rename(c('targetGene'='V1', 'Symbol'='V2')),
    by=c('V4'='Symbol')
  ) %>% select(-c('V2', 'V4'))
  return(res)
}


#' Extract the co-expression list containing the input lncRNA and the co-expressed genes.
#' @importFrom  igraph graph_from_data_frame V E all_shortest_paths
#' @param lncID A ensembl id of lncRNA.
#' @param endpoints The list of mRNAs co-expressed with this lncRNA.
#' @param corrlist The co-expressed correlation list including lncRNA-mRNA and mRNA-mRNA pairs.
#'
#' @return A co-expression list between lncRNAs and target genes, including both lncRNA-mRNA and mRNA-mRNA pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(CoExpreList_lm_mm)
#' ID = 'ENSG00000268066'
#' end = unique((corrlist %>% dplyr::filter(simple_type.x == 'lncRNA',
#' simple_type.y == 'mRNA', row == lncID))$column)
#' corrlist_shortPath = getShortPath(lncID=ID, endpoints=end,
#' corrlist=CoExpreList_lm_mm)
#'}
#'
getShortPath = function(lncID, endpoints, corrlist){


  g <- graph_from_data_frame(corrlist[, c("row", "column", "cor")] %>% filter(row %in% c(lncID, endpoints), column %in% c(lncID, endpoints)), directed = FALSE, vertices = NULL)
  E(g)$weight <- 1 - E(g)$cor

  short_path <- all_shortest_paths(g, from = V(g)[name == lncID], to = V(g)[name %in% endpoints])


  res = short_path$res[lapply(short_path$res, length) > 2]
  if(length(res) == 0){
    return(NULL)
  }
  df_res = do.call('rbind', lapply(res, function(x){
    xx = names(x)
    do.call('rbind', lapply(1:(length(xx)-1), function(i){
      data.frame(from = xx[i], to=xx[i+1])
    }))
  })) %>% unique %>% remove_rownames

  df_res = left_join(
    df_res,
    unique(rbind(
      corrlist[, c("row", "column", "cor")],
      corrlist[, c("column", "row", "cor")] %>% rename(c('column'='row', 'row'='column'))
    )),
    by=c('from'='row', 'to'='column')
  )
  return(df_res)
}

#' Draw shortPath plot
#'
#' @param corrlist_shortPath A co-expression list between lncRNAs and target genes, including both lncRNA-mRNA and mRNA-mRNA pairs.
#'
#' @importFrom  igraph graph_from_data_frame graph_attr V E layout_with_dh

draw_shortPath = function(corrlist_shortPath){
  if(is.null(corrlist_shortPath)){
    return(NULL)
  }


  lncID = unique(corrlist_shortPath$from)[!unique(corrlist_shortPath$from) %in% unique(corrlist_shortPath$to)]
  gres <- graph_from_data_frame(corrlist_shortPath, directed = FALSE, vertices = NULL)
  E(gres)$weight <- 1 - E(gres)$cor
  V(gres)$color <- ifelse(
    names(V(gres)) %in% lncID, "#EE7D6F",
    ifelse(
      names(V(gres)) %in% corrlist_shortPath[corrlist_shortPath$from==lncID, ]$to, "#9ECC68",
      '#BCA9D0'
    )
  )
  V(gres)$label <- ifelse(
    names(V(gres)) %in% lncID, names(V(gres)),
    ifelse(
      names(V(gres)) %in% corrlist_shortPath$from, names(V(gres)),
      ''
    )
  )
  V(gres)$size <- ifelse(
    names(V(gres)) %in% lncID, 10,
    ifelse(
      names(V(gres)) %in% corrlist_shortPath[corrlist_shortPath$from==lncID, ]$to, 8,
      5
    )
  )
  graph_attr(gres, "layout") <- layout_with_dh
  plot(gres, vertex.label.dist=1, edge.curved=0, vertex.frame.color=NA)
}

#' Identify the indirect correlated mRNAs with input lncRNA.
#' @import dplyr
#' @param corrlist_shortPath A co-expression list between lncRNAs and target genes, including both lncRNA-mRNA and mRNA-mRNA pairs.
#'
#' @return A vector of non-directly correlated mRNAs for lncRNA.
#' @export
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(CoExpreList_lm_mm)
#' ID = 'ENSG00000204054'
#' end = unique((corrlist %>% dplyr::filter(simple_type.x == 'lncRNA',
#' simple_type.y == 'mRNA', row == lncID))$column)
#' corrlist_shortPath = getShortPath(lncID=ID,
#' endpoints=end, corrlist=CoExpreList_lm_mm)
#' getIndirectNodes(corrlist_shortPath)
#' }
#'
getIndirectNodes = function(corrlist_shortPath){
  if(is.null(corrlist_shortPath)){
    return(NULL)
  }
  lncID = unique(corrlist_shortPath$from)[!unique(corrlist_shortPath$from) %in% unique(corrlist_shortPath$to)]
  return(unique(c(corrlist_shortPath[corrlist_shortPath$from != lncID, 'to'])))
}

#' Identify the indirect and direct correlated mRNAs with lncRNAs.
#' @import dplyr
#' @import ggplot2
#' @param corrlist The co-expressed correlation list including lncRNA-mRNA and mRNA-mRNA pairs.
#'
#' @return List of lncRNA-mRNA pairs with indirect and direct correlation and statistical diagram.
#' @export
#'
#' @examples
#' \dontrun{
#' data(CoExpreList_lm_mm)
#' direct.list=stats_indirect_pairs(corrlist=CoExpreList_lm_mm)
#'}
#'
stats_indirect_pairs = function(corrlist){
  lnc_mRNA_pairs = unique((corrlist %>% filter(simple_type.x == 'lncRNA', simple_type.y == 'mRNA')))[c('row', 'column', 'cor')] %>% remove_rownames
  lnc_indirectNodes = lapply(unique(lnc_mRNA_pairs$row), function(lncID){
    endpoints          = unique((corrlist %>% dplyr::filter(simple_type.x == 'lncRNA', simple_type.y == 'mRNA', row == lncID))$column)
    corrlist_shortPath = getShortPath(lncID, endpoints, corrlist)
    indirectNodes      = getIndirectNodes(corrlist_shortPath)
    indirectNodes
  })
  names(lnc_indirectNodes) = unique(lnc_mRNA_pairs$row)

  x = lnc_indirectNodes[lapply(lnc_indirectNodes, length) > 0]
  indirect_lnc_mRNA_pairs = do.call('rbind', lapply(names(x), function(i){
    data.frame(lncRNA=i, mRNA=x[[i]], corr_type='indirect')
  }))

  lnc_mRNA_pairs_corrType = left_join(
    lnc_mRNA_pairs,
    indirect_lnc_mRNA_pairs,
    by=c('row'='lncRNA', 'column'='mRNA'),
  )
  lnc_mRNA_pairs_corrType[is.na(lnc_mRNA_pairs_corrType$corr_type), 'corr_type'] = 'direct'

  ggdata = lnc_mRNA_pairs_corrType %>% select(c('row', 'corr_type')) %>% table %>% as.data.frame(stringsAsFactors = F) %>%
    pivot_wider(names_from = 'corr_type', values_from = 'Freq') %>%
    filter(indirect > 0) %>%
    pivot_longer(names_to = 'corr_type', values_to = 'Freq', cols = c('direct', 'indirect'))

  xturns = ggdata %>% pivot_wider(names_from = 'corr_type', values_from = 'Freq') %>% column_to_rownames('row') %>% rowSums() %>% sort(decreasing = T) %>% names


  p = ggplot(ggdata) + geom_col(aes(x=row, y=Freq, fill=corr_type)) + xlim(xturns) + theme(axis.ticks = element_blank(),
                                                                                           panel.background = element_blank(),
                                                                                           axis.text.x = element_text(angle = 30,vjust = 1.5,hjust = 1),
  ) + scale_fill_manual(values=c("#a1dc65", "#c8adf8")) + labs(x=NULL, y='mRNA counts', fill=NULL)
  return(list(p=p,lnc_mRNA_pairs_corrType=lnc_mRNA_pairs_corrType))
}



#' Annotation co-expressed relationship pairs.
#' @import dplyr
#' @import utils
#' @param object A scLNC object.
#' @param TF Boolean values determining if mRNA should be annotated the TF information. Default is TRUE.
#' @param cytokine Boolean values determining if mRNA should be annotated the cytokine information. Default is TRUE.
#' @param LMpairs Boolean values determining if co-expressed relationship pairs should be annotated the lncRNA-mRNA database information. Default is TRUE.
#' @param Seq Boolean values determining if co-expressed relationship pairs should be annotated the sequence-based pairs information. Default is FALSE.
#' @param Enhancer Boolean values determining if co-expressed relationship pairs should be annotated the enhancer information. Default is FALSE.
#' @param Promoter Boolean values determining if co-expressed relationship pairs should be annotated the promoter information. Default is FALSE.
#' @param direct Boolean values determine whether direct or indirect coexpression relation pairs are distinguished. The default value is FALSE.
#' @param mycorrlist  The co-expressed correlation list including lncRNA-mRNA and mRNA-mRNA pairs.
#'
#' @return Return A scLNC object with annotational pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject)
#' scLNC <- PairsAnnotation(object=LNCobject,TF=TRUE,cytokine=TRUE,LMpairs=TRUE)
#'}
#'
PairsAnnotation=function(object,Seq=FALSE,Enhancer=FALSE, Promoter=FALSE,
              TF=TRUE,cytokine=TRUE,LMpairs=TRUE,direct=FALSE,mycorrlist=NULL){

  if(Seq==TRUE){
    path_res_lnctar = run_LncTar(object@ link.data$pairs[,c('row','column')])
    seq_all=read.table(path_res_lnctar,sep='\t',header=FALSE)
    seq_all=seq_all[,c(1,3,6)]
    colnames(seq_all)=c('Query','Target','ndG')
    seq_all$relation=paste(seq_all$Query,seq_all$Target,sep='_')
    object@ link.data$pairs=left_join(object@ link.data$pairs,seq_all[,c('relation','ndG')],by='relation')
    object@ link.data$pairs$withSeq='true'
    object@ link.data$pairs$withSeq[is.na(object@ link.data$pairs$ndG)]='false'
  }

  if(Enhancer==TRUE){
    path_res_enh = run_TDF(pairList=object@ link.data$pairs[,c('row','column')], tdf_type = 'enhancer')
    enhancer=load_tdf_enhancer(path_res_enh)
    enhancer=unique(enhancer[,c('lncRNA','targetGene')])
    colnames(enhancer)=c('V1','V3')
    enhancer$relation=paste(enhancer$V1,enhancer$V3,sep='_')

    object@ link.data$ pairs=left_join(object@ link.data$ pairs,enhancer,by='relation')
    object@ link.data$ pairs$withEnhancer='true'
    object@ link.data$ pairs$withEnhancer[is.na(object@ link.data$ pairs$V1)]='false'
  }

  if(Promoter==TRUE){
    path_res_pro = run_TDF(pairList=object@ link.data$pairs[,c('row','column')], tdf_type = 'promoter')
    promoter=load_tdf_promoter(path_res_pro)
    promoter=unique(promoter[,c('lncRNA','targetGene')])
    colnames(promoter)=c('V1','V3')
    promoter$relation=paste(promoter$V1,promoter$V3,sep='_')

    object@ link.data$ pairs=left_join(object@ link.data$ pairs,promoter,by='relation')
    object@ link.data$ pairs$withPromotor='true'
    object@ link.data$ pairs$withPromotor[is.na(object@ link.data$ pairs$V1)]='false'

  }


  if(cytokine==TRUE){
    #data(KE_cyto)
    object@ link.data$pairs=left_join(object@ link.data$pairs,KE_cyto[,c('gene','cytokine_anno')],c('columnname'='gene'))
    object@ link.data$pairs$withcyto='true'
    object@ link.data$pairs$withcyto[is.na(object@ link.data$pairs$cytokine_anno)]='false'
  }
  if(TF==TRUE){

    object@ link.data$pairs=left_join(object@ link.data$pairs,KE_TF[,c('gene','TF_anno')],c('columnname'='gene'))
    object@ link.data$pairs$withTF='true'
    object@ link.data$pairs$withTF[is.na(object@ link.data$pairs$TF_anno)]='false'
  }

  if(LMpairs==TRUE){

    #data(Alllinks_database)
    Database_lm=Alllinks_database
    object@ link.data$pairs=left_join(object@ link.data$pairs,Database_lm,by='relation')
    object@ link.data$pairs$withPairs='true'
    object@ link.data$pairs$withPairs[is.na(object@ link.data$pairs$score)]='false'
  }

  if(direct){
    direct.list=stats_indirect_pairs(corrlist=mycorrlist)
    temp=direct.list[['lnc_mRNA_pairs_corrType']]

    object@ link.data$pairs$tm=paste(object@ link.data$pairs$row,object@ link.data$pairs$column,sep='_')
    temp$tm=paste(temp$row,temp$column,sep='_')
    object@ link.data$pairs=left_join(object@ link.data$pairs,temp[,c('tm','corr_type')],by=tm)

  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# unit
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculates lncRNA units activity score.
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @param object A scLNC object.
#' @param lnclist lncRNA list. Only keep a subset of lncRNA units, defaults to all lncRNA units.
#' @param min.target Selected units with higher number of mRNA. Default is 0.
#'
#' @return Return a scLNC object with AUC slot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject2)
#' LNCobject <- AUCell_score(object = LNCobject2,lnclist = NULL,min.target = 0)
#' }
#'
AUCell_score <- function(object,lnclist=NULL,min.target=0){
  unit<-object@ link.data$pairs
  colnames(unit)[1:2]=c('lncRNA','mRNA')
  countsMat_filtered<-as.matrix(object@counts$normalized)
  unit$lncRNA <- as.character(unit$lncRNA)
  unit_list <- split(unit[,"mRNA"], unit$lncRNA)
  unit_list <- unit_list[lengths(unit_list)>=min.target]
  if(!is.null(lnclist)){
    unit_list<-unit_list[lnclist]
  }
  unit<-subset(unit,lncRNA%in% names(unit_list) & mRNA %in% unlist(unit_list))
  inputMat<-countsMat_filtered
  set.seed(1)
  aucellRankings <- AUCell_buildRankings(inputMat, plotStats=FALSE)
  unitAUC <- AUCell_calcAUC(unit_list, aucellRankings,
                            aucMaxRank=aucellRankings@nGenesDetected["1%"])
  object@unit$AUC<-unitAUC@assays@data@ listData$AUC
  return(object)
}





#' Visualization of functional enrichment results of metascape multiple gene lists.
#' @import pheatmap
#' @import utils
#' @param gores Metascape results.
#' @param cutoff_p P-value cut off. Default is 0.05.
#' @param cutoff_top Top markers cut off.
#'
#' @return Return a dotplot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(MetascapeGoMultG)
#' CompareGo(gores=MetascapeGoMultG,cutoff_p = 0.05,cutoff_top = 10)
#' }
#'
CompareGo=function(gores,cutoff_p = 0.05,cutoff_top = 10){
  gores = gores[gores$Category == 'GO Biological Processes', ]
  gores = do.call('rbind', lapply(split(gores, gores$GROUP_ID), function(x){head(arrange(x, x$Log.q.value.),1)}))

  gores = head(gores,20)
  tmp = gores[colnames(gores)[grepl('X_LogP_', colnames(gores))]]
  rownames(tmp) = paste0(gores$GO, ": ", gores$Description)
  colnames(tmp) = gsub('X_LogP_', '', colnames(tmp))

  tmp = tmp %>% rownames_to_column(var = 'Term') %>% pivot_longer(names_to = 'geneset', values_to = 'log10P', cols = colnames(tmp)) %>% filter(log10P < log10(cutoff_p))
  #tmp = do.call('rbind', lapply(split(tmp, tmp$geneset), function(m){ head(arrange(m, m$log10P), cutoff_top) }))
  tmp$Term=gsub('.*:','',tmp$Term)
  tmp_mtx = tmp %>% pivot_wider(names_from = 'geneset', values_from = log10P, values_fill = 0) %>% remove_rownames %>% column_to_rownames(var = 'Term')

  p = pheatmap::pheatmap(-tmp_mtx, silent=TRUE)
  options(repr.plot.width=8, repr.plot.height = 5+as.integer(length(unique(tmp$Term)) / 10))

  genelist=rev(unique(tmp$geneset))
  #
  p2 = ggplot(tmp) +
    geom_point(aes(factor(geneset, levels=genelist), Term, size=-log10P, color =-log10P)) +
    scale_color_gradient(low = 'blue', high = 'red') +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + scale_y_discrete(position = "right", limits=p$tree_row$labels[p$tree_row$order]) +
    labs(title='C04', x=NULL)
  return(p2)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# unit_comparision
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' The sames and differences of mRNAs coexpressed with lncRNA in an unit between two groups.
#' @description The network diagram of units between the two groups included the same lncRNA and differences mRNAs co-expressed with them.
#' @importFrom  igraph graph_from_data_frame add_shape V E
#'
#' @param object.list A list of objects including two scLNC objects that need to be compared.
#' @param myunit Interesting lncRNA unit name.
#' @param corcut1 Top pairs with high correlation coefficient in group1.
#' @param corcut2 Top pairs with high correlation coefficient in group2.
#'
#'
#' @return Return a network diagram.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject2)
#' display_unit(object.list=list(LNCobject2,LNCobject2),
#' myunit='LINC00861',corcut1=0.7,corcut2=0.5)
#'}
#'
display_unit=function(object.list,myunit='LINC010273',corcut1=0.7,corcut2=0.5){
  object1=object.list[[1]]
  object2=object.list[[2]]
  object1@ link.data$pairs=subset(object1@ link.data$pairs,cor>corcut1)
  object2@ link.data$pairs=subset(object2@ link.data$pairs,cor>corcut2)
  N_unit1_BMS1P4=object1@ link.data$pairs[object1@ link.data$pairs$rowname==myunit,][,c('rowname','columnname','withEnhancer','withPromotor','withSeq','ndG','withcyto','withTF','withPairs','score')]
  T_unit1_BMS1P4=object2@ link.data$pairs[object2@ link.data$pairs$rowname==myunit,][,c('rowname','columnname','withEnhancer','withPromotor','withSeq','ndG','withcyto','withTF','withPairs','score')]
  if(nrow(N_unit1_BMS1P4)>0){
    N_unit1_BMS1P4$Tissue='N'
    temp_inter=N_unit1_BMS1P4}
  if(nrow(T_unit1_BMS1P4)>0){
    T_unit1_BMS1P4$Tissue='T'
    temp_inter=T_unit1_BMS1P4}
  if(nrow(N_unit1_BMS1P4)>0 & nrow(T_unit1_BMS1P4)>0){
    temp_inter=rbind(N_unit1_BMS1P4,T_unit1_BMS1P4)#!!!!!!!!!!
  }
  temp_inter$Tissue[temp_inter$columnname%in%names(table(temp_inter$columnname))[table(temp_inter$columnname)==2]]="T&N"
  temp_inter=base::unique(temp_inter)


  if(!(nrow(unique(temp_inter[,c('withEnhancer','withPromotor','withSeq','withcyto','withTF','withPairs')]))==1&&
       unique(temp_inter[,c('withEnhancer')])=='false'&&unique(temp_inter[,c('withPromotor')])=='false'&&unique(temp_inter[,c('withSeq')])=='false'&&unique(temp_inter[,c('withcyto')])=='false'&&unique(temp_inter[,c('withTF')])=='false'&&unique(temp_inter[,c('withPairs')])=='false')){
    temp_add=NULL
    temp_reduce=NULL
    temp_index=NULL
    for(i in 1:nrow(temp_inter)){
      if(!(length(unique(as.character(temp_inter[i,c('withEnhancer','withPromotor','withSeq','withcyto','withTF','withPairs')]))=='false')==1&&unique(as.character(temp_inter[i,c('withEnhancer','withPromotor','withSeq','withcyto','withTF','withPairs')]))=='false')){


        temp_reduce=rbind(temp_reduce,temp_inter[i,])
        temp_index=c(temp_index,i)


        if(temp_inter[i,'withEnhancer']=='true'){
          temp_add_new=temp_inter[i,]

          temp_add_new['withPromotor']='false'
          temp_add_new['withSeq']='false'
          temp_add_new['withcyto']='false'
          temp_add_new['withTF']='false'
          temp_add_new['withPairs']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }
        if(temp_inter[i,'withPromotor']=='true'){
          temp_add_new=temp_inter[i,]
          temp_add_new['withEnhancer']='false'

          temp_add_new['withSeq']='false'
          temp_add_new['withcyto']='false'
          temp_add_new['withTF']='false'
          temp_add_new['withPairs']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }


        if(temp_inter[i,'withSeq']=='true'){
          temp_add_new=temp_inter[i,]
          temp_add_new['withEnhancer']='false'
          temp_add_new['withPromotor']='false'
          temp_add_new['withcyto']='false'
          temp_add_new['withTF']='false'
          temp_add_new['withPairs']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }
        if(temp_inter[i,'withcyto']=='true'){
          temp_add_new=temp_inter[i,]
          temp_add_new['withEnhancer']='false'
          temp_add_new['withPromotor']='false'
          temp_add_new['withSeq']='false'
          temp_add_new['withTF']='false'
          temp_add_new['withPairs']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }
        if(temp_inter[i,'withTF']=='true'){
          temp_add_new=temp_inter[i,]
          temp_add_new['withEnhancer']='false'
          temp_add_new['withPromotor']='false'
          temp_add_new['withSeq']='false'
          temp_add_new['withcyto']='false'
          temp_add_new['withPairs']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }
        if(temp_inter[i,'withPairs']=='true'){
          temp_add_new=temp_inter[i,]
          temp_add_new['withEnhancer']='false'
          temp_add_new['withPromotor']='false'
          temp_add_new['withSeq']='false'
          temp_add_new['withcyto']='false'
          temp_add_new['withTF']='false'
          temp_add=rbind(temp_add,temp_add_new)
        }
      }

    }
    temp_inter=rbind(temp_inter[-temp_index,],temp_add)

  }
  temp_vertices=unique(rbind(c(unique(temp_inter[,1]),'lncRNA'),temp_inter[,c('columnname','Tissue')]))
  temp_inter_edge=temp_inter[,c('rowname','columnname')]


  mygraph<-graph_from_data_frame(
    d=temp_inter_edge,vertices=temp_vertices,
    directed=FALSE)




  vcolor<-list(Tissue=c("lncRNA"=alpha("white",0.7),"T"=alpha('#EA746A',0.5),"T&N"=alpha('#A593E0',0.5),"N"=alpha('#1DB4B8',0.5)))

  V(mygraph)$color<-vcolor$Tissue[V(mygraph)$Tissue]
  V(mygraph)$size =0.1



  E(mygraph)$color <- alpha('grey50',0.5)

  E(mygraph)[which(temp_inter$withEnhancer=='true')]$color='purple'
  E(mygraph)[which(temp_inter$withPromotor=='true')]$color='#5f7dd8'

  E(mygraph)[which(temp_inter$withSeq=='true')]$color='orange'

  E(mygraph)[which(temp_inter$withPairs=='true')]$color='green'

  E(mygraph)$width<-1


  de=rep(1,nrow(temp_inter)+1)
  for (x in 1:nrow(temp_inter)){
    if(all(temp_inter[x,c('withEnhancer','withPromotor','withSeq','withcyto','withTF','withPairs')]=='false')){
      de[x+1]=1
      V(mygraph)$name[x+1]=""
    }}
  de[1]=2

  mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }

    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
  }

  add_shape("triangle", clip=shapes("circle")$clip,
            plot=mytriangle)

  shapes=rep('circle',length(V(mygraph)$name))
  if(length(unique(temp_inter$'withcyto'))!=1){
    shapes[grep(paste(unique(temp_inter$columnname[which(temp_inter$'withcyto'=='true')]),collapse = '|'),V(mygraph)$name)]='square'
  }
  if(length(unique(temp_inter$'withTF'))!=1){
    shapes[grep(paste(unique(temp_inter$columnname[which(temp_inter$'withTF'=='true')]),collapse = '|'),V(mygraph)$name)]="triangle"
  }

pdf(paste0('display_unit_',myunit,'.pdf'))
  plot(mygraph,vertex.size=15*de,vertex.shape=shapes,
       vertex.label.cex=.6,vertex.label.dist=0,
       edge.curved=0,vertex.label.color='black',vertex.frame.color='grey')
  legend( x=1,y=-0.8, c("T","N",'T&N'), pch=21,
          col="grey", pt.bg=c(alpha('#EA746A',0.5), alpha('#1DB4B8',0.5),alpha('#A593E0',0.5)), bty="n",
          pt.cex=2, cex=.8,  ncol=1,title="Tissue")

  legend(x=1,y=-0.3, inset=.05, title="Pairs", c("Enhancer","Promotor",'Seq','Database'),bty="n",
         cex=.8,lty=c(1),  col=c("purple", "#5f7dd8",'orange','green'))
  legend( x=1,y=-1.2, c('TF','cytokine'), pch=c(2,0),
          col="grey", pt.bg=c(alpha('#C65146',0.5), alpha('#4F86C6',0.5),alpha('#A593E0',0.5)), bty="n",
          pt.cex=2, cex=.8,  ncol=1,title="Target annotation")
  dev.off()
}



#' Differential activity of units between two groups.
#' @import Seurat
#'
#' @param object A scLNC object.
#' @param item Attribute for comparison.
#' @param FC Limit testing to units which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.1.
#' @param pvalue P-value cut off. Default is 0.05.
#' @param padj P-value adjust cut off. Default is 1.
#' @param min.pct Only test units that are obtained activity score in a minimum fraction of min.pct cells in either of the two populations. Default is 0.3.
#' @param ... Other options used to control matching behavior between duplicate strings. Passed on to [Seurat::FindAllMarkers()].
#'
#'
#' @return Return a scLNC object with Differential activity of units slot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(LNCobject2)
#' LNCobject <- AUCell_score(object = LNCobject2,lnclist = NULL,min.target = 0)
#' LNCobject<- DeActivity(object=LNCobject,item='majorCluster',
#' FC=0.1,pvalue=0.05,min.pct=0.3,padj=1)
#'}
DeActivity=function(object,item,FC=0.1,pvalue=0.05,min.pct=0.3,padj=1,...){
  data.matrix=object@unit$ AUC
  names(dimnames(data.matrix))=NULL
  rownames(data.matrix)=gsub('\\(.*','',rownames(data.matrix))

  Seu=CreateSeuratObject(counts = as.matrix(data.matrix),meta.data = object@cell.info)
  Idents(Seu)=as.factor(Seu@meta.data[item][,1])

  markers = FindAllMarkers(Seu, assay = 'RNA', logfc.threshold = 0,only.pos=TRUE,base=exp(1),...)
  markers$genename=gtfid2genename(id =markers$gene,gtf.info = scLNCgencode)

  DE=subset(markers,avg_logFC>FC)
  DE=subset(DE,p_val<pvalue)
  DE=subset(DE,p_val_adj<padj)
  DE$minpct=apply(DE,1,function(x){max(x['pct.1'],x['pct.2'])})

  DE=subset(DE,minpct>min.pct)
  object@ gene.list$DEAUC=DE
  return(object)
}




#' Differences in functional enrichment between the two groups of units.
#' @importFrom stringr str_detect str_extract
#' @importFrom dplyr distinct
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend par symbols
#' @importFrom igraph graph_from_data_frame V E layout_nicely
#'
#' @param GO_file Metascape results.
#' @param genelist Two mRNA lists from two units.
#' @param lncRNA Interesting lncRNA name.
#'
#' @return A network diagram, including differences and same in biological functions and mRNAs.
#' @export
#'
#' @examples
#' \dontrun{
#' data(Go2Group_LINC00861)
#' data(Gene2Group_LINC00861)
#' lnc_network(GO_file=Go2Group_LINC00861,genelist=Gene2Group_LINC00861,lncRNA='LINC00861')
#' }
#'
lnc_network<-function(GO_file,genelist,lncRNA){
  group1<-colnames(genelist)[1]
  group2<-colnames(genelist)[2]
  colnames(GO_file)[which(str_detect(colnames(GO_file),paste0("MEMBER_",group1)))]<-"MEMBER_group1"
  colnames(GO_file)[which(str_detect(colnames(GO_file),paste0("MEMBER_",group2)))]<-"MEMBER_group2"
  colnames(GO_file)[which(str_detect(colnames(GO_file),paste0("LogP_",group1)))]<-"LogP_group1"
  colnames(GO_file)[which(str_detect(colnames(GO_file),paste0("LogP_",group2)))]<-"LogP_group2"
  gores = subset(GO_file,Category == "GO Biological Processes")
  gores=subset(gores,gores$Log.q.value.>'-1')
  gores$GO = paste0(gores$GO, ": [", gores$Description, "]")
  goress<-subset(gores,MEMBER_group1==1&MEMBER_group2==0)
  goress<-goress[order(goress$LogP_group1),]
  group1_top5<-goress[c(1:5),]
  goresss<-subset(gores,MEMBER_group1==0&MEMBER_group2==1)
  goresss<-goresss[order(goresss$LogP_group2),]
  group2_top5<-goresss[c(1:5),]
  Go<-subset(gores,MEMBER_group1==1&MEMBER_group2==1)
  gore<-rbind(group1_top5,group2_top5,Go)
  gore<-gore[apply(gore, 1, function(x) !all(is.na(x))),]
  mRNA= do.call('rbind', apply(gore[c('Hits', 'GO','MEMBER_group1','MEMBER_group2')], 1, function(m){
    data.frame(
      mRNA = unlist(strsplit(m[['Hits']], '\\|')),
      GO = m[['GO']],
      group1 = m[['MEMBER_group1']],
      group2 = m[['MEMBER_group2']]
    )
  }))

  new_group<-data.frame(GO=unique(mRNA$GO),new_group=NA)

  for (i in 1:nrow(new_group)){
    sub<-subset(mRNA,GO==new_group$GO[i])
    if (all(sub$group1=='0')){
      new_group$new_group[i]<-group1
    } else if (all(sub$group2=='0')) {
      new_group$new_group[i]<-group2
    } else {
      new_group$new_group[i]<-paste(group1,group2)
    }
  }
  mRNA<-merge(mRNA,new_group,by='GO',all.x=T,sort=F)

  ####分割线
  data2<-data.frame(mRNA=unique(mRNA$mRNA),group=NA)
  data2$group[data2$mRNA %in% genelist[[1]]]<-group1
  data2$group[data2$mRNA %in% genelist[[2]]]<-group2
  data2$group[data2$mRNA %in% intersect(genelist[[1]],genelist[[2]])]<-paste(group1,group2)
  data2$group1<-1
  data2$group2<-1
  data2$group1[which(str_detect(data2$group, group2))] <- "0"
  data2$group2[which(str_detect(data2$group, group1))] <- "0"
  data2$group1[which(str_detect(data2$group, paste(group1,group2)))] <- "1"
  data2$group2[which(str_detect(data2$group, paste(group1,group2)))] <- "1"

  ##node and edge
  node1<-distinct(data.frame(name=mRNA$GO,lc=mRNA$new_group,size=10,label=str_extract(mRNA$GO,"\\[.+?\\]")))

  node1$lc[which(node1$lc==group2)]<-'#EA746A'
  node1$lc[which(node1$lc==group1)]<-'#1DB4B8'
  node1$lc[which(node1$lc==paste(group1,group2))]<-'black'

  node2<-data.frame(name=data2$mRNA,lc=data2$group,size=2)
  node2$lc[which(node2$lc==group2)]<-'#EA746A'
  node2$lc[which(node2$lc==group1)]<-'#1DB4B8'
  node2$lc[which(node2$lc==paste(group1,group2))]<-'black'
  node2$label<-''
  node2$label[node2$name %in% KE_cyto$gene]<-'cyto'
  node2$label[node2$name %in% KE_TF$gene]<-'TF'
  for (i in 1:nrow(node2)) {
    if (node2$label[i]=='TF') {
      node2$label[i]=node2$name[i]
    }
    if (node2$label[i]=='cyto') {
      node2$label[i]=node2$name[i]
    }
  }

  node<-rbind(node1,node2)

  node$dist<-1
  node$dist[which(str_detect(node$size, "2"))] <- 0.3


  edge<-data.frame(from=mRNA$mRNA,to=mRNA$GO,color='gray30')


  network <- graph_from_data_frame(
    vertices = node,
    d = edge,
    directed = FALSE
  )


  ###correct_pie_values
  mrna_go<-mRNA[,c(1,3,4)]
  mrna_go<-unique(mrna_go)
  mrna_go1<-t(mrna_go)
  mrna_go1<-as.data.frame(mrna_go1)
  colnames(mrna_go1)=mrna_go1[1,]
  mrna_go1<-mrna_go1[-1,]

  mrna_go2=as.data.frame(lapply(mrna_go1,as.numeric))

  go_values<-as.list(mrna_go2)

  mrna_group <-data2[,-2]
  mrna_group1<-t(mrna_group)
  mrna_group1<-as.data.frame(mrna_group1)
  colnames(mrna_group1)=mrna_group1[1,]
  mrna_group1<-mrna_group1[-1,]
  mrna_group2=as.data.frame(lapply(mrna_group1,as.numeric))
  mrna_values<-as.list(mrna_group2)
  values<-c(go_values,mrna_values)

  V(network)$size<-node$size
  V(network)$pie.color=list(c("#1DB4B8", "#EA746A"))
  V(network)$label<-node$label
  V(network)$label.color<-node$lc
  E(network)$color<-edge$color
  ###try pie vertex

  pdf(paste0(lncRNA,'_DEGO.pdf'),width=20,height=20)
      plot(network,vertex.shape="pie",vertex.frame.color = NA,vertex.pie=values,vertex.label.dist=node$dist, vertex.label.degree=-pi/2, ,label.cex=0.5,
           layout =layout_nicely,
           edge.arrow.size=0.3,
           edge.curved=0.3)
      legend("bottomright", legend = c(group2,group1), col = c('#EA746A','#1DB4B8'),bty = "n", pch=19 , pt.cex = 1,
             cex = 1, text.col="black" , horiz = FALSE,
             inset = c(0.1, 0.1))
      dev.off()
}























