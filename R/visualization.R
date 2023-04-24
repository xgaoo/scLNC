pie_plot <- function(df, item, colorder = 1, cols = NULL, title = "Classification") {


  col.default <- c("#FEC510", "#0B775E", "#DA2828", "#000E92", "#F1BB7B", "#5B1A18", "#446455", "#FDD262", "#85D4E3", "#F4B5BD", "#6AACD0", "#E58267")
  if (is.null(cols) & colorder == 1) {
    cols <- col.default
  } else if (is.null(cols) & colorder == -1) {
    cols <- rev(col.default)
  } else if (!is.null(cols) & colorder == -1) {
    cols <- rev(cols)
  } else {
    cols <- cols
  }
  label_num <- as.character(t(table(df[, item]))[, levels(df[, item])[levels(df[, item]) %in% df[, item]]])
  p <- ggplot(df, aes(x = plot)) + geom_bar(aes(fill = df[, item])) + coord_polar(theta = "y") + labs(x = "", y = "", title = title) + theme(axis.text = element_blank(),
                                                                                                                                             axis.ticks = element_blank(), panel.grid = element_line(colour = NA), panel.background = element_rect(fill = "white", colour = NA), panel.border = element_blank(),
                                                                                                                                             legend.text = element_text(size =38, margin = margin(l = 46, unit = "pt")), legend.direction = "vertical", legend.title = element_text(size = 0,
                                                                                                                                                                                                                                                                                    face = "bold"), legend.position = "bottom", legend.key.size = unit(30, "pt"), plot.title = element_text(size = 30, hjust = 0.5))
  if (length(label_num) > length(cols)) {
    p <- p + scale_fill_manual(labels = paste(levels(df[, item]), "(", label_num, ")", sep = ""))
  } else {
    cols <- cols[1:length(label_num)]
    p <- p + scale_fill_manual(values = cols, labels = paste(levels(df[, item]), "(", label_num, ")", sep = ""))
  }

  return(p)
}


#' Pie chart of gene type statistics.
#' @import ggplot2
#'
#' @param object A scLNC object.
#' @param genetype_range The range of gene types, must be one of "all" (default), "lncRNA" or "lncRNA_overlap".
#' @return A plot of gene type statistics
#' @export
#'
#'

classification <- function(object, genetype_range = "all") {

  expMatr <- object@counts$filtered
  gtf <- object@gtf.info$info
  gtf.simple <- unique(gtf[, c("gene_ID", "gene_type", "simple_type")])
  rownames(gtf.simple) <- gtf.simple$gene_ID
  geneid <- rownames(expMatr)

  df <- data.frame(geneID = geneid, geneType = ifelse(geneid %in% gtf$gene_ID, as.character(gtf.simple[geneid, "gene_type"]), "not_annotate"), geneSimpletype = ifelse(geneid %in%
                                                                                                                                                                         gtf$gene_ID, as.character(gtf.simple[geneid, "simple_type"]), "not_annotate"), plot = genetype_range)
  df$geneSimpletype <- factor(df$geneSimpletype, ordered = TRUE, levels = c(object@gtf.info$mRNA.type, "long_non_coding", "others", "not_annotate"))
  df$geneType <- factor(df$geneType, levels = c(object@gtf.info$lncRNA.type, levels(df$geneType)[!levels(df$geneType) %in% object@gtf.info$lncRNA.type]))
  if (genetype_range == "all") {
    p <- pie_plot(df, item = "geneSimpletype", colorder = 1, title = "Classification of genes")
  } else if (genetype_range == "lncRNA") {
    p <- pie_plot(df[df$geneSimpletype == "long_non_coding", ], item = "geneType", colorder = -1, title = "Classification of lncRNAs")
  } else if (genetype_range == "lncRNA_overlap") {


    q1 <- subset(cisfile, locus_distance == 0)
    q1$lncRNA <- gsub("_.*", "", q1$relation)
    lncRNA_0 <- unique(q1$lncRNA)
    df$lncRNA_0 <- "NoneOverlap"
    df$lncRNA_0[df$geneID %in% lncRNA_0] <- "Overlap"
    df$lncRNA_0 <- factor(df$lncRNA_0, levels = unique(df$lncRNA_0))
    p <- pie_plot(df[df$geneSimpletype == "long_non_coding", ], item = "lncRNA_0", colorder = -1, title = "Classification of lncRNAs on the Genome")
  }


  return(p)
}



#' Identify highly variable genes.
#' @import ggplot2
#' @importFrom edgeR rpkm
#' @param object A scLNC object.
#' @param limx The X-axis range of the fitted curve.
#' @param labelgene A vector containing the gene symbols labeled in the plot.
#'
#' @return A curve to the relationship between gene expression and variance of genes.
#' @export
#'
#'
scCV2=function(object,limx = c(-2, 2),labelgene=''){
  #gene_transcript_length = read.table('/home/leiyang/test/pan-cancer/workdir/Item_iter3/00.data/mapGene/gene_transcript.length', sep = ' ')
  genelen = sapply(split(gene_transcript_length$V3, gene_transcript_length$V1), max)



  geneinfo=as.data.frame(rbind(cbind(object@gene.list[[1]],'mRNA'),cbind(object@gene.list[[2]],'lncRNA')))
  colnames(geneinfo)=c('gene_ID','simple_type')
  rownames(geneinfo)=geneinfo$gene_ID
  geneinfo$show.name=gtfid2genename(geneinfo$gene_ID,gtf.info = gtf_read()$info)
  genelist =rownames(geneinfo)
  gene_rpkm = rpkm(object@ counts$ filtered[genelist, ], gene.length=genelen[genelist])


  meansGenes <- rowMeans(gene_rpkm, na.rm = TRUE)
  rowVars <- function(x) {
    unlist(apply(x, 1, var, na.rm = TRUE))
  }

  varsGenes <- rowVars(gene_rpkm)
  cv2Genes <- varsGenes/meansGenes^2

  geneinfo$cv2Genes = cv2Genes
  geneinfo$meansGenes = meansGenes

  geneinfo_alll = geneinfo %>% filter(simple_type %in% c('lncRNA'))
  geneinfo_allm = geneinfo %>% filter(simple_type %in% c('mRNA'))

  geneinfo_fitl = geneinfo %>% filter(simple_type %in% c('lncRNA'), meansGenes < 10**limx[2], meansGenes > 10**limx[1])
  geneinfo_fitm = geneinfo %>% filter(simple_type %in% c('mRNA'), meansGenes < 10**limx[2], meansGenes > 10**limx[1])

  xx = data.frame(logMean = log10(geneinfo_fitl$meansGenes), logcv2 = log10(geneinfo_fitl$cv2Genes))
  fit = loess(logcv2 ~ logMean, xx)
  pred <- predict(fit, log10(geneinfo_alll$meansGenes), se=TRUE)
  tmp <- data.frame(
    fit=pred$fit,
    max = pred$fit + pred$se.fit * qt(0.99 / 2 + .5, pred$df),
    min = pred$fit - pred$se.fit * qt(0.99 / 2 + .5, pred$df)
  )
  geneinfo_alll_info = cbind(geneinfo_alll, tmp)
  geneinfo_alll_info$color = ifelse(
    !is.na(geneinfo_alll_info$fit),
    ifelse(
      log10(geneinfo_alll_info$cv2Genes) > geneinfo_alll_info$max,
      'highVar',
      ifelse(
        log10(geneinfo_alll_info$cv2Genes) < geneinfo_alll_info$min,
        'lowVar',
        'mid'
      )
    ),
    'other'
  )


  p2b=ggplot(data=NULL) +
    geom_point(data=geneinfo_allm, aes(log10(meansGenes), log10(cv2Genes)), color='#7d7878', size=0.1, alpha=0.2) +
    geom_point(data=geneinfo_alll_info, aes(log10(meansGenes), log10(cv2Genes), color=color), size=0.2) +
    geom_vline(xintercept = c(limx, 0), linetype='dashed') +
    geom_smooth(data=geneinfo_fitl, aes(log10(meansGenes), log10(cv2Genes)), method='loess', level=0.99, inherit.aes = F, color='red', size=0.5, fill=alpha('red',0.2)) +
    scale_color_manual(breaks = c('highVar', 'lowVar', 'mid', 'other', 'mRNA'), values = c('red', 'blue', 'black', '#666666', '#7d7878')) +
    geom_smooth(data=geneinfo_fitm, aes(log10(meansGenes), log10(cv2Genes)), method='loess', level=0, inherit.aes = F, color='#009900', size=0.5) +
    scale_x_continuous(breaks = seq(-5, 5, 1)) +
    geom_text_repel(
      data=geneinfo_alll_info,
      aes(x=log10(meansGenes), y=log10(cv2Genes), label=ifelse(show.name %in% labelgene, show.name, ''), color=color),
      min.segment.length = 0,
      segment.color='#666666',
      max.overlaps = 1000 )+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", fill = NA))
  return(p2b)
}



#' Identify the enrichment of cell types under different experimental conditions.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @param object A scLNC object.
#'
#' @return A heatmap of cell type enrichment.
#' @export
#'
#'
celltype_enrich=function(object=LRAT.object){

  per_dataset_Celltype=table(object@ cell.info$majorCluster,object@ cell.info$Tissue)
  per_dataset_Celltype=per_dataset_Celltype[,1:2]

  p_mx=matrix(rep(1,nrow(per_dataset_Celltype)*ncol(per_dataset_Celltype)),nrow(per_dataset_Celltype),ncol(per_dataset_Celltype))
  colnames(p_mx)=colnames(per_dataset_Celltype)
  rownames(p_mx)=rownames(per_dataset_Celltype)

  OR=matrix(rep(0,nrow(per_dataset_Celltype)*ncol(per_dataset_Celltype)),nrow(per_dataset_Celltype),ncol(per_dataset_Celltype))
  colnames(OR)=colnames(per_dataset_Celltype)
  rownames(OR)=rownames(per_dataset_Celltype)

  for(j in 1:length(rownames(per_dataset_Celltype))){
    for(i in 1:length(colnames(per_dataset_Celltype))){
      mycelltype=c(per_dataset_Celltype[j,i],sum(per_dataset_Celltype[j,-i]))
      others=c(sum(per_dataset_Celltype[-j,i]),sum(per_dataset_Celltype[-j,-i]))
      mychisq_matrix=cbind(mycelltype,others)
      chi=fisher.test(mychisq_matrix)
      p_value=chi$p.value
      p_mx[j,i]=p_value=chi$p.value
      OR[j,i]=(mychisq_matrix[1,1]*mychisq_matrix[2,2])/(mychisq_matrix[1,2]*mychisq_matrix[2,1])
    }}
  p_o=p_mx#
  p_o[p_o<0.01]="**"
  p_o[p_o<0.05&p_o!='**']="*"
  p_o[p_o>=0.05&p_o!='**'&p_o!='*']=''


  col_fun <- colorRamp2(breaks = c(0,1,2,3,4),
                        colors = c("#FFEDA0", "white","#41B6C4",'#225EA8',"#081D58"))

  TextFunc <- function(dat, col = "black", fontsize = 12, numdat = TRUE,
                       digit = 2){
    if(numdat == TRUE){
      function(j, i, x, y, width, height, fill){
        grid.text(round(dat, digit)[i, j], x, y, gp = gpar(fontsize = fontsize, col  = col))
      }
    }else{
      function(j, i, x, y, width, height, fill){
        grid.text(dat[i, j], x, y, gp = gpar(fontsize = fontsize, col  = col))
      }
    }}
  pph=Heatmap(OR, name = "Ratio", col = col_fun,
              cluster_rows=TRUE,cluster_columns=FALSE ,
              # rect_gp = gpar(col = "white", lwd = 1),
              cell_fun = TextFunc(p_o, numdat = F),row_names_side = c("left"))
  return(pph)
}


#' The difference in the number of genes across different cell types in different experimental condition.
#'
#' @param object A scLNC object.
#' @param genetype 'lncRNA' or 'mRNA'.
#' @param item A experimental condition variables to group cells by.
#' @param item.level The order of experimental condition variables to display.
#' @param split.by A variable to split the violin plots by.
#' @param disorder The order of cell types to display.
#'
#' @return A violin plot.
#' @export
#'
#'
StatGeneNum=function(object,genetype='lncRNA',item='Tissue',item.level=c('T','N','P'),split.by='majorCluster', disorder=c('C04_CD8-LAYN','C08_CD4-CTLA4','C10_CD4-CXCL13','C05_CD8-GZMK','C11_CD4-GNLY',"C09_CD4-GZMA",'C07_CD4-FOXP3','C03_CD8-SLC4A10','C01_CD8-LEF1','C06_CD4-CCR7','C02_CD8-CX3CR1')){
  object.majorTissue=splitLRAT(object,split.item=item)
  gene.num.list=lapply(names(object.majorTissue),function(x){
    df=object.majorTissue[[x]]@ counts$ filtered
    genes=object.majorTissue[[x]]@ gene.list[[genetype]]
    count=object.majorTissue[[x]]@ counts$ filtered[genes,]
    gn<-sapply(colnames(df), function(x){y<-count[,x];return(length(y[y>1]))})
    gnlist<-sapply(colnames(df), function(x){y<-count[,x];return(names(y[y>1]))})
    temp=data.frame(group=x,num=gn,item.name=gsub('.*__','',x))
  })

  gene.num1=do.call(rbind,gene.num.list)

  object@cell.info$majorCluster3=paste(object@cell.info[split.by][,1],object@cell.info[item][,1],sep="__")
  object.majorTissue=splitLRAT(object,split.item='majorCluster3')
  gene.num.list=lapply(names(object.majorTissue),function(x){
    df=object.majorTissue[[x]]@ counts$ filtered
    genes=object.majorTissue[[x]]@ gene.list[[genetype]]
    count=object.majorTissue[[x]]@ counts$ filtered[genes,]
    gn<-sapply(colnames(df), function(x){y<-count[,x];return(length(y[y>1]))})
    gnlist<-sapply(colnames(df), function(x){y<-count[,x];return(names(y[y>1]))})
    temp=data.frame(group=x,num=gn,sp.name=gsub('__.*','',x),item.name=gsub('.*__','',x))
  })

  gene.num2=do.call(rbind,gene.num.list)
  gene.num1$sp.name='all'
  gene.num=rbind(gene.num1,gene.num2)
  gene.num$sp.name=factor(gene.num$sp.name,levels=unique(gene.num$sp.name))
  gene.num$item.name=factor(gene.num$item.name,levels=item.level)



  gene.num2=subset(gene.num,sp.name!='all')
  gene.num2$sp.name=factor(gene.num2$sp.name,levels=disorder)


  p2=ggplot(gene.num2, aes(x=sp.name, y=num,fill=item.name)) +
    geom_violin(trim=FALSE,color="white") +

    geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.shape = NA)+

    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",family="Times",size=10),
          axis.text.y=element_text(family="Times",size=16,face="plain"),
          axis.title.y=element_text(family="Times",size = 20,face="plain"),
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
          legend.text=element_text(face="italic", family="Times", colour="black",
                                   size=16),
          legend.title=element_text(face="italic", family="Times", colour="black",
                                    size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ylab("num")+xlab(" ")+
  facet_wrap(~item.name,ncol=1)+labs(title = paste0('Counts of ',genetype))+
    stat_compare_means(aes(group = sp.name),ref.group = ".all.",label = "p.signif",hide.ns=TRUE, symnum.args=list(
      cutpoints = c(0,  0.001, 0.01, 0.05, 1),
      symbols = c( "***", "**", "*", "ns")),method.args = list(alternative = "two.sided"))+theme(strip.text.x = element_blank())+
    geom_hline(aes(yintercept=round(mean(num),2)),linetype="dashed",color='grey33',size=0.1)+
    guides(fill=guide_legend(title=item))+scale_fill_manual(values = colVars2$Tissue)
  return(p2)
}



#' Dotplot of gene expression.
#'
#' @param object A scLNC object.
#' @param features Interesting gene ensemble ids, such as the calculated differentially expressed lncRNAs.
#' @param item A variable to split the dotplots by.
#' @param mytitle The title of dotplot.
#' @param split.by A experimental condition variables to group cells by.
#' @param mincell.peritem Minimum number of cells per cell type.
#'
#' @return A dotplot.
#' @export
#'
#'
 DotFeatures=function(object=LRAT.object,features=unique(DE_l$gene),item="",mytitle= NULL,split.by=NULL,mincell.peritem=15){
  if(is.null(split.by)){
    bubble.df=as.matrix(object@ counts$ normalized[unique(features),])
    bubble.df=t(bubble.df)
    bubble.df=as.data.frame(scale(bubble.df))
    bubble.df$CB=rownames(bubble.df)
    object@ cell.info$CB=rownames(object@ cell.info)
    bubble.df=merge(bubble.df,object@ cell.info[,c("CB",item)],by = "CB")
    bubble.df$CB=NULL
    celltype_v=c()
    gene_v=c()
    mean_v=c()
    ratio_v=c()
    bubble.df[,item]=as.character(bubble.df[,item])
    colnames(bubble.df)[ncol(bubble.df)]='tempItem'
    for (i in unique(bubble.df$tempItem)) {
      bubble.df_small=bubble.df%>%filter(tempItem==i)
      for (j in features) {
        exp_mean=mean(bubble.df_small[,j])
        exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
        celltype_v=append(celltype_v,i)
        gene_v=append(gene_v,j)
        mean_v=append(mean_v,exp_mean)
        ratio_v=append(ratio_v,exp_ratio)  }}
    plotdf=data.frame(celltype=celltype_v,  gene=gene_v,  exp=mean_v,  ratio=ratio_v)

    plotdf$gene=factor(plotdf$gene,levels = rev(unique(as.character(features))))
    plotdf$exp=ifelse(plotdf$exp>1,1,plotdf$exp)


    plotdf$genename=gtfid2genename(id =as.character(plotdf$gene),gtf.info = gtf_read()$info)
    plotdf$genename= factor(plotdf$genename,levels=unique(plotdf$genename))
    plotdf$celltype=as.character(plotdf$celltype)

    p1=ggplot( plotdf,aes(x=celltype,y=genename,size=ratio,color=exp))+geom_point()+  scale_x_discrete("")+scale_y_discrete("")+
      scale_color_gradientn(colours = rev(c(brewer.pal(11, "Spectral")[1:11])))+  scale_size_continuous(limits = c(0,1))+theme_bw()+
      theme(legend.position="right",    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)  )+
      labs( title =mytitle)
    return(p1)
  }else{
    object@ cell.info$majorCluster2=paste(object@ cell.info[,split.by],object@ cell.info[,item],sep='&_')
    bubble.df=as.matrix(object@ counts$ normalized[unique(features),])
    bubble.df=t(bubble.df)
    bubble.df=as.data.frame(scale(bubble.df))
    bubble.df$CB=rownames(bubble.df)
    object@ cell.info$CB=rownames(object@ cell.info)
    bubble.df=merge(bubble.df,object@ cell.info[,c("CB","majorCluster2")],by = "CB")
    bubble.df$CB=NULL
    majorCluster2_v=c()
    gene_v=c()
    mean_v=c()
    ratio_v=c()
    bubble.df$majorCluster2=as.character(bubble.df$majorCluster2)
    for (i in unique(bubble.df$majorCluster2)) {
      bubble.df_small=bubble.df%>%filter(majorCluster2==i)
      for (j in features) {
        exp_mean=mean(bubble.df_small[,j])
        exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
        majorCluster2_v=append(majorCluster2_v,i)
        gene_v=append(gene_v,j)
        mean_v=append(mean_v,exp_mean)
        ratio_v=append(ratio_v,exp_ratio)  }}
    plotdf=data.frame(  majorCluster2=majorCluster2_v,  gene=gene_v,  exp=mean_v,  ratio=ratio_v)
    plotdf$gene=factor(plotdf$gene,levels = rev(unique(as.character(features))))
    plotdf$Tissue=gsub('&_.*','',plotdf$majorCluster2)
    plotdf$Tissue=factor( plotdf$Tissue,levels=c('T','N','P'))
    plotdf$celltype=gsub('.*&_','',plotdf$majorCluster2)
    plotdf$genename=gtfid2genename(id =as.character(plotdf$gene),gtf.info = gtf_read()$info)
    plotdf$genename= factor(plotdf$genename,levels=unique(plotdf$genename))
    plotdf$majorCluster2=as.character(plotdf$majorCluster2)
    plotdf$exp=ifelse(plotdf$exp>1,1,plotdf$exp)
    plotdf[plotdf$majorCluster2%in%names(table(object@ cell.info$majorCluster2)[table(object@ cell.info$majorCluster2)<mincell.peritem]),]$ratio=NA

    p1=ggplot(plotdf,aes(x=celltype,y=genename,size=ratio,color=exp))+geom_point()+  scale_x_discrete("")+scale_y_discrete("")+ facet_grid(~Tissue,scales='free',space='free')+
      scale_color_gradientn(colours = rev(c(brewer.pal(11, "Spectral")[1:11])))+  scale_size_continuous(limits = c(0,1))+theme_bw()+
      theme(legend.position="right", axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)  )+labs( title =mytitle)
    return(p1)
  }
 }



#' Probability density plot of the paired co-expression correlation coefficients of lncRNAs and their nearby mRNAs.
#' @import scales
#' @param CoexpPairs A data.frame where the columns are correlation coefficient and distance of lncRNAs and their nearby genes.
#'
#' @return A probability density plot.
#' @export
#'
#'
#'
#CoexpPairs_cis_cor=read.csv("/home/mengqq/scLNC/20221019figure/cis250kb/new_T_lm_250k_4.csv",head=TRUE,stringsAsFactors=FALSE)
#CoexpPairs=CoexpPairs_cis_cor[,c('cor','locus_distance_range')]
 CorDensity <- function(CoexpPairs){

   CoexpPairs$group=paste0("T_",CoexpPairs$locus_distance_range)
   CoexpPairs$Tissue='CoexpPairs'

   CoexpPairs$group=factor(CoexpPairs$group,levels=c("T_0-10kb",'T_10-100kb','T_100-250kb',"T_>=250kb"))


   Tpalette<- colorRampPalette(c("#EA746A" ,"white" ))
   Tcolors<- Tpalette(5)
   #show_col(Tcolors)

   sample1 <- subset(CoexpPairs,group=='T_0-10kb')$cor
   sample2 <- subset(CoexpPairs,group=='T_100-250kb')$cor
   group <- c(rep("T_0-10kb", length(sample1)), rep("T_100-250kb", length(sample2)))
   dat <- data.frame(KSD = c(sample1,sample2), group = group)
   # create ECDF of data
   cdf1 <- ecdf(sample1)
   cdf2 <- ecdf(sample2)

   # find min and max statistics to draw line between points of greatest distance
   minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1))
   x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )]
   y0 <- cdf1(x0)
   y1 <- cdf2(x0)
   y0
   y1
   pT=ggplot(CoexpPairs, aes(x = cor,color=group)) +stat_ecdf()+
     theme(plot.title = element_text(hjust = 0.5,size=15),
           axis.title = element_text(size=15),
           axis.text = element_text(size=13),
           axis.text.x = element_text(angle = 0,hjust = 1),
           legend.title = element_text(size = 15),
           legend.text = element_text(size = 12))+
     scale_color_manual(values =c(Tpalette(5)[1:4]),breaks=levels(CoexpPairs$group),labels=paste(names(table(CoexpPairs$group)),'(',table(CoexpPairs$group),')'))+
     labs(x="cor",y="cumulative fraction")+
     labs(title ='CoexpPairs')+theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.border = element_rect(colour = "black",fill = NA))+
     geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
                  linetype = "solid", color = "black")



   CoexpPairs.ls=split(CoexpPairs,CoexpPairs$group)

   x=combn(levels(CoexpPairs$group),2)
   RE=NULL
   for(i in 1:ncol(x)){
     tt=ks.test(CoexpPairs.ls[[x[1,i]]]$cor,CoexpPairs.ls[[x[2,i]]]$cor)

     temp= c(round(tt$p.value,4),round(tt$statistic,2)   )
     RE=rbind(RE,temp)
   }
   rownames(RE)=paste(x[1,],x[2,],sep='&')

   return(list(pT,RE))
 }




 #' Draw heatmap plot for units activity score.
 #'
 #' @import pheatmap
 #'
 #' @param object A scLNC object.
 #' @param items Attribute for comparison.
 #' @param features Interesting lncRNA unit name.
 #' @param mytitle The title of this heatmap.
 #'
 #' @return A heatmap.
 #' @export
 #'

 HeatmapPlot= function(object,items,features,mytitle){
   bubble.df=as.data.frame(object@unit$ AUC)
   df.sub=bubble.df
   rownames(df.sub)=gtfid2genename(id =gsub('\\(.*','',rownames(df.sub)),gtf.info = gtf_read()$info)
   anno_col<-data.frame(object@cell.info[colnames(df.sub),items])#
   anno_col<-anno_col[order(anno_col[,items[1]],anno_col[,items[2]],anno_col[,items[3]],decreasing=FALSE),] #
   df.sub=df.sub[features,rownames(anno_col)]
   range=c(-2,2)
   bk<-unique(c(seq(range[1],range[2], length=100)))

   pheatmap(df.sub,show_colnames = FALSE,
                      annotation_col = anno_col,
                      cluster_cols = FALSE,
                      cluster_rows = TRUE,
                      annotation_colors =colVars2,
                      fontsize = 8,
                      show_rownames = TRUE,scale='row',breaks = bk,main=mytitle)
 }





#' Similarities and differences of mRNAs co-expressed with lncRNAs between two experimental conditions.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import patchwork
#' @param object.list A list of objects including two scLNC objects that need to be compared.
#'
#'
#' @return A graph with combined bar and dot plots.
#' @export
#'
#'
 bar_stat_units=function(object.list){
   object1=object.list[[1]]
   object2=object.list[[2]]
   lo_T_new=object1@ link.data$ pairs
   lo_N_new=object2@ link.data$ pairs


   N_ls=split(lo_N_new$columnname,lo_N_new$rowname)
   T_ls=split(lo_T_new$columnname,lo_T_new$rowname)

   #target个数
   q1=union(names(N_ls),names(T_ls))
   q2=lapply(q1,function(x){target=length(union(N_ls[[x]],T_ls[[x]]));T=length(setdiff(T_ls[[x]],N_ls[[x]]));N=length(setdiff(N_ls[[x]],T_ls[[x]]));
   TN=length(intersect(N_ls[[x]],T_ls[[x]]));return(c(target,T,N,TN))})
   names(q2)=q1
   q3=do.call(rbind,q2)
   colnames(q3)=c('Target','T','N','T&N')
   q4=melt(q3)
   target1=subset(q4,Var2=='Target')
   target1$Tissue[target1$Var1%in%setdiff(names(T_ls),names(N_ls))]='T'
   target1$Tissue[target1$Var1%in%setdiff(names(N_ls),names(T_ls))]='N'
   target1$Tissue[target1$Var1%in%intersect(names(T_ls),names(N_ls))]='T&N'
   target1$Tissue=factor(target1$Tissue,levels=c('T&N','N','T'))

   target1=target1[order(target1$value,decreasing=TRUE),]#target1$Tissue,
   geneOrder=as.character(target1$Var1)

   RE=subset(q4,Var2!='Target')
   RE=subset(RE,value!=0)
   colnames(RE)=c( 'lncRNA', 'TargetTissue','target' )
   RE$lncRNA=factor(RE$lncRNA,levels=geneOrder)
   RE$TargetTissue=factor(RE$TargetTissue,levels=c('T','N','T&N'))


   colVar_unit=list(locus=c('#EA746A','#1DB4B8',"#A593E0"))
   names(colVar_unit$locus)=c('T','N','T&N')


   g22<- ggplot(RE, aes(x=lncRNA, y=target,fill=TargetTissue)) +
     geom_bar(position="stack", stat="identity") +
     geom_text(aes(label = target),size = 2, position = position_stack(vjust = 0.5))+
     theme(
       panel.grid.major = element_blank(), axis.ticks.x = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank(),axis.text.x = element_blank() )+scale_fill_manual(values =colVar_unit$locus)+labs(x = "", title = "")


   #locus
   n1=melt(table(lo_N_new$rowname,lo_N_new$locus2))
   t1=melt(table(lo_T_new$rowname,lo_T_new$locus2))
   colnames(n1)=c('gene','locus','value')
   colnames(t1)=c('gene','locus','value')
   n1$Tissue='N'
   t1$Tissue='T'
   t1$temp=paste(t1$locus,t1$gene,sep='+')
   n1$temp=paste(n1$locus,n1$gene,sep='+')
   tn=as.data.frame(rbind(t1,n1))
   tn=tn[which(tn$value!=0),]
   tn=subset(tn,select=-c(value))


   tn[tn$temp%in%names(table(tn$temp))[table(tn$temp)==2],]$Tissue="T&N"

   tn=unique(tn)
   tn$gene=factor(tn$gene,levels=geneOrder)

   tn$Tissue=factor(tn$Tissue,levels=c('T','N',"T&N"))

   tn$locus=gsub('cis_Antisense Head-to-head_withoutOverlap','XH',tn$locus)
   tn$locus=gsub('cis_Antisense Tail-to-tail_withoutOverlap','XT',tn$locus)
   tn$locus=gsub('cis_Sense Downstream_withoutOverlap','SD',tn$locus)
   tn$locus=gsub('cis_Sense Upstream_withoutOverlap','SU',tn$locus)


   tn$locus=factor(tn$locus,levels=c('SD','SU','XH','XT','trans'))

   g4=ggplot(data =tn, aes(gene, locus,colour = Tissue)) +
     geom_point(size=1) +
     #geom_line(aes(group = unit))+
     scale_color_manual(values =colVar_unit$locus)+#
     theme(
       axis.text.x.bottom = element_text(size=8,angle = 90,hjust = 1),
       axis.text.y = element_text(size=8),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank())+labs(x = "unit")
   p1 <- g22 + g4 + plot_layout(ncol = 1, heights = c(4,1.5))


   return(p1)
 }



