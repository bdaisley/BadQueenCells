

interactive_strip_plot <- function(df= NULL, eff= NULL, pvalue=NULL, group=NULL, labs=NULL, effect_cut=0.3, pval_cut=0.05, jitter.height= 0.2, point.size=2, text.size=0.7, min.x=NULL, max.x=NULL, plot.labs=TRUE){
  
  aldex.graphing <- df %>% mutate(col.list =case_when(
    .[,pvalue] < pval_cut & abs(.[,eff]) < effect_cut ~ "sig_p",
    abs(.[,eff]) > effect_cut ~ "sig_eff",
    #.[,pvalue] < pval_cut & abs(.[,eff]) > effect_cut ~ "sig_p_eff",
    TRUE ~ "ns")) %>%
    mutate(sig.labels = ifelse(col.list!="sig_eff" & col.list!="sig_p" , "", get(labs)))
  #min.vol <- min(df %>% dplyr::select(eff)) - 0.2
  #max.vol <- max(df %>% dplyr::select(eff)) + 0.2
  if(is.null(min.x) == TRUE){
    min.vol <- (max(abs(df %>% dplyr::select(eff))) + 0.5) * -1
    max.vol <- max(abs(df %>% dplyr::select(eff))) + 0.5
  }
  else if(is.null(min.x) == FALSE){
    min.vol <- min.x
    max.vol <- max.x
  }
  unique.rect <- length(unique(aldex.graphing[,group]))/2
  aldex.graphing[,group] <- factor(aldex.graphing[,group])
  
  #unique.rect <- 10
  unique.rect.seq1 <- seq(1.5, 100, by=2)
  unique.rect.seq2 <- seq(2.5, 100, by=2)
  unique.rect.seq1 <- unique.rect.seq1[1:unique.rect]
  unique.rect.seq2 <- unique.rect.seq2[1:unique.rect]
  unique.rect.df <- as.data.frame(cbind(paste(-Inf),paste(Inf), unique.rect.seq1, unique.rect.seq2))
  colnames(unique.rect.df)<- c("xmin", "xmax", "ymin", "ymax")
  unique.rect.df <- unique.rect.df %>% mutate_at(vars(xmin,xmax,ymin,ymax), as.numeric)
  aldex.graphing
  if(plot.labs==TRUE){
    ggplot(aldex.graphing, aes(x=get(eff), y=get(group), label=sig.labels)) +
      geom_text(size = text.size, hjust=0.3, vjust=0, check_overlap = FALSE, position=position_jitter(width=0.1,height=0.75)) +
      ggiraph::geom_point_interactive(aes(color = col.list,
                                          tooltip =  paste0(
                                            "Ortholog: ", .data[[labs]], #stringr::str_trunc(.data[[labs]], width = 70),
                                            "\nEffect size: ", round(.data[[eff]], 3),
                                            "\nP-value: ", round(.data[[pvalue]], 3)),
                                          data_id =sig.labels), 
                                      hover_nearest = TRUE,
                                      size = point.size,
                                      position = ggbeeswarm::position_quasirandom(
                                        width = jitter.height,
                                        method = "pseudorandom",
                                        orientation = "y")) +
      # ggbeeswarm::geom_quasirandom(aes(), groupOnX=FALSE, size=point.size, width= jitter.height, method = "pseudorandom") + #quasirandom, pseudorandom
      scale_x_continuous(limits=c(min.vol, max.vol), breaks=c(-6,-4,-3,-2,-1,0,1,2,3,4,6)) +
      scale_y_discrete(limits = rev(levels(aldex.graphing[,group]))) +
      geom_vline(xintercept=c(0), lwd=0.1, linetype = "longdash", colour="grey50") +
      geom_vline(xintercept=c(-1,1), lwd=0.1, linetype = "dashed", colour="grey50") +
      scale_color_manual(values=c(ns="#0000001A", sig_p="#FE4242B3", sig_eff="#3333FF"),
                         labels = c(
                           ns = "Not significant",
                           sig_p = "p<0.05",
                           sig_eff = "p<0.05 and |effect|>1"
                         )) +
      geom_rect(data=unique.rect.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha=0.10,  fill = "black", inherit.aes=FALSE) +
      xlab(eff) +
      theme(panel.border = element_rect(color ="grey", fill = NA, size = 0.5, linetype = 1), axis.line = element_line()) +
      theme(axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)) + # axis titles
      theme(axis.title.y=element_blank()) +#+ geom_point() + coord_flip()
      theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
      theme(axis.text.x=element_text(angle=0, hjust=0.5, size=8), axis.text.y=element_text(size=10))  # axis tick labels
  }
  else if(plot.labs==FALSE){
    ggplot(aldex.graphing, aes(x=get(eff), y=get(group), label=sig.labels)) +
      geom_jitter(aes(color=col.list), height = jitter.height, size=point.size) +
      scale_x_continuous(limits=c(min.vol, max.vol), breaks=c(-6,-4,-3,-2,-1,0,1,2,3,4,6)) +
      scale_y_discrete(limits = rev(levels(aldex.graphing[,group]))) +
      geom_vline(xintercept=c(0), lwd=0.1, linetype = "longdash", colour="grey") +
      geom_vline(xintercept=c(-1,1), lwd=0.1, linetype = "dashed", colour="grey") +
      scale_color_manual(values=c(ns="#0000001A", sig_p="#FE4242B3", sig_eff="#3333FF")) +
      geom_rect(data=unique.rect.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha=0.10,  fill = "black", inherit.aes=FALSE) +
      xlab(eff) +
      theme(panel.border = element_rect(color ="grey", fill = NA, size = 0.5, linetype = 1), axis.line = element_line()) +
      theme(axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)) + # axis titles
      theme(axis.title.y=element_blank()) +#+ geom_point() + coord_flip()
      theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
      theme(axis.text.x=element_text(angle=0, hjust=0.5, size=8), axis.text.y=element_blank()) # axis tick labels
  }
}
