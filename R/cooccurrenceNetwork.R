#' Plots a co-occurrence network visualizing the co-occurrence of taxa (at genus level) in spots. Edges are drawn if SparCC correlation > threshold.
#'
#' @param object List. Output from decontaminate() function.
#' @param threshold double. SparCC correlation threshold (value between -1 and 1) to be used to define edges between taxa on co-occurrence network. Defaults to 0.1.
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#'
#' @return two plots: first plot is an igraph plot, second plot is the corresponding legend.
#' @export
#'
#' @examples
#' \dontrun{
#' cooccurrenceNetwork(CRC_16, threshold=0.1)
#' cooccurrenceNetwork(CRC_16, threshold=0.5)
#' }
#'
cooccurrenceNetwork <- function(object, threshold=0.1,taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

  devtools::install_github("zdk123/SpiecEasi")

  #create a matrix with rows=taxa, cols=spots for sample
  s1 <- object$matrix
  #s1 <- as.matrix(object[["seurat_object"]]@assays[["TAXA_G"]]@counts)
  #remove taxa with UMI counts < 10, remove spots with no taxa counts
  s1 <- s1[Matrix::rowSums(s1)>10, Matrix::colSums(s1)>0]
  colnames(s1) <- NULL
  rownames <- taxonomizr::getTaxonomy(ids=rownames(s1), taxonomizrDB, desiredTaxa = c("genus"))
  combined_matrix <- as.matrix(s1)

  #run sparcc
  sparcc_example <- SpiecEasi::sparcc(t(combined_matrix))
  ## Define arbitrary threshold for SparCC correlation matrix for the graph
  sparcc.graph <- abs(sparcc_example$Cor) >= threshold
  diag(sparcc.graph) <- 0
  sparcc.graph <- Matrix::Matrix(sparcc.graph, sparse=TRUE)
  n_sig_taxa <- sum(Matrix::rowSums(sparcc.graph)>0)
  ## Create igraph object
  ig.sparcc <- SpiecEasi::adj2igraph(sparcc.graph)
  #only keep nodes with degree >= 1
  rownames <- rownames[igraph::degree(ig.sparcc)>=1]
  ig.sparcc <- igraph::delete.vertices(ig.sparcc, igraph::degree(ig.sparcc)<1)
  vsize    <- 4
  am.coord <- igraph::layout.fruchterman.reingold(ig.sparcc)
  #color palette with 25 colors
  c25 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black","gold1","skyblue2",
    "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon",
    "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1",
    "yellow4", "yellow3","darkorange4", "brown")
  plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, vertex.color=c25, main="taxa co-occurrence network at spot level")
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  graphics::legend("bottomleft", inset=0, fill=c25[1:n_sig_taxa], legend=rownames, border=NA)
}


