#' Creates a microbiome relative abundance bar plot of the decontaminated sample (at genus level, showing the top 10 taxa).
#'
#' @param sampleName String. Sample name.
#' @param object List. Output from decontaminate() function.
#'
#' @return ggplot2 plot.
#' @export
#'
#' @examples
#' pseudoBulkProfile(sampleName="CRC_16", object=CRC_16)
pseudoBulkProfile <- function(sampleName, object){

  ##check parameters
  if (is.null(sampleName)){
    stop("please provide a sample name in sampleName")
  }

  #selecting the top10 taxa
  #checking if at least 10 taxa are present
  if (nrow(object$taxid_counts)<10){
    top = nrow(object$taxid_counts)-1
  } else {
    top=10
  }
  top_10 <- utils::head(object$taxid_counts, n=top) %>% dplyr::select(counts, genus)
  others <- object$taxid_counts[(top+1):nrow(object$taxid_counts),]
  add <- data.frame(counts=sum(others$counts), genus="others")
  top_10 <- rbind(top_10, add)
  total <- sum(top_10$counts)
  top_10$counts <- top_10$counts/total
  top_10$sample <- sampleName

  ggplot2::ggplot(top_10, ggplot2::aes(x=sample, y=counts, fill=genus)) + ggplot2::geom_bar(stat="identity", width=0.8) +
    ggplot2::xlab("") + ggplot2::ylab("relative abundance") +
    ggplot2::scale_fill_brewer(palette="Paired") +
    ggplot2::scale_y_continuous(labels = scales::percent) + ggplot2::theme_classic() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=11), axis.title=ggplot2::element_text(size=12),
                   legend.text=ggplot2::element_text(size=12),legend.title=ggplot2::element_text(size=12))
}
