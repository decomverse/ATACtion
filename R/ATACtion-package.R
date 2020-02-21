#' ATACtion: A multiresolution framework for interactive visualization, clustering, and annotation of single-cell epigenomic profiles
#'
#' Extenstion of ACTIONet package for single-cell epigenomic dataset
#'
#' @section Usage: 
#' 
#' \enumerate{
#' \item ?reduce.peaks.ACTION to reduce raw count matrix in an SingleCellExperiment objects.
#' \item ?run.ACTIONet to run ACTIONet on SingleCellExperiment objects.
#' }
#' @section Useful links:
#' 
#' \enumerate{
#' \item Report bugs at \url{https://github.com/shmohammadi86/ATACtion/issues}
#' }
#' 
#'
#' @name ATACtion
#' @docType package
#' @importFrom Rcpp evalCpp
#' @exportPattern ^[[:alpha:]]+ 
#' @useDynLib ATACtion, .registration=TRUE
NULL

