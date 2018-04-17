#' Reads headers and sample ID list from SAV file.
#' @param path SAV file path.
#' @return A list of headers and sample IDs.
#' @export

read_sav_header <- function(path){

}

#' Reads the specified region from a SAV file.
#' @param path SAV file path.
#' @param chrom Chromosome to query.
#' @param beg Start position (Default: 0).
#' @param end End position (Default: 2147483647).
#' @param fmt Whether to read data as genotypes, allele counts, haplotype dosages or dosages (GT, AC, HDS, DS, Default: GT).
#' @return A data frame of site info and a matrix of genotype data.
#' @export
 
read_sav_region <- function(path, chrom, beg=0, end=2147483647, fmt="GT"){
    
}
