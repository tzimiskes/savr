
#include <savvy/sav_reader.hpp>
#include <Rcpp.h>

#include <vector>
#include <cstdint>
#include <string>
#include <list>

//' Reads headers and sample ID list from SAV file.
//' @param path SAV file path.
//' @return A list of headers and sample IDs.
//' @export
// [[Rcpp::export]]
bool read_sav_header(const std::string& file_path)
{
  return false;
}


//' Reads the specified region from a SAV file.
//' @param path SAV file path.
//' @param chrom Chromosome to query.
//' @param beg Start position.
//' @param end End position.
//' @param fmt Whether to read data as genotypes, allele counts, haplotype dosages or dosages (GT, AC, HDS, DS, Default: GT).
//' @return A data frame of site info and a matrix of genotype data.
//' @export
// [[Rcpp::export]]
Rcpp::List read_sav_region(const std::string& file_path, const std::string& chrom, std::int32_t beg, std::int32_t end, const std::string& fmt_str = "GT")
{ 
  savvy::fmt fmt = savvy::fmt::allele;
  if (fmt_str == "AC")
    fmt = savvy::fmt::genotype;
  else if (fmt_str == "HDS")
    fmt = savvy::fmt::haplotype_dosage;
  else if (fmt_str == "DS")
    fmt = savvy::fmt::dosage;
    
  std::list<savvy::variant<savvy::compressed_vector<float>>> tmp_vars(1);
  savvy::sav::indexed_reader file(file_path, {chrom, static_cast<std::uint64_t>(beg), static_cast<std::uint64_t>(end)}, fmt);

  while (file >> tmp_vars.back()) 
    tmp_vars.emplace_back();

  tmp_vars.pop_back();
  std::size_t nrows = tmp_vars.size();
  std::size_t ncols = file.samples().size() * file.ploidy();

  Rcpp::List result;
  
  std::vector<std::string> chromosomes(nrows);
  std::vector<std::int32_t> positions(nrows);
  std::vector<std::string> ref_alleles(nrows);
  std::vector<std::string> alt_alleles(nrows);

  Rcpp::NumericMatrix geno_data(nrows, ncols);

  std::size_t i = 0;
  for (auto it = tmp_vars.begin(); it != tmp_vars.end(); ++it,++i)
  {
    chromosomes[i] = it->chromosome();
    positions[i] = it->position();
    ref_alleles[i] = it->ref();
    alt_alleles[i] = it->alt();
    
    auto jt_end = it->data().end();
    for (auto jt = it->data().begin(); jt != jt_end; ++jt)
    {
      geno_data(i, jt.offset()) = *jt;
    }
  }
  
  result["data"] = std::move(geno_data);
  result["variants"] = Rcpp::DataFrame::create(
    Rcpp::Named("chromosome") = std::move(chromosomes),
    Rcpp::Named("position") = std::move(positions),
    Rcpp::Named("ref") = std::move(ref_alleles),
    Rcpp::Named("alt") = std::move(alt_alleles));

  return result;
}
