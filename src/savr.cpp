
#include <savvy/sav_reader.hpp>
#include <Rcpp.h>

#include <vector>
#include <cstdint>
#include <string>
#include <list>

//Rcpp::DataFrame data_frame_from_list( Rcpp::List obj ){
//             bool use_default_strings_as_factors = true ;
//             bool strings_as_factors = true ;
//             int strings_as_factors_index = -1 ;
//             R_xlen_t n = obj.size() ;
//             CharacterVector names = obj.attr( "names" ) ;
//             if( !names.isNULL() ){
//                 for( int i=0; i<n; i++){
//                     if( names[i] == "stringsAsFactors" ){
//                         strings_as_factors_index = i ;
//                         use_default_strings_as_factors = false ;
//                         if( !as<bool>(obj[i]) ) strings_as_factors = false ;
//                         break ;
//                     }
//                 }
//             }
//             if( use_default_strings_as_factors )
//                 return DataFrame_Impl(obj) ;
//             SEXP as_df_symb = Rf_install("as.data.frame");
//             SEXP strings_as_factors_symb = Rf_install("stringsAsFactors");
// 
//             obj.erase(strings_as_factors_index) ;
//             names.erase(strings_as_factors_index) ;
//             obj.attr( "names") = names ;
//             Shield<SEXP> call( Rf_lang3(as_df_symb, obj, wrap( strings_as_factors ) ) ) ;
//             SET_TAG( CDDR(call),  strings_as_factors_symb ) ;
//             Shield<SEXP> res( Rcpp_eval( call ) ) ;
//             DataFrame_Impl out( res ) ;
//             return out ;
// 
//         }

//' Reads headers and sample ID list from SAV file.
//' @param path SAV file path.
//' @return A list of headers and sample IDs.
//' @export
// [[Rcpp::export]]
Rcpp::List read_sav_header(const std::string& file_path)
{
  savvy::sav::reader r(file_path, savvy::fmt::genotype);
  
  std::vector<std::string> header_names(r.headers().size());
  std::vector<std::string> header_values(r.headers().size());
  
  std::size_t i = 0;
  for (auto it = r.headers().begin(); it != r.headers().end(); ++it,++i)
  {
    header_names[i] = it->first;
    header_values[i] = it->second;
  }

  Rcpp::List ret;
  ret["headers"] = Rcpp::DataFrame::create( 
    Rcpp::Named("name") = std::move(header_names),
    Rcpp::Named("value") = std::move(header_values),
    Rcpp::Named("stringsAsFactors") = false);

  ret["sample_ids"] = r.samples();
  return ret;
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
  
  Rcpp::List result;
  
  Rcpp::List columns(4);
  std::vector<std::string> column_names = {"chrom", "pos", "ref", "alt"};

  // Hack thanks to Matthew Flickinger
  columns.attr("names") = column_names;
  columns.attr("class") = "data.frame";
  columns.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -nrows);

  columns[0] = std::move(chromosomes);
  columns[1] = std::move(positions);
  columns[2] = std::move(ref_alleles);
  columns[3] = std::move(alt_alleles);

  result["variants"] = columns;
  
  result["data"] = std::move(geno_data);

  return result;
}
