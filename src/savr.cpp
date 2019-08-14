
#include <savvy/sav_reader.hpp>
#include <Rcpp.h>

#include <vector>
#include <cstdint>
#include <string>
#include <list>

//' Get statistics about SAV file.
//' @param path SAV file path.
//' @return A data frame of statistics about file.
//' @export
// [[Rcpp::export]]
Rcpp::List stat_sav_file(const std::string& file_path)
{
  savvy::s1r::reader index_file(file_path + ".s1r");

  if (!index_file.good())
  {
    Rcpp::stop("Could not open index file (" + file_path + ")");
  }

  int nrows = index_file.tree_names().size();
  std::vector<std::string> chromosomes(nrows);
  std::vector<std::int32_t> variant_counts(nrows);
  std::vector<std::int32_t> min_positions(nrows);
  std::vector<std::int32_t> max_positions(nrows);

  std::vector<std::string> column_names(4);
  column_names[0] = "chromosome";
  column_names[1] = "variant_count";
  column_names[2] = "min_position";
  column_names[3] = "max_position";

  std::size_t i = 0;
  for (auto it = index_file.trees_begin(); it != index_file.trees_end(); ++it,++i)
  {
    // chromosome
    chromosomes[i] = it->name();

    // marker count
    std::uint64_t cnt = 0;
    auto q = it->create_query(0, std::numeric_limits<std::uint64_t>::max());
    for (auto e = q.begin(); e != q.end(); ++e)
    {
      cnt += std::uint32_t(0x000000000000FFFF & e->value()) + 1;
    }
    variant_counts[i] = cnt;

    std::uint32_t min;
    std::uint32_t max;
    std::tie(min, max) = it->range();

    // min pos
    min_positions[i] = min;

    // max pos
    max_positions[i] = max;
  }

  Rcpp::List ret(4);
  ret[0] = std::move(chromosomes);
  ret[1] = std::move(variant_counts);
  ret[2] = std::move(min_positions);
  ret[3] = std::move(max_positions);

  ret.attr("names") = column_names;
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -nrows);

  return ret;
}


//' Reads headers and sample ID list from SAV file.
//' @param path SAV file path.
//' @return A list of headers and sample IDs.
//' @export
// [[Rcpp::export]]
Rcpp::List read_sav_header(const std::string& file_path)
{
  savvy::sav::reader r(file_path, savvy::fmt::gt);
  if (!r.good())
  {
    Rcpp::stop("Could not open SAV file (" + file_path + ")");
  }

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
//' @param samples Which samples to extract.
//' @param chrom Chromosome to query.
//' @param beg Start position.
//' @param end End position.
//' @param tranpose Whether or not to transpose the genotype data.
//' @param fmt Whether to read data as genotypes, allele counts, haplotype dosages, dosages or genotype probabilities (GT, AC, HDS, DS, GP, Default: GT).
//' @return A data frame of site info and a matrix of genotype data.
//' @export
// [[Rcpp::export]]
Rcpp::List read_sav_region(const std::string& file_path, std::vector<std::string>& samples, const std::string& chrom, std::int32_t beg, std::int32_t end, bool transpose = false, const std::string& fmt_str = "GT")
{
  savvy::fmt fmt = savvy::fmt::gt;
  if (fmt_str == "AC")
    fmt = savvy::fmt::ac;
  else if (fmt_str == "HDS")
    fmt = savvy::fmt::hds;
  else if (fmt_str == "DS")
    fmt = savvy::fmt::ds;
  else if (fmt_str == "GP")
    fmt = savvy::fmt::gp;

  std::list<savvy::variant<savvy::compressed_vector<float>>> tmp_vars(1);
  savvy::sav::indexed_reader file(file_path, {chrom, static_cast<std::uint64_t>(beg), static_cast<std::uint64_t>(end)}, fmt);
  if (!file.good())
  {
    Rcpp::stop("Could not open indexed SAV file (" + file_path + ")");
  }

  // extract only the samples the user wants
  auto sav_samples = file.subset_samples({samples.begin(), samples.end()});

  // check to see the both vectors contain the same elements
  if (sav_samples.size() != samples.size())
    Rcpp::stop("samples contains IDs that are not in the sav file.");

  // sav file probably is not storing the samples in the same order as the
  // input, so generate an order vector that will be used to fill in the
  // genotype matrix acoording to the user's order.
  std::vector<int> order(samples.size());
  for (int j = 0; j < sav_samples.size(); ++j) {
      for (int i = 0; i < samples.size(); ++i)
        if (sav_samples[j] == samples[i])
            order[j] = i;
  }

  while (file >> tmp_vars.back())
    tmp_vars.emplace_back();

  if (file.bad())
  {
    Rcpp::stop("I/O error (" + file_path + ")");
  }

  tmp_vars.pop_back();
  std::size_t nrows = tmp_vars.size();
  std::size_t ncols;
  if (fmt == savvy::fmt::gt || fmt == savvy::fmt::hds)
    ncols = samples.size() * (file.ploidy() ? file.ploidy() : 2);
  else if (fmt == savvy::fmt::gp)
    ncols = samples.size() * ((file.ploidy() ? file.ploidy() : 2) + 1);
  else
    ncols = samples.size();

  std::vector<std::string> chromosomes(nrows);
  std::vector<std::int32_t> positions(nrows);
  std::vector<std::string> ref_alleles(nrows);
  std::vector<std::string> alt_alleles(nrows);

  Rcpp::NumericMatrix geno_data;
  if (transpose) {
      Rcpp::Rcout << "r: " << ncols << " c: " << nrows << "\n";
      geno_data = Rcpp::NumericMatrix(ncols, nrows);
  }
  else {
      Rcpp::Rcout << "r: " << nrows << " c: " << ncols << "\n";
      geno_data = Rcpp::NumericMatrix(nrows, ncols);
  }

  std::vector<std::vector<std::string>> info_columns(file.info_fields().size(), std::vector<std::string>(nrows));

  std::size_t i = 0;
  for (auto it = tmp_vars.begin(); it != tmp_vars.end(); ++it,++i)
  {
    chromosomes[i] = it->chromosome();
    positions[i] = it->position();
    ref_alleles[i] = it->ref();
    alt_alleles[i] = it->alt();

    std::size_t j = 0;
    for (auto jt = file.info_fields().begin(); jt != file.info_fields().end(); ++jt,++j)
    {
      info_columns[j][i] = it->prop(*jt);
    }

    if (ncols != it->data().size())
    {
      Rcpp::stop("Logic Error: variant vector size (" + std::to_string(it->data().size()) + ") does not match matrix row size (" + std::to_string(ncols) + ")");
    }

    j = 0;
    auto jt_end = it->data().end();
    for (auto jt = it->data().begin(); jt != jt_end; ++jt)
    {
      if (transpose)
        geno_data(order[jt.offset()], i) = *jt;
      else
        geno_data(i, order[jt.offset()]) = *jt;
    }
  }

  Rcpp::List result;

  std::vector<std::string> column_names;
  column_names.reserve(4 + file.info_fields().size());
  column_names = {"chrom","pos", "ref", "alt"};
  for (auto it = file.info_fields().begin(); it != file.info_fields().end(); ++it)
    column_names.push_back(*it);

  Rcpp::List columns(4 + file.info_fields().size());

  // Hack thanks to Matthew Flickinger
  columns.attr("names") = column_names;
  columns.attr("class") = "data.frame";
  columns.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -nrows);

  columns[0] = std::move(chromosomes);
  columns[1] = std::move(positions);
  columns[2] = std::move(ref_alleles);
  columns[3] = std::move(alt_alleles);

  for (std::size_t i = 0; i < file.info_fields().size(); ++i)
  {
    columns[i + 4] = std::move(info_columns[i]);
  }

  result["variants"] = columns;

  Rcpp::CharacterVector sample_ids = Rcpp::wrap(samples);
  if (transpose)
    geno_data.attr("dimnames") = Rcpp::List::create(sample_ids, columns["ID"]);
  else
    geno_data.attr("dimnames") = Rcpp::List::create(columns["ID"], sample_ids);
  result["data"] = std::move(geno_data);

  return result;
}
