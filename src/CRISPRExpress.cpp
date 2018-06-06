#include <Rcpp.h>
#include "AdaptiveHash.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List quant(std::string ref_path, std::vector<std::string> fastq_path) {
  gRNA_Reference ref(ref_path.c_str());
  DataFrame df = DataFrame();
  std::vector<long long> sgRNA_hash;
  std::vector<std::string> sgRNA_name;
  for(auto &s : ref.lib) {
    sgRNA_hash.push_back(s.first);
    sgRNA_name.push_back(s.second);
  }
  
  const int N = sgRNA_name.size();
  NumericMatrix sgRNA_count(N, fastq_path.size());

  int j = 0;
  for(auto &f : fastq_path) {
    sgRNA_MAP smap(ref);
    smap.run_MAP(f.c_str());
    for(int i = 0; i < N ; ++i) {
      sgRNA_count(i,j) = smap.cnt[sgRNA_hash[i]];
    }
    j++;
  }
  return List::create(_["sgRNA"] = sgRNA_name, 
                      _["count"] = sgRNA_count);
}