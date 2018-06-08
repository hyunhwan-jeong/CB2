#include "RcppArmadillo.h"
#include "AdaptiveHash.h"
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace arma;

// [[Rcpp::export]]
Rcpp::List quant(std::string ref_path, std::vector<std::string> fastq_path) {
  gRNA_Reference ref(ref_path.c_str());
  Rcpp::DataFrame df = Rcpp::DataFrame();
  std::vector<long long> sgRNA_hash;
  std::vector<std::string> sgRNA_name;
  for(auto &s : ref.lib) {
    sgRNA_hash.push_back(s.first);
    sgRNA_name.push_back(s.second);
  }
  
  const int N = sgRNA_name.size();
  Rcpp::NumericMatrix sgRNA_count(N, fastq_path.size());

  int j = 0;
  for(auto &f : fastq_path) {
    sgRNA_MAP smap(ref);
    smap.run_MAP(f.c_str());
    for(int i = 0; i < N ; ++i) {
      sgRNA_count(i,j) = smap.cnt[sgRNA_hash[i]];
    }
    j++;
  }
  return Rcpp::List::create(Rcpp::_["sgRNA"] = sgRNA_name, 
                            Rcpp::_["count"] = sgRNA_count);
}

// [[Rcpp::export]]
Rcpp::List fit_ab(const arma::mat &xvec, const arma::mat &nvec) {
  const int R = xvec.n_rows;
  const int C = xvec.n_cols;
  const int MAX_ITER = 15;
  arma::mat pvec = xvec / nvec; 
  arma::mat wvec = nvec / arma::repelem(sum(nvec, 1), 1, C);
  arma::mat phat = arma::sum(wvec % pvec, 1);
  arma::mat vmin = phat % (1-phat) / arma::sum(nvec,1);
  
  arma::mat kvec = arma::zeros<mat>(R, 2);
  arma::mat ktmp = arma::zeros<mat>(R, 1);
  arma::mat bet = arma::zeros<mat>(R, 1);
  arma::mat alp = arma::zeros<mat>(R, 1);
  arma::mat vhat = arma::zeros<mat>(R, 1);
  
  for(int i = 0 ; i < MAX_ITER ; ++i) {
    
    phat = arma::sum(wvec % pvec, 1);
    arma::mat w2 = arma::sum(wvec % wvec, 1);
    arma::mat w2n = arma::sum(wvec % wvec / nvec, 1);
    vhat = (arma::sum(wvec % wvec % pvec % pvec, 1) - w2 % phat % phat) / (1-w2);
    arma::umat valid = arma::find(vhat > vmin);
    arma::umat in_valid = arma::find(vhat <= vmin);
    bet.elem(valid) = (vhat.elem(valid) - phat.elem(valid) % 
      (1-phat.elem(valid)) % w2.elem(valid)) / 
      (phat.elem(valid) % w2n.elem(valid) - vhat.elem(valid) / (1-phat.elem(valid)));

    alp.elem(valid) = phat.elem(valid) / (1-phat.elem(valid)) % bet.elem(valid);
    wvec.rows(valid) = arma::repelem(alp.elem(valid)+bet.elem(valid), 1, C) % nvec.rows(valid) /
    (arma::repelem(alp.elem(valid)+bet.elem(valid), 1, C) + nvec.rows(valid));
    
    wvec.rows(valid) = wvec.rows(valid) / 
      arma::repelem(arma::sum(wvec.elem(valid), 1), 1, C);
    
    alp.elem(in_valid).fill(0);
    bet.elem(in_valid).fill(0);
    vhat.elem(in_valid) = vmin.elem(in_valid);
    
    kvec.col(i&1) = alp + bet;
    if(i&1) {
      ktmp.elem(valid) = mean(kvec.elem(valid), 1);
      wvec.rows(valid) = arma::repelem(ktmp.elem(valid), 1, C) % nvec.rows(valid) /
        (arma::repelem(ktmp.elem(valid), 1, C) + nvec.rows(valid));
      wvec.rows(valid) = wvec.rows(valid) / 
        arma::repelem(arma::sum(wvec.rows(valid), 1), 1, C);
    }
  }
  return Rcpp::List::create(Rcpp::_["phat"] = phat,
                            Rcpp::_["vhat"] = vhat);
}