#ifndef __ADAPTIVE_HASH__
#define __ADAPTIVE_HASH__
#include <Rcpp.h>
#include <unordered_map>
#include <map>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <string> 
#include <cctype> 
#include <cstring>
#include <algorithm>
#include <cassert>
#include <ctime>
using namespace std;

struct gRNA_Reference {
  int lib_seq_len;
  unordered_map<long long, string> lib;
  vector<string> seq;
  
  gRNA_Reference(const char *f_lib) {
    ifstream inp(f_lib);
    string name, se;
    lib_seq_len = 20;
    while(inp>>name) {
      if(name[0]=='>') name = name.substr(1);
      inp >> se;
      long long num = 0;
      for(auto &x : se) {
        num *= 4;
        num += (toupper(x)>>1)&3;
      }
      lib_seq_len = se.size();
      lib[num] = name;
      seq.push_back(se);
    }
    inp.close();
    Rcpp::Rcerr << "Detects the length of guided RNA is " << lib_seq_len << endl;
    Rcpp::Rcerr << lib.size() << " gRNAs were found from the gRNA_Reference library." << endl;
  } 
  
};

struct sgRNA_MAP {
  unordered_map<long long, int> cnt;
  unordered_map<long long, int> cnt_rc;
  
  unordered_map<int, int> pos;
  unordered_map<int, int> pos_rc;
  int num_proc_line;
  int num_hits;
  int num_hits_rc;
  long long tot_reads_len;
  bool is_rc;
  gRNA_Reference &ref;
  void run_MAP(const char *f_seq, bool need_subsample, double subsample_ratio) {
    Rcpp::Rcerr << "Reading " << f_seq << endl;
    string line;
    int num_line = 0;
    is_rc = false;
    num_proc_line = 0;
    num_hits = 0;
    num_hits_rc = 0;
    tot_reads_len = 0;
    long long mod = 1LL<<(2*ref.lib_seq_len);
    const int rc[] = {2, 3, 0, 1};
    string msg;	

    ifstream inp(f_seq);
    
    if(need_subsample) {
      REprintf("Sub-samping has been enabled. Only %.2f%% of the read will be counted.\n", subsample_ratio * 100.0);
    }
    while(getline(inp, line)) {
      if(num_line++%4!=1) continue;
      
      if(need_subsample && R::runif(0, 1) > subsample_ratio) continue;
      
      tot_reads_len += line.size();
      if(++num_proc_line%int(1e6)==0) {
        Rcpp::Rcerr << "Processing " << num_proc_line << "lines..." << endl;
        Rcpp::Rcerr << "Current Mappability: " << 100.0*max(num_hits,num_hits_rc)/(num_proc_line-1) << "%" << endl;
      }
      int cur_len = 0;
      long long num = 0;
      int i = 0;
      bool is_found = false;
      for(auto &c: line) {
        i++; 
        c = toupper(c);
        if(c=='N') {
          cur_len = 0;
          num = 0;
          continue;
        }
        num *= 4;
        num += (c>>1)&3;
        num %= mod;
        if(++cur_len==ref.lib_seq_len) {
          if(ref.lib.count(num)>0) {
            pos[i-ref.lib_seq_len]++; 
            ++num_hits;
            cnt[num]++;
            is_found = true;  
            break;
          }
          --cur_len;
        }
      }
      is_found = false;
      cur_len = 0;
      num = 0;
      reverse(line.begin(), line.end());
      i = 0;
      for(auto &c: line) {
        c = toupper(c);
        if(c=='N') {
          cur_len = 0;
          num = 0;
          continue;
        }
        num *= 4;
        num += rc[(int(c)>>1)&3];
        num %= mod;
        if(++cur_len==ref.lib_seq_len) {
          if(ref.lib.count(num)>0) {
            pos_rc[i]++;
            ++num_hits_rc;
            cnt_rc[num]++;
            is_found = true;
            break;
          }
          --cur_len;
        }
        i++;
      }
      reverse(line.begin(), line.end());
    }
    Rcpp::Rcerr << "Total " << num_proc_line << " were proceed!" << endl;
    Rcpp::Rcerr << "Final Mappability : " << 100.0*max(num_hits,num_hits_rc)/num_proc_line << "%" << endl;
    
    inp.close();
    
    if(num_hits < num_hits_rc) {
      is_rc = true;
    }
  }
  sgRNA_MAP(gRNA_Reference &r) : ref(r) {}
};

#endif