#ifndef __ADAPTIVE_HASH__
#define __ADAPTIVE_HASH__
#include <Rcpp.h>
#include <unordered_map>
#include <map>
#include <set>
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
  unordered_map<string, string> id2seq;
  vector<string> seq;
  
  gRNA_Reference(const char *f_lib, bool verbose = false) {
    ifstream inp(f_lib);
    
    string name, se;
    
    unordered_map<string, int> sgrna_count;
    while(inp>>name) {
      if(name[0]=='>') name = name.substr(1);
      inp >> se;
      sgrna_count[se]++;
    }
    
    int tot_dups = 0;
    for(auto &it: sgrna_count) {
      if(it.second > 1) {
        ++tot_dups;
      }
    } 
   
    if(verbose) {
      Rcpp::Rcerr << tot_dups << " sgRNA sequences were repetitive and will be discarded." << endl;   
    }   
    
    lib_seq_len = 20;
    
    inp.close();
    inp.clear();
    inp.open(f_lib);
    
    while(inp>>name) {
      if(name[0]=='>') name = name.substr(1);
      inp >> se;
      if(sgrna_count[se]!=1) continue;
      long long num = 0;
      for(auto &x : se) {
        num *= 4;
        num += (toupper(x)>>1)&3;
      }
      
      lib_seq_len = se.size();
      lib[num] = name;
      id2seq[name] = se;
      seq.push_back(se);
    }
    inp.close();
    
    if(verbose) {
      Rcpp::Rcerr << "CB2 Detects the length of guided RNA is " << lib_seq_len << endl;
      Rcpp::Rcerr << lib.size() << " gRNAs were found from the gRNA_Reference library." << endl;
    }
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
  long long mod;
  gRNA_Reference &ref;
  
  bool verbose;
  
  void search(string &line) {
    const static int rc[] = {2, 3, 0, 1};
    
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
    if(is_found) {
      return;
    }
    
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
          break;
        }
        --cur_len;
      }
      i++;
    }
    reverse(line.begin(), line.end());
  }
  
  void run_MAP(const char *f_seq) {
    
    if(verbose) {
      Rcpp::Rcerr << "Reading " << f_seq << endl;
    }
    string line;
    int num_line = 0;
    is_rc = false;
    num_proc_line = 0;
    num_hits = 0;
    num_hits_rc = 0;
    tot_reads_len = 0;
    mod = 1LL<<(2*ref.lib_seq_len);
    string msg;	

    ifstream inp(f_seq);
    
    while(getline(inp, line)) {
      if(num_line++%4!=1) continue;
      
      tot_reads_len += line.size();
      if(++num_proc_line%int(1e6)==0 && verbose) {
        Rcpp::Rcerr << "Processing " << num_proc_line << "th line..." << endl;
        Rcpp::Rcerr << "Current Mappability: " << 100.0*max(num_hits,num_hits_rc)/(num_proc_line-1) << "%" << endl;
      }
      
      search(line);
    }
    
    if(verbose) {
      Rcpp::Rcerr << "Total " << num_proc_line << " were proceed!" << endl;
      Rcpp::Rcerr << "Final Mappability: " << 100.0*max(num_hits,num_hits_rc)/num_proc_line << "%" << endl;
    }
    
    inp.close();
    
    if(num_hits < num_hits_rc) {
      is_rc = true;
    }      
  }

  sgRNA_MAP(gRNA_Reference &r, bool v) : ref(r), verbose(v) {}
};

#endif