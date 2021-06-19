#ifndef GENER_DECLARATIONS_H
#define GENER_DECLARATIONS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

// prepare.cpp
void read_fasta(std::unordered_map<std::string, std::string> & ref_map,
                std::unordered_map<std::string, std::unordered_set<std::string> > & long_map,
                std::vector<std::string> & id_vec);

void make_hashmap(std::unordered_map<std::string, std::string> & src_map,
                  std::unordered_map<std::string, std::unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len);

void make_inv(std::string const &, std::string &);

void make_inv_hashmap(std::unordered_map<std::string, std::string> & src_map,
                  std::unordered_map<std::string, std::unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len);

unsigned long long get_mask_from_hash_len(int);

#endif //GENER_DECLARATIONS_H
