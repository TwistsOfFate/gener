#ifndef GENER_DECLARATIONS_H
#define GENER_DECLARATIONS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

// prepare.cpp
void read_fasta(std::unordered_map<std::string, std::string> & ref_map,
                std::unordered_map<std::string, std::unordered_map<std::string, std::string> > & long_map,
                std::vector<std::string> & id_vec);

void make_hashmap(std::string & base, std::unordered_multimap<unsigned long long, unsigned long long> & res_map, int hash_len);

void make_hashmap(std::unordered_map<std::string, std::string> & src_map,
                  std::unordered_map<std::string, std::unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len);

void make_inv(std::string const &, std::string &);

void make_inv_hashmap(std::unordered_map<std::string, std::string> & src_map,
                  std::unordered_map<std::string, std::unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len);


// commons.cpp
unsigned long long get_mask_from_hash_len(int);
std::string inverse(const std::string & src);


// alignment.cpp
align_res_t align(const std::string & id1, const std::string & id2,
                  std::unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                  std::string const & query, int hash_len, unsigned long long base_len);

align_res_t align(const std::string & id1, const std::string & id2,
                  std::string const & str1, std::string const & str2, int hash_len);


// concat.cpp
std::vector<std::string> concat_longs(std::string const & ref, std::unordered_map<std::string, std::string> & long_map,
                                      std::unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                                      std::string const & id, unsigned long long base_len);


// compare.cpp
std::vector<std::string> find_answers(std::string & ref, std::vector<std::string> & concats, std::string const & id);

#endif //GENER_DECLARATIONS_H
