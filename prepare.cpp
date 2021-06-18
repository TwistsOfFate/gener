#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <unordered_set>
#include "defines.h"

using std::cout;        using std::endl;
using std::string;      using std::unordered_map;
using std::vector;      using std::unordered_set;
using std::unordered_multimap;
using std::ifstream;    using std::set;


void read_fasta(unordered_map<string, string> & ref_map,
                unordered_map<string, unordered_set<string> > & long_map, vector<string> & id_vec)
{
    ifstream ref_if, long_if;
    string line;

    set<string> id_set;

    ref_if.open(REF_FASTA);
    while (getline(ref_if, line)) {
        string title = line.substr(1);
        id_set.insert(title);
        if (!getline(ref_if, line)) {
            break;
        }
        ref_map[title] = line;
    }
    ref_if.close();

    for (auto const & id: id_set) {
        id_vec.push_back(id);
        unordered_set<string> new_set;
        long_map.insert({id, new_set});
    }

    long_if.open(LONG_FASTA);
    while (getline(long_if, line)) {
        string title = line.substr(1);
        if (!getline(long_if, line)) {
            break;
        }
        long_map[id_vec[title[1] - '1']].insert(line);
    }
    long_if.close();
}

unsigned long long get_mask_from_hash_len(int hash_len)
{
    unsigned long long mask = 0;
    if (hash_len > 32) {
        mask = ~0;
    } else {
        for (auto i = 0; i < hash_len * 2; ++i) {
            mask = (mask << 1) + 1;
        }
    }
    return mask;
}

void make_hashmap(unordered_map<string, string> & src_map,
                  unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len)
{
    unsigned long long mask = get_mask_from_hash_len(hash_len);
    cout << "make_hashmap: mask=" << mask << endl;

    unordered_map<char, unsigned long long> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};

    res_map.clear();
    for (auto const & item: src_map) {
        string id = item.first;
        string base = item.second;
        auto base_len = base.length();

        if (base_len < hash_len) {
            cout <<"make_hashmap: base_len < hash_len" << endl;
            return;
        }

        unsigned long long hash = 0;
        auto i = 0;
        for (; i < hash_len; ++i) {
            hash = ((hash << 2) | dict[base[i]]) & mask;
        }
        if (res_map.find(id) == res_map.end()) {
            unordered_multimap<unsigned long long, unsigned long long> new_map;
            res_map.insert({id, new_map});
        }
        res_map[id].insert({hash, i - 1});
        for (; i < base_len; ++i) {
            hash = ((hash << 2) | dict[base[i]]) & mask;
            res_map[id].insert({hash, i});
        }
    }

    cout << "make_hashmap: done." << endl;
}