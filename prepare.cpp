#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "defines.h"
#include "declarations.h"

using std::cout;        using std::endl;
using std::string;      using std::unordered_map;
using std::vector;      using std::unordered_set;
using std::unordered_multimap;
using std::ifstream;    using std::set;
using std::sort;

/*
 * long_map: {ref_id: {long_id: long}}
 * id_vec: IDs of ref, sorted
 */
void read_fasta(unordered_map<string, string> & ref_map,
                unordered_map<string, unordered_map<string, string> > & long_map, vector<string> & id_vec)
{
    ifstream ref_if, long_if;
    string line;

    id_vec.clear();

    ref_if.open(REF_FASTA);
    while (getline(ref_if, line)) {
        string title = line.substr(1);
        id_vec.emplace_back(title);
        if (!getline(ref_if, line)) {
            break;
        }
        ref_map[title] = line;
    }
    ref_if.close();

    sort(id_vec.begin(), id_vec.end());
    for (auto const & id: id_vec) {
        unordered_map<string, string> new_map;
        long_map.insert({id, new_map});
    }

    long_if.open(LONG_FASTA);
    while (getline(long_if, line)) {
        string title = line.substr(1);
        if (!getline(long_if, line)) {
            break;
        }
        long_map[id_vec[title[1] - '1']].insert({title, line});
    }
    long_if.close();
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

void make_inv(string const & src, string & res)
{
    unordered_map<char, char> comp = {{'A', 'T'}, {'G', 'C'}, {'C', 'G'}, {'T', 'A'}};
    auto len = src.length();
    res.resize(len);
    for (auto i = 0; i < len; ++i) {
        res[len - 1 - i] = comp[src[i]];
    }
}

void make_inv_hashmap(unordered_map<string, string> & src_map,
                  unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > & res_map,
                  int hash_len)
{
    unsigned long long mask = get_mask_from_hash_len(hash_len);
    cout << "make_inv_hashmap: mask=" << mask << endl;

    unordered_map<char, unsigned long long> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};

    res_map.clear();
    for (auto const & item: src_map) {
        string id = item.first;
        string base = item.second;
        auto base_len = base.length();

        if (base_len < hash_len) {
            cout <<"make_inv_hashmap: base_len < hash_len" << endl;
            return;
        }

        for (auto i = hash_len - 1; i < base_len; ++i) {
            string inv;
            make_inv(base.substr(i - hash_len + 1, hash_len), inv);
            unsigned long long hash = 0;
            for (auto ch: inv) {
                hash = ((hash << 2) | dict[ch]) & mask;
            }
            if (res_map.find(id) == res_map.end()) {
                unordered_multimap<unsigned long long, unsigned long long> new_map;
                res_map.insert({id, new_map});
            }
            res_map[id].insert({hash, i});
        }
    }

    cout << "make_inv_hashmap: done." << endl;
}
