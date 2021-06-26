#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include "defines.h"
#include "declarations.h"

using std::cout;        using std::endl;
using std::string;      using std::unordered_map;
using std::unordered_multimap;
using std::vector;      using std::to_string;
using std::set;         using std::map;
using std::pair;        using std::unordered_set;
using std::make_tuple;  using std::tuple;
using std::get;         using std::max;
using std::min;         using std::make_pair;

typedef unsigned long long ull;

void print_hash_res(unordered_multimap<ull, id_pos_t> const & hashmap, ull key, string const & id, ull pos)
{
    cout << key << ":" << endl;
    auto it = hashmap.equal_range(key);
    for (auto y = it.first; y != it.second; ++y) {
        if (y->second.id != id) {
            continue;
        }
        cout << "\t" << y->second.pos - pos << endl;
    }
}

void count_pos_in_ref(unordered_multimap<ull, id_pos_t> const & ref_hashmap, map<ull, int> & pos_count, ull hash,
                      string const & id, ull pos_in_long)
{
    auto it = ref_hashmap.equal_range(hash);
    for (auto y = it.first; y != it.second; ++y) {
        if (y->second.id == id) {
            auto pos = y->second.pos - pos_in_long;
            pos_count.insert({pos, 0});
            pos_count[pos]++;

        }
    }
}

void dispatch(unordered_map<string, string> & ref_map, unordered_map<string, unordered_map<string, string> > & long_map,
              unordered_map<string, unordered_multimap<ull, ull> > & ref_hashmap, vector<string> const & id_vec,
              int hash_len)
{
    for (auto const & id: id_vec) {
        concat_longs(ref_map[id], long_map[id], ref_hashmap[id], id, ref_map[id].length());
    }
}

int main()
{
    unordered_map<string, string> ref_map;
    unordered_map<string, unordered_map<string, string> > long_map;
    unordered_map<string, unordered_multimap<ull, ull> > ref_hashmap;
    vector<string> id_vec;

    read_fasta(ref_map, long_map, id_vec);
    make_hashmap(ref_map, ref_hashmap, LEN_A);

    dispatch(ref_map, long_map, ref_hashmap, id_vec, LEN_C);

    cout << "Done" << endl;

    return 0;
}
