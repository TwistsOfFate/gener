#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include <fstream>
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
using std::ofstream;

typedef unsigned long long ull;


void dispatch(unordered_map<string, string> &ref_map,
              unordered_map<string, unordered_map<string, string> > &long_map,
              unordered_map<string, unordered_multimap<ull, ull> > &ref_hashmap,
              vector<string> const &id_vec)
{
    vector<string> answers;
    for (auto const & id: id_vec) {
        vector<string> concats = concat_longs(ref_map[id], long_map[id], ref_hashmap[id], id, ref_map[id].length());
        find_answers(ref_map[id], concats, id, answers);
    }



    ofstream ofs;
    ofs.open("../logs/answers.txt");
    for (const auto& x: answers) {
        ofs << x << endl;
    }
    ofs.close();
}

int main()
{
    unordered_map<string, string> ref_map;
    unordered_map<string, unordered_map<string, string> > long_map;
    unordered_map<string, unordered_multimap<ull, ull> > ref_hashmap;
    vector<string> id_vec;

    read_fasta(ref_map, long_map, id_vec);
    make_hashmap(ref_map, ref_hashmap, SEED_LEN_A);

    dispatch(ref_map, long_map, ref_hashmap, id_vec);

    cout << "Done" << endl;

    return 0;
}
