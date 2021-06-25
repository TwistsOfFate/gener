#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include <algorithm>
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
using std::sort;        using std::ofstream;

typedef unsigned long long ull;


void concat_longs(string const & ref, unordered_map<string, string> & long_map, unordered_multimap<ull, ull> const & ref_hashmap,
                  string const & id, int hash_len, ull base_len)
{
    cout << endl << id << endl;

    vector<int> ref_coverage(ref.size(), 0);

    vector<align_res_t> aligns;

    auto count = 0;
    for (auto const & item: long_map) {
        string seq = item.second;
        cout << "Matching long to ref: " << ++count << "/" << long_map.size() << "\r";
        cout.flush();
        auto align_res = align(id, item.first, ref_hashmap, seq, hash_len, base_len);
        if (align_res.r1 > 0) {
            aligns.emplace_back(align_res);
        }

        /* DEBUG */
        for (auto i = align_res.l1; i < align_res.r1; ++i) {
            ref_coverage[i]++;
        }
    }
    cout << endl;
    /* DEBUG */
    auto cover_count = 0;
    for (auto x: ref_coverage) {
        if (x) {
            cover_count++;
        }
    }
    cout << "Coverage on ref: " << cover_count << "/" << ref.size() << endl;
    cout << "Valid matches: " << aligns.size() << "/" << long_map.size() << endl;




    sort(aligns.begin(), aligns.end());
    vector<string> concats;
    concats.emplace_back(long_map[aligns[0].id2]);
    string sv = concats.back();
    ull sv_begin = aligns[0].l1;
    string result(base_len * 2, ' ');

    auto broken_num = 0;
    for (ull i = 1; i < aligns.size() - 1;) {
        cout << "Concatenating longs: " << i << "/" << aligns.size() << "\r";
        bool broken = true;

        string idi = aligns[i].id2;
        align_res_t best_res;

        // use a suffix of sv as comparing base
        string base;
        if (sv.length() > 20000) {
            base = sv.substr(sv.length() - 20000);
        } else {
            base = sv;
        }

        // find a best match to concat
        ull j = min(aligns.size() - 1, i+20);
        ull best_j = i;
        for (; j > i; --j) {
            string idj = aligns[j].id2;

            align_res_t res = align(idi, idj, base, long_map[idj], hash_len);
            if (res.r1 && res.extend_len() > best_res.extend_len() && res.extend_len() > sv.length()) {
                best_res = res;
                best_j = j;
                broken = false;
            }
        }

        if (!broken) {
            // erase and concat
            auto r1 = best_res.r1 + (sv.length() > 20000 ? sv.length() - 20000 : 0);
            sv.erase(r1);
            sv.append(long_map[aligns[best_j].id2].substr(best_res.r2));
            i = best_j + 1;
        } else {
            for (auto k = 0; k < sv.size(); ++k) {
                if (result[sv_begin + k] == ' ') {
                    result[sv_begin + k] = sv[k];
                }
            }
            concats[concats.size() - 1] = sv;

            sv = long_map[aligns[j].id2];
            sv_begin = aligns[j].l1;
            concats.emplace_back(sv);
            i = j + 1;
            best_res.clear();
        }


        /* DEBUG */
        if (broken) {
            broken_num++;
        }
    }
    cout << endl;

    cout << "broken_num: " << broken_num << "/" << aligns.size() << endl;
    cout << "number after concats: " << concats.size() << endl;

    auto result_cover = 0;
    auto result_left = 0, result_right = 0;
    for (auto i = 0; i < result.size(); ++i) {
        if (result[i] != ' ') {
            result_cover++;
            if (!result_left) {
                result_left = i;
            }
            result_right = i;
        }
    }
    cout << "result covers: " << result_cover << ", span: " << result_left << "-" << result_right << endl;

    ofstream ofs;
    ofs.open("blanks.txt");
    for (auto i = 0; i < result_right; ++i) {
        if (result[i] == ' ') {
            ofs << i << endl;
        } else {
            ofs << endl;
        }
    }
    ofs.close();
}