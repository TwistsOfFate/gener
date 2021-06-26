#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include <algorithm>
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
using std::sort;

typedef unsigned long long ull;


vector<align_res_t> align_longs_to_ref(string const & ref, unordered_map<string, string> & long_map,
                                       unordered_multimap<ull, ull> const & ref_hashmap, string const & id, ull base_len)
{
    vector<int> ref_coverage(ref.size(), 0);

    vector<align_res_t> aligns;

    auto count = 0;
    double score_sum = 0.0;
    ull query_len_sum = 0;
    for (auto const & item: long_map) {
        string seq = item.second;
        cout << "Matching long to ref: " << ++count << "/" << long_map.size() << "\r";
        cout.flush();
        auto align_res = align(id, item.first, ref_hashmap, seq, SEED_LEN_A, base_len);
        if (align_res.r1 > 0) {
            aligns.emplace_back(align_res);
        }

        /* DEBUG */
        for (auto i = align_res.l1; i < align_res.r1; ++i) {
            ref_coverage[i]++;
        }
        score_sum += align_res.score;
        query_len_sum += seq.length();
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
    cout << "Average matching accuracy: " << score_sum / query_len_sum << endl;
    cout << "Average matching score: " << score_sum / aligns.size() << endl;

    return aligns;
}


vector<string> concat_longs(string const & ref, unordered_map<string, string> & long_map,
                                 unordered_multimap<ull, ull> const & ref_hashmap, string const & id, ull base_len)
{
    cout << endl << id << endl;

    vector<align_res_t> aligns = align_longs_to_ref(ref, long_map, ref_hashmap, id, base_len);
    sort(aligns.begin(), aligns.end());

     for (auto i = 0; i < aligns.size(); ++i) {
//         cout << "align l1: " << aligns[i].l1 << endl;
     }



    vector<string> concats;
    vector<align_res_t> concat_begin;
    string sv = long_map[aligns[0].id2];
    concats.emplace_back(sv);
    concat_begin.emplace_back(aligns[0]);

    for (ull i = 1; i < aligns.size() - 1;) {
        cout << "Concatenating longs: " << i << "/" << aligns.size() << "\r";
        cout.flush();
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
        unordered_multimap<ull, ull> hashmap;
        make_hashmap(base, hashmap, SEED_LEN_B);

        // find a best match to concat
        ull j = min(aligns.size() - 1ull, i+30);
        ull best_j = i;
        for (; j > i; --j) {
            string idj = aligns[j].id2;

            align_res_t res = align(idi, idj, hashmap, long_map[idj], SEED_LEN_B, base.length());
            if (res.r1 && res.extend_len() > best_res.extend_len() && res.extend_len() > base.length()) {
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
            // save and start a new sequence
            concats[concats.size() - 1] = sv;

            sv = long_map[aligns[j].id2];
            concats.emplace_back(sv);
            concat_begin.emplace_back(aligns[j]);
            i = j + 1;
            best_res.clear();
        }

    }
    cout << endl;

    vector<string> valid_concats;
    ull not_covered = 0;
    ull valid_length = 0;

    for (auto i = 0; i < concat_begin.size(); ++i) {
        auto l = concat_begin[i].l1, r = concat_begin[i].l1 + concats[i].length();
//        if (concats[i].length() < 10000) {
//            continue;
//        }
        valid_concats.emplace_back(concats[i]);
        if (i > 0 && l > concat_begin[i-1].l1 + concats[i-1].length()) {
            not_covered += l - (concat_begin[i-1].l1 + concats[i-1].length());
        }
        valid_length = max(valid_length, concat_begin[i].l1 + concats[i].length());
    }


    /* DEBUG */
    cout << "valid concats: " << valid_concats.size() << "/" << concats.size() << endl;
    cout << "covered: " << valid_length - not_covered << "/" << valid_length << endl;

    return valid_concats;
}