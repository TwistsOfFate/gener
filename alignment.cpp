#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include <cmath>
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


unordered_multimap<ull, ull> make_hashmap_from_str(string const & base, int hash_len)
{
    ull mask = get_mask_from_hash_len(hash_len);
    unordered_map<char, ull> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};
    unordered_multimap<ull, ull> res_map;
    auto base_len = base.length();

    if (base_len < hash_len) {
        cout << "make_hashmap: base_len < hash_len" << endl;
        return res_map;
    }

    ull hash = 0;
    auto i = 0;
    for (; i < hash_len; ++i) {
        hash = ((hash << 2) | dict[base[i]]) & mask;
    }
    res_map.insert({hash, i - 1});
    for (; i < base_len; ++i) {
        hash = ((hash << 2) | dict[base[i]]) & mask;
        res_map.insert({hash, i});
    }

    return res_map;
}

void insert_anchor(unordered_multimap<ull, ull> const & ref_hashmap, vector<tuple<ull, ull, ull> > & anchor_vec,
                   ull hash, ull y, ull w)
{
    auto range = ref_hashmap.equal_range(hash);
    for (auto it = range.first; it != range.second; ++it) {
        anchor_vec.emplace_back(make_tuple(it->second, y, w));
    }
}

double calc_chain_score(double fj, double xi, double xj, double yi, double yj, double wi, double wj, double seed_len)
{
    const double INF = 9999999;
    const double G = LEN_B;

    fj += min(min(yi - yj, xi - xj), wi);

    double beta;
    if (yj > yi || max(yi - yj, xi - xj) > G) {
        beta = INF;
    } else {
        auto l = (yi - yj) - (xi - xj);
        if (l < 0) {
            l = -l;
        }
        if (l >= 0.1) { // l != 0
            beta = 0.01 * seed_len * l + 0.5 * log2(l);
        } else {
            beta = 0;
        }
    }
    fj -= beta;

    return fj;
}

align_res_t perform_chaining(const string & id1, const string & id2, vector<tuple<ull, ull, ull> > const & anchors,
                             int seed_len, ull base_len, ull seq_len)
{
    if (anchors.empty()) {
        return align_res_t(id1, id2, base_len, seq_len, 0, 0, 0, 0, 0.0);
    }

    vector<ull> pred;
    ull best_score_id = 0;
    vector<double> chain_scores;

    chain_scores.push_back(seed_len);
    pred.push_back(0);

    for (auto i = 1; i < anchors.size(); ++i) {
        auto xi = get<0>(anchors[i]), yi = get<1>(anchors[i]), wi = get<2>(anchors[i]);
        chain_scores.push_back(double(wi));
        pred.push_back(i);
        for (auto j = i - 1; j >= i - 50 && j >= 0; --j) {
            auto xj = get<0>(anchors[j]), yj = get<1>(anchors[j]), wj = get<2>(anchors[j]);
            double candidate = calc_chain_score(chain_scores[j], xi, xj, yi, yj, wi, wj, seed_len);
            if (chain_scores[i] < candidate) {
                chain_scores[i] = candidate;
                pred[i] = j;
            }
        }
        if (chain_scores[best_score_id] < chain_scores[i]) {
            best_score_id = i;
        }
    }

    ull r1 = get<0>(anchors[best_score_id]), r2 = get<1>(anchors[best_score_id]);
    ull l1 = 0, l2 = 0;

    for (auto i = best_score_id; i > 0; i = pred[i]) {
        if (i == pred[i]) {
            l1 = get<0>(anchors[i]) - seed_len + 1;
            l2 = get<1>(anchors[i]) - seed_len + 1;
            break;
        }
    }

    if (chain_scores[best_score_id] < seq_len / 10.0 && chain_scores[best_score_id] < LEN_B || r1 - l1 > seq_len * 2) {
        return align_res_t(id1, id2, base_len, seq_len, 0, 0, 0, 0, chain_scores[best_score_id]);
    } else {
        return align_res_t(id1, id2, base_len, seq_len, l1, r1, l2, r2, chain_scores[best_score_id]);
    }
}

align_res_t align(const string & id1, const string & id2, unordered_multimap<ull, ull> const & ref_hashmap,
                  string const & query, int hash_len, ull base_len)
{
    unordered_map<char, ull> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};
    map<ull, pair<ull, ull> > anchor_map;
    vector<tuple<ull, ull, ull> > anchor_vec;
    auto query_len = query.length();
    ull mask = get_mask_from_hash_len(hash_len);

    ull hash = 0;
    auto i = 0;
    for (; i < hash_len; ++i) {
        hash = ((hash << 2) | dict[query[i]]) & mask;
    }
    insert_anchor(ref_hashmap, anchor_vec, hash, i - 1, hash_len);
    for (; i < query_len; ++i) {
        hash = ((hash << 2) | dict[query[i]]) & mask;
        insert_anchor(ref_hashmap, anchor_vec, hash, i, hash_len);
    }
    sort(anchor_vec.begin(), anchor_vec.end());


    return perform_chaining(id1, id2, anchor_vec, hash_len, base_len, query_len);
}

align_res_t align(const string & id1, const string & id2, string const & str1, string const & str2, int hash_len)
{
    return align(id1, id2, make_hashmap_from_str(str1, hash_len), str2, hash_len, str1.length());
}
