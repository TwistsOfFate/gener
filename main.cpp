#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
#include <cmath>

#include "defines.h"
#include "declarations.h"

//#include <cstdlib>
//#include <ctime>

using std::cout;        using std::endl;
using std::string;      using std::unordered_map;
using std::unordered_multimap;
using std::vector;      using std::to_string;
using std::set;         using std::map;
using std::pair;        using std::unordered_set;
using std::make_tuple;  using std::tuple;
using std::get;         using std::max;
using std::min;         using std::make_pair;

void print_hash_res(unordered_multimap<unsigned long long, id_pos_t> const & hashmap, unsigned long long key,
                    string const & id, unsigned long long pos)
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

void count_pos_in_ref(unordered_multimap<unsigned long long, id_pos_t> const & ref_hashmap,
                      map<unsigned long long, int> & pos_count, unsigned long long hash, string const & id,
                      unsigned long long pos_in_long)
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

// XuYipei code ------------------------------------

//inline int randint_generate(int a, int b) {
//    srand((unsigned)time(NULL));
//    return rand() % (b - a + 1) + a;
//}
//
//inline vector<unsigned long long> set_join(vector<unsigned long long> &a, vector<unsigned long long> &b) {
//    unordered_map<unsigned long long, unsigned long long> key_count;
//    vector<unsigned long long> result;
//    for (unsigned long long & i : a) {
//        key_count.insert({i, 1});
//    }
//    for (unsigned long long & i : b) {
//        key_count.insert({i, 0});
//        key_count[i]++;
//    }
//    for (auto x: key_count) {
//        if (x.second == 2) {
//            result.push_back(x.first);
//        }
//    }
//    return result;
//}
//
//void rand_probe(unordered_multimap<unsigned long long, id_pos_t> const & ref_hashmap, vector<unsigned long long> hash_vec,
//                string const & id, unsigned long long seq_len, int hash_len)
//{
//    int looks = 4;
//    vector<unsigned long long> best_result;
//    for (auto j = 0; j < 50; j++) {
//        vector<unsigned long long> poss[50];
//        vector<unsigned long long> result;
//        for (auto i = 0; i < looks; ++i) {
//            auto r = randint_generate(0, seq_len - hash_len - 5);
//            auto it = ref_hashmap.equal_range(hash_vec[r]);
//            for (auto y = it.first; y != it.second; ++y) {
//                if (y->second.id != id) {
//                    continue;
//                }
//                poss[i].push_back(y->second.pos - r);
//            }
//            if (i == 0) {
//                result = poss[i];
//            } else {
//                result = set_join(result, poss[i]);
//            }
//        }
//        if (best_result.empty() || best_result.size() > result.size()) {
//            best_result = result;
//        }
//        if (best_result.size() == 1) {
//            break;
//        }
//    }
//
//    for (auto x: best_result) {
//        cout << x << " ";
//    }
//    cout << endl;
//}

// -------------------------------------------------

void insert_anchors(unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                    map<unsigned long long, pair<unsigned long long, unsigned long long> > & anchor_map,
                    unsigned long long hash, unsigned long long y, unsigned long long w)
{
    auto range = ref_hashmap.equal_range(hash);
    for (auto it = range.first; it != range.second; ++it) {
        anchor_map.insert({it->second, {y, w}});
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

void perform_chaining(map<unsigned long long, pair<unsigned long long, unsigned long long> > const & anchor_map,
                      vector<double> & chain_scores, vector<int> & ref_coverage, int seed_len, unsigned long long seq_len)
{
    if (anchor_map.empty()) {
        return;
    }

    vector<tuple<unsigned long long, unsigned long long, unsigned long long> > anchors(anchor_map.size());
    vector<unsigned long long> pred;
    pair<unsigned long long, double> best_score;
    for (auto const & anchor: anchor_map) {
        anchors.emplace_back(anchor.first, anchor.second.first, anchor.second.second);
    }

    // Initialize first score
    chain_scores.push_back(seed_len);
    pred.push_back(0);
    best_score = {0, seed_len};

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
        if (best_score.second < chain_scores[i]) {
            best_score = {i, chain_scores[i]};
        }
    }

    unsigned long long last_occur = get<0>(anchors[best_score.first]);
    unsigned long long first_occur = 0;

    for (auto i = best_score.first; i > 0; i = pred[i]) {
//        cout << get<0>(anchors[i]) << "\t" << get<1>(anchors[i]) << endl;
        if (i == pred[i]) {
            first_occur = get<0>(anchors[i]) - seed_len + 1;
            break;
        }
    }

    if (best_score.second < seq_len / 5.0 || last_occur - first_occur > seq_len * 2) {
        return;
    }

    for (auto i = first_occur; i <= last_occur; ++i) {
        ref_coverage[i]++;
    }

//    static auto total_occurs = 0;
//    static auto total_len = 0;
//    total_occurs += last_occur - first_occur;
//    total_len += seq_len;
//    cout << total_occurs << "/" << total_len << "="
//         << 1.0 * total_occurs / total_len << "\r";
//    cout.flush();
}

void query_vs_reference(string const & query, unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                        int hash_len, vector<int> & ref_coverage)
{
    unordered_map<char, unsigned long long> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};
    map<unsigned long long, pair<unsigned long long, unsigned long long> > anchor_map;
    auto query_len = query.length();
    unsigned long long mask = get_mask_from_hash_len(hash_len);

    unsigned long long hash = 0;
    auto i = 0;
    for (; i < hash_len; ++i) {
        hash = ((hash << 2) | dict[query[i]]) & mask;
    }
    insert_anchors(ref_hashmap, anchor_map, hash, i - 1, hash_len);
    for (; i < query_len; ++i) {
        hash = ((hash << 2) | dict[query[i]]) & mask;
        insert_anchors(ref_hashmap, anchor_map, hash, i, hash_len);
    }
//    cout << "query_vs_reference: anchor_map size = " << anchor_map.size() << endl;

    vector<double> chain_scores;
    perform_chaining(anchor_map, chain_scores, ref_coverage, hash_len, query_len);
}

void process_longs(string const & ref, unordered_set<string> const & long_set,
                   unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                   string const & id, int hash_len, vector<pair<unsigned long long, unsigned long long> > & suspects)
{
    cout << endl << id << endl;

    vector<int> ref_coverage(ref.size(), 0);

    auto count = 0;
    for (auto const & seq: long_set) {
//        cout << ++count << "/" << long_set.size() << "\n";
        cout.flush();
        query_vs_reference(seq, ref_hashmap, hash_len, ref_coverage);
    }

    auto covers = 0;
    for (auto x: ref_coverage) {
        if (x) {
            covers++;
        }
    }
    cout << "ref coverage: " << covers << "/" << ref_coverage.size() << endl;

    auto error_series = 0;
    suspects.clear();
    for (auto i = 0; i < ref_coverage.size(); ++i) {
        if (!ref_coverage[i]) {
            error_series++;
        } else {
            if (error_series > LEN_B) {
                suspects.emplace_back(make_pair(i - error_series, i - 1));
            }
            error_series = 0;
        }
    }
    cout << "Total " << suspects.size() << " intervals" << endl;
}

void check_inv(string const & ref_base, vector<pair<unsigned long long, unsigned long long> > & suspects,
               unordered_set<string> const & long_set, int hash_len)
{
    unordered_map<string, string> ref_map;
    ref_map.insert({"A", ref_base});
//    for (auto i = 0; i < suspects.size(); ++i) {
//        auto l = suspects[i].first;
//        auto r = suspects[i].second;
//        auto old_str = ref_base.substr(l, r - l + 1);
//        string new_str;
//        make_inv(old_str, new_str);
//        ref_map.insert({to_string(i), new_str});
//    }
    unordered_set<string> long_inv_set;
    for (const auto & it : long_set) {
        string new_str;
        make_inv(it, new_str);
        long_inv_set.insert(new_str);
    }

    unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > ref_hashmap;
    make_hashmap(ref_map, ref_hashmap, hash_len);

    auto ref_len = ref_base.size();
    vector<pair<unsigned long long, unsigned long long> > misses;
    process_longs(ref_base, long_inv_set, ref_hashmap["A"], "A", hash_len, misses);

    vector<int> ref_coverage(ref_len, 1);
    for (auto x: misses) {
        for (auto i = x.first; i <= x.second; ++i) {
            ref_coverage[i] = 0;
        }
    }
    for (auto x: suspects) {
        auto covers = 0;
        for (auto i = x.first; i <= x.second; ++i) {
            covers += ref_coverage[i];
        }
        cout << x.first << " " << x.second << ":\t" << covers << "/" << x.second - x.first + 1 << endl;
    }

//    for (auto i = 0; i < suspects.size(); ++i) {
//        bool found_inv = false;
//        auto ref_len = ref_map[to_string(i)].size();
//
//        vector<int> ref_coverage(ref_len, 0);
//
//        auto count = 0;
//        for (auto & query: long_set) {
//            query_vs_reference(query, ref_hashmap[to_string(i)], hash_len, ref_coverage);
//        }
//
//        auto cover_count = 0;
//        for (auto x: ref_coverage) {
//            cover_count += (x > 0);
//        }
//        cout << suspects[i].first << " " << suspects[i].second << ":\t" << cover_count << "/" << ref_len << endl;
//        if (cover_count > ref_len / 2) {
////            cout << "INV " << suspects[i].first << " " << suspects[i].second << endl;
//        }
//
//    }
}

void dispatch(unordered_map<string, string> & ref_map, unordered_map<string, unordered_set<string> > & long_map,
              unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > & ref_hashmap,
              vector<string> const & id_vec, int hash_len)
{
    for (auto const & id: id_vec) {
        vector<pair<unsigned long long, unsigned long long> > suspects;
        process_longs(ref_map[id], long_map[id], ref_hashmap[id], id, hash_len, suspects);
        check_inv(ref_map[id], suspects, long_map[id], hash_len);
    }
}

int main()
{
    unordered_map<string, string> ref_map;
    unordered_map<string, unordered_set<string> > long_map;
    unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > ref_hashmap;
    vector<string> id_vec;

    read_fasta(ref_map, long_map, id_vec);
    make_hashmap(ref_map, ref_hashmap, LEN_A);

    dispatch(ref_map, long_map, ref_hashmap, id_vec, LEN_A);


    return 0;
}
