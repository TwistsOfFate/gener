#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <bitset>
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
using std::ifstream;    using std::vector;
using std::ofstream;    using std::bitset;
using std::set;         using std::map;
using std::pair;        using std::unordered_set;
using std::make_tuple;  using std::tuple;
using std::get;         using std::max;
using std::min;

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

//void exam_hashmap(unordered_multimap<unsigned long long, id_pos_t> const & ref_hashmap)
//{
//    unordered_map<unsigned long long, int> key_count;
//
//    int multi_occur_count = 0;
//    int total_keys = 0;
//
//    for (auto const & x: ref_hashmap) {
//        if (key_count.find(x.first) == key_count.end()) {
//            total_keys++;
//        }
//        key_count.insert({x.first, 0});
//        key_count[x.first]++;
//    }
//
//
//    for (auto const & x: key_count) {
//        if (x.second > 2) {
//            multi_occur_count++;
////            cout << x.first << ":\t" << x.second << endl;
////            auto it = ref_hashmap.equal_range(x.first);
////            for (auto y = it.first; y != it.second; ++y) {
////                y->second.print(cout);
////                cout << endl;
////            }
//        }
//    }
//
//    cout << "multiple occurs: " << multi_occur_count << "/" << total_keys << "=" << double(multi_occur_count) / total_keys << endl;
//    cout << "average occurs: " << double(ref_hashmap.size()) / total_keys << endl;
//}

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

void perform_chaining(map<unsigned long long, pair<unsigned long long, unsigned long long> > const & anchor_map,
                      vector<double> & chain_scores, vector<int> & ref_coverage, int seed_len, unsigned long long seq_len)
{
    const double INF = 9999999;
    const unsigned long long G = 200;
    vector<tuple<unsigned long long, unsigned long long, unsigned long long> > anchors;
    vector<unsigned long long> pred;
    pair<unsigned long long, double> best_score;
    for (auto const & anchor: anchor_map) {
        anchors.push_back(make_tuple(anchor.first, anchor.second.first, anchor.second.second));
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
            double candidate = chain_scores[j];

            candidate += min(double(min(yi - yj, xi - xj)), double(wi));
            double beta;
            if (yj > yi || max(yi - yj, xi - xj) > G) {
                beta = INF;
            } else {
                auto l = (yi - yj) - (xi - xj);
                if (l < 0) {
                    l = -l;
                }
                if (l != 0) {
                    beta = 0.01 * seed_len * double(l) + 0.5 * log2(double(l));
                } else {
                    beta = 0;
                }
            }
            candidate -= beta;

            if (chain_scores[i] < candidate) {
                chain_scores[i] = candidate;
                pred[i] = j;
            }
        }
        if (best_score.second < chain_scores[i]) {
            best_score = {i, chain_scores[i]};
        }
    }

//    cout << "best score: " << best_score.second << endl;

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

void process_longs(string const & ref, unordered_set<string> const & long_set,
                   unordered_multimap<unsigned long long, unsigned long long> const & ref_hashmap,
                   string const & id, int hash_len)
{
    unsigned long long mask = get_mask_from_hash_len(hash_len);
    cout << id << endl;

    vector<int> ref_coverage(ref.size(), false) ;
    auto count = 0;

    for (auto const & seq: long_set) {
        auto seq_len = seq.length();
//        cout << ++count << "/" << long_set.size() << " ";
        cout.flush();
        unordered_map<char, unsigned long long> dict = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};
        map<unsigned long long, pair<unsigned long long, unsigned long long> > anchor_map;

        unsigned long long hash = 0;
        auto i = 0;
        for (; i < hash_len; ++i) {
            hash = ((hash << 2) | dict[seq[i]]) & mask;
        }
        insert_anchors(ref_hashmap, anchor_map, hash, i - 1, hash_len);
        for (; i < seq_len; ++i) {
            hash = ((hash << 2) | dict[seq[i]]) & mask;
            insert_anchors(ref_hashmap, anchor_map, hash, i, hash_len);
        }

        vector<double> chain_scores;
        perform_chaining(anchor_map, chain_scores, ref_coverage, hash_len, seq_len);
    }

    auto covers = 0;
    for (auto x: ref_coverage) {
        if (x) {
            covers++;
        }
    }

    cout << "\nref coverage: " << covers << "/" << ref_coverage.size() << endl;

    auto error_series = 0;
    auto cover_series = 0;
    auto cover_count = 0;
    auto error_count = 0;
    for (auto i = 0; i < ref_coverage.size(); ++i) {
        if (!ref_coverage[i]) {
            error_series++;
        } else {
            if (error_series > 100) {
                cout << i - error_series << "\t" << i - 1 << endl;
                ++error_count;
            }
            error_series = 0;
        }
        if (cover_count == 0 && ref_coverage[i]) {
            cover_count = ref_coverage[i];
        }
        if (ref_coverage[i] == cover_count && cover_count != 0) {
            cover_series++;
        } else if (cover_count) {
//            cout << i - cover_series << "-" << i - 1 << ": " << cover_count << endl;
            cover_count = 0;
            cover_series = 0;
        }
    }
    cout << "Total " << error_count << " intervals" << endl;
}

void dispatch(unordered_map<string, string> & ref_map, unordered_map<string, unordered_set<string> > & long_map,
              unordered_map<string, unordered_multimap<unsigned long long, unsigned long long> > & ref_hashmap,
              vector<string> const & id_vec, int hash_len)
{
    for (auto const & id: id_vec) {
        process_longs(ref_map[id], long_map[id], ref_hashmap[id], id, hash_len);
    }
}

void exam_data(unordered_map<string, string> & ref_map, unordered_map<string, string> & long_map)
{
    unordered_map<string, int> long_visited;

    for (const auto& item: long_map) {
        long_visited.insert({item.first, 0});
    }

    auto long_size = long_map.size();
    auto long_matched = 0;

    for (const auto& item: ref_map) {
        string base = item.second;
        auto base_len = base.length();
        auto last_perc = base_len;
        last_perc = 0;
        for (auto i = 0; i + LEN_A - 1 < base_len; i += LEN_A) {
            string base_excerpt = base.substr(i, LEN_A);
            bool printed = false;

            auto perc = 100 * i / base_len;
            if (perc >= last_perc + 1) {
                cout << "base: " << perc << "%" << endl;
                last_perc = perc;
            }

            for (const auto& x: long_map) {
                if (x.second.find(base_excerpt) != string::npos) {
                    if (!printed) {
                        cout << base_excerpt << " matched: " << endl;
                        printed = true;
                    }
                    cout << "\tin " << x.first << endl;
                    auto old_size = long_visited[x.first];
                    long_visited[x.first]++;
                    if (old_size == 0) {
                        cout << ++long_matched << "/" << long_size << endl;
                    }
                }
            }
        }
    }

//    for (const auto& item: long_map) {
//        long_visited.insert({item.first, 0});
//    }
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

//    exam_data(ref_map, long_map);
//    work(ref_map, sv_map);


    return 0;
}
