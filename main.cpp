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

//#include <cstdlib>
//#include <ctime>

#define REF_FASTA "../data/ref.fasta"
#define SV_FASTA "../data/sv.fasta"
#define LONG_FASTA "../data/long.fasta"
#define SVS_BED "../sv.bed"
#define IDENTICAL_LEN 50
#define MAX_SV_LEN 1005
#define MIN_SV_LEN 50
#define LEN_A 12

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

struct ans_t {
    std::string type;
    std::string id1;
    int pos1, len1;
    std::string id2;
    int pos2, len2;

    void update(const std::string & id1, const int & pos1, const int & len1)
    {
        this->id1 = id1;
        this->pos1 = pos1;
        this->len1 = len1;
    }

    void update(const std::string & id1, const int & pos1, const int & len1, const std::string & id2, const int & pos2, const int & len2)
    {
        this->id1 = id1;
        this->pos1 = pos1;
        this->len1 = len1;
        this->id2 = id2;
        this->pos2 = pos2;
        this->len2 = len2;
    }

    void print(std::ostream & os = std::cout) const
    {
        os << type << " " << id1 << " " << pos1 << " " << pos1 + len1;
        if (type == "UNK" || type == "TRA") {
            os << " " << id2 << " " << pos2 << " " << pos2 + len2;
        }
        os << std::endl;
    }
};

struct id_pos_t {
    std::string id;
    unsigned long long pos;

    id_pos_t(std::string id, unsigned long long pos): id(std::move(id)), pos(pos) {}

    void print(std::ostream & os = std::cout) const
    {
        os << id << "\t" << pos;
    }
};

inline bool base_complement(const char & x, const char & y)
{
    static const string dict = "ACTG";

    for (int i = 0; i < 4; ++i) {
        if (x == dict[i]) {
            return y == dict[(i + 2) % 4];
        }
    }

    return false;
}

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

inline unsigned long long get_mask_from_hash_len(int hash_len)
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

bool check_ins(const string & ref, const string & sv, int & i, int & j, int & pos, int & len)
{
    auto ref_len = ref.size();
    auto sv_len = sv.size();
    int delta_j = MIN_SV_LEN;

    // too far back
    if (i + IDENTICAL_LEN >= ref_len) {
        return false;
    }

    string ref_base = ref.substr(i, IDENTICAL_LEN);

    for (; j + delta_j + IDENTICAL_LEN < sv_len && delta_j <= MAX_SV_LEN; ++delta_j) {
        if (sv.compare(j + delta_j, IDENTICAL_LEN, ref_base) == 0) {
            j = j + delta_j;
            pos = i;
            len = delta_j;
            return true;
        }
    }

    return false;
}

// j, pos and len are results of check_ins()
bool check_dup(const string & ref, const string & sv, int & i, int & j, int & pos, int & len)
{
    int dup_base = j - len;
    for (; dup_base > j - MAX_SV_LEN && dup_base > 0; --dup_base) {
        int dup_len = j - dup_base;
        if (sv.compare(dup_base, dup_len, sv, dup_base - dup_len, dup_len) == 0) {
            pos = i - dup_len - (j - len - dup_base);
            len = dup_len;
            return true;
        }
    }
    return false;
}

bool check_del(const string & ref, const string & sv, int & i, int & j, int & pos, int & len)
{
    auto ref_len = ref.size();
    auto sv_len = sv.size();
    int delta_i = MIN_SV_LEN;

    // too far back
    if (j + IDENTICAL_LEN >= sv_len) {
        return false;
    }

    string sv_base = sv.substr(j, IDENTICAL_LEN);

    for (; i + delta_i + IDENTICAL_LEN < ref_len && delta_i <= MAX_SV_LEN; ++delta_i) {
        if (ref.compare(i + delta_i, IDENTICAL_LEN, sv_base) == 0) {
            pos = i;
            i = i + delta_i;
            len = delta_i;
            return true;
        }
    }

    return false;
}

bool check_inv(const string & ref, const string & sv, int & i, int & j, int & pos, int & len)
{
    int inv_len = MIN_SV_LEN;
    auto ref_len = ref.size();
    auto sv_len = sv.size();
    bool inv_exists = false;

    for (; inv_len <= MAX_SV_LEN && i + inv_len < ref_len && j + inv_len < sv_len; ++inv_len) {
        // check if ref[i:i+inv_len] and sv[j:j+inv_len] are inversions
        int delta = 0;
        for (; delta < inv_len; ++delta) {
            if (!base_complement(ref[i + delta], sv[j + inv_len - delta - 1])) {
                break;
            }
        }
        if (delta == inv_len) {
            pos = i;
            len = inv_len;
            inv_exists = true;
        }
    }

    if (inv_exists) {
        i = i + len;
        j = j + len;
    }

    return inv_exists;
}

void check_unk(const string & ref, const string & sv, int & i, int & j, int & pos, int & len, int & pos2, int & len2)
{
    auto ref_len = ref.size();
    auto sv_len = sv.size();
    int delta_i = MAX_SV_LEN;
    int delta_j = 0;

    // too far back
    if (i + delta_i + IDENTICAL_LEN >= ref_len) {
        return;
    }

    string ref_base = ref.substr(i + delta_i, IDENTICAL_LEN);

    for (; j + delta_j + IDENTICAL_LEN < sv_len && delta_j <= 3 * MAX_SV_LEN; ++delta_j) {
        if (sv.compare(j + delta_j, IDENTICAL_LEN, ref_base) == 0) {
            break;
        }
    }

    if (delta_j > 3 * MAX_SV_LEN) {
        // something went wrong
        return;
    }

    while (delta_i > 0 && delta_j > 0 && ref[i + delta_i - 1] == sv[j + delta_j - 1])  {
        delta_i--;
        delta_j--;
    }

    pos = i;
    pos2 = j;
    len = delta_i;
    len2 = delta_j;
    i = i + delta_i;
    j = j + delta_j;
}

void find_svs(const string & ref, const string & sv, const string & id, vector<ans_t> & ans_vec, vector<ans_t> & unk_vec)
{
    int i, j, pos, len;
    i = j = 0;

    while (i < ref.size() && j < sv.size()) {
        if (ref[i] == sv[j]) {
            i++;
            j++;
            continue;
        } else {
            pos = len = -1;
            ans_t ans;
            ans.type = "UNK";
            if (check_ins(ref, sv, i, j, pos, len)) {
                if (check_dup(ref, sv, i, j, pos, len)) {
                    ans.type = "DUP";
                } else {
                    ans.type = "INS";
                }
            } else if (check_del(ref, sv, i, j, pos, len)) {
                ans.type = "DEL";
            } else if (check_inv(ref, sv, i, j, pos, len)) {
                ans.type = "INV";
            }

            if (ans.type == "UNK") {
                int pos2, len2;
                check_unk(ref, sv, i, j, pos, len, pos2, len2);
                ans.update(id, pos, len, id, pos2, len2);
                unk_vec.push_back(ans);
            } else {
                ans.update(id, pos, len);
                ans_vec.push_back(ans);
            }
        }
    }
}

void check_tra(unordered_map<string, string> & ref_map, unordered_map<string, string> & sv_map, vector<ans_t> & ans_vec,
               vector<ans_t> & unk_vec)
{
    for (auto i = 0; i < unk_vec.size(); ++i) {
        for (auto j = i + 1; j < unk_vec.size(); ++j) {
            if (unk_vec[i].id1 == unk_vec[j].id1) {
                // same chain
                continue;
            }
            ans_t ei = unk_vec[i];
            ans_t ej = unk_vec[j];
            string i_ref = ref_map[ei.id1].substr(ei.pos1, ei.len1);
            string j_sv = sv_map[ej.id1].substr(ej.pos2, ei.len2);
            if (i_ref == j_sv) {
                ans_t ans;
                ans.type = "TRA";
                ans.update(ei.id1, ei.pos1, ei.len1, ej.id1, ej.pos1, ej.len2);
                ans_vec.push_back(ans);
                break;
            }
        }
    }
}

void work(unordered_map<string, string> & ref_map, unordered_map<string, string> & sv_map)
{
    vector<ans_t> ans_vec, unk_vec;

    // find svs in each DNA chain
    for (const auto& x : ref_map) {
        string ref = x.second;
        string sv = sv_map[x.first];

        find_svs(ref, sv, x.first, ans_vec, unk_vec);
    }

    // check for TRA
    check_tra(ref_map, sv_map, ans_vec, unk_vec);

    // write output to file
    ofstream ofs;
    ofs.open(SVS_BED);
    for (const auto& x: ans_vec) {
        x.print(ofs);
    }
    ofs.close();
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
