#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <bitset>

#define REF_FASTA "../data/ref.fasta"
#define SV_FASTA "../data/sv.fasta"
#define LONG_FASTA "../data/long.fasta"
#define SVS_BED "../sv.bed"
#define IDENTICAL_LEN 50
#define MAX_SV_LEN 1005
#define MIN_SV_LEN 50
#define LEN_A 15

using std::cout;        using std::endl;
using std::string;      using std::unordered_map;
using std::unordered_multimap;
using std::ifstream;    using std::vector;
using std::ofstream;    using std::bitset;

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

void read_fasta(unordered_map<string, string> & ref_map, unordered_map<string, string> & long_map)
{
    ifstream ref_if, long_if;
    string line;

    ref_if.open(REF_FASTA);
    while (getline(ref_if, line)) {
        string title = line.substr(1);
        if (!getline(ref_if, line)) {
            break;
        }
        ref_map[title] = line;
    }
    ref_if.close();

    long_if.open(LONG_FASTA);
    while (getline(long_if, line)) {
        string title = line.substr(1);
        if (!getline(long_if, line)) {
            break;
        }
        long_map[title] = line;
    }
    long_if.close();
}

void make_hashmap(unordered_map<string, string> & src_map, unordered_multimap<unsigned long long, id_pos_t> & res_map, int hash_len)
{
    unsigned long long mask = 0;
    if (hash_len > 32) {
        cout << "make_hashmap: possible collisions" << endl;
        mask = ~0;
    } else {
        for (auto i = 0; i < hash_len * 2; ++i) {
            mask = (mask << 1) + 1;
        }
    }
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
        auto i = base_len;
        for (i = 0; i < hash_len; ++i) {
            hash = ((hash << 2) | dict[base[i]]) & mask;
        }
        for (; i < base_len; ++i) {
            res_map.insert({hash, id_pos_t(id, i - hash_len)});
            hash = ((hash << 2) | dict[base[i]]) & mask;
        }
        res_map.insert({hash, id_pos_t(id, i - hash_len)});
    }

    cout << "make_hashmap: done." << endl;
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

void check_tra(unordered_map<string, string> & ref_map, unordered_map<string, string> & sv_map, vector<ans_t> & ans_vec, vector<ans_t> & unk_vec)
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
    unordered_map<string, string> ref_map, long_map;
    unordered_multimap<unsigned long long, id_pos_t> ref_hashmap;

    read_fasta(ref_map, long_map);
    make_hashmap(ref_map, ref_hashmap, 11);

    unordered_map<unsigned long long, int> key_count;

    int multi_occur_count = 0;
    int total_keys = 0;

    for (auto const & x: ref_hashmap) {
        if (key_count.find(x.first) == key_count.end()) {
            total_keys++;
        }
        key_count.insert({x.first, 0});
        key_count[x.first]++;
    }


    for (auto const & x: key_count) {
        if (x.second > 2) {
            multi_occur_count++;
//            cout << x.first << ":\t" << x.second << endl;
//            auto it = ref_hashmap.equal_range(x.first);
//            for (auto y = it.first; y != it.second; ++y) {
//                y->second.print(cout);
//                cout << endl;
//            }
        }
    }

    cout << "multiple occurs: " << multi_occur_count << "/" << total_keys << "=" << double(multi_occur_count) / total_keys << endl;
    cout << "average occurs: " << double(ref_hashmap.size()) / total_keys << endl;

//    exam_data(ref_map, long_map);
//    work(ref_map, sv_map);


    return 0;
}
