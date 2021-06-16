#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>

#define REF_FASTA "../ref.fasta"
#define SV_FASTA "../task2_sv.fasta"
#define LONG_FASTA "../long.fasta"
#define SVS_BED "../sv.bed"
#define IDENTICAL_LEN 50
#define MAX_SV_LEN 1005
#define MIN_SV_LEN 50

using namespace std;

class ans_t {
public:
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

void init(unordered_map<string, string> & ref_map, unordered_map<string, string> & sv_map)
{
    ifstream ref_if, long_if;
    string line;
    int count[10005], ct2[15];
    int tot_lines = 0;

    for (int i = 0; i < 10001; ++i) {
        count[i] = 0;
    }
    for (int i = 0; i < 15; ++i) {
        ct2[i] = 0;
    }

    long_if.open(LONG_FASTA);
    while (getline(long_if, line)) {
        string title = line.substr(1);
        if (!getline(long_if, line)) {
            break;
        }
        count[line.size()]++;
        tot_lines++;
    }
    long_if.close();

    for (int i = 0; i < 10001; ++i) {
        if (count[i]) {
            ct2[i / 1000] += count[i];
        }
    }

    for (int i = 0; i < 11; ++i) {
        cout << i * 1000 << "-" << (i + 1) * 1000 << ": " << ct2[i] << endl;
    }
    cout << "Total lines: " << tot_lines << endl;
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

int main() {
    // for timing ------------------------------------------------------------------------------------------------------
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();
    //------------------------------------------------------------------------------------------------------------------

    unordered_map<string, string> ref_map, sv_map;
    vector<ans_t> ans_vec, unk_vec;

    init(ref_map, sv_map);

    auto t2 = high_resolution_clock::now();


    // for timing ------------------------------------------------------------------------------------------------------

    auto t3 = high_resolution_clock::now();

    duration<double, std::milli> ms1 = t2 - t1;
    duration<double, std::milli> ms2 = t3 - t2;

//    cout << ms1.count() << "ms + " << ms2.count() << "ms" << endl;
    //------------------------------------------------------------------------------------------------------------------

    return 0;
}
