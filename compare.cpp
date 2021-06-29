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
using std::ofstream;    using std::abs;

typedef unsigned long long ull;

void find_dup(vector<double> const & ref_cover, const align_res_t& alignment, int COMPARE_LEN, int SCAN_TIMES,
              vector<ans_t> & suspects)
{
    const double LO = COMPARE_LEN * 0.02, HI = COMPARE_LEN * SCAN_TIMES * 1.2;
    ull ans_l = 0, ans_r = 0;
    string ans_stat = "NONE";
    for (auto i = 0; i < ref_cover.size(); ++i) {
        if (ref_cover[i] > HI) {
            ans_r = i;
            if (ans_stat == "NONE") {
                ans_stat = "DUP";
                ans_l = i;
            }
        } else if (ref_cover[i] < LO) {
            ans_r = i;
            if (ans_stat == "NONE") {
                ans_stat = "DEL";
                ans_l = i;
            }
        } else {
            if (ans_r - ans_l > 100 && ans_r - ans_l < 1100 && ans_l > 0) {
                bool valid = false;
                double score = 0.0;
                for (auto j = ans_l; j <= ans_r; ++j) {
                    score += ref_cover[j];
                }
                if (ans_stat == "DUP" && score / (ans_r - ans_l) > HI) {
                    valid = true;
                } else if (ans_stat == "DEL" && score / (ans_r - ans_l) < LO) {
                    valid = true;
                }
                if (valid) {
                    for (auto k = ans_l + alignment.pos; k <= ans_r + alignment.pos; ++k) {
                        suspects[k].votes[ans_stat]++;
                    }
                }
            }
            ans_l = ans_r = 0;
            ans_stat = "NONE";
        }
    }
}

void identify_misses(vector<double> const & ref_cover, vector<double> const & query_cover,
                     const align_res_t& alignment, int COMPARE_LEN, int SCAN_TIMES, ull query_len,
                     vector<align_res_t> & misses)
{
    pair<ull, ull> ref_interval{0, 0}, query_interval{0, 0};
    double ref_score{99999.0}, query_score{99999.0};

    ofstream ofs;
    ofs.open("../logs/id_misses.log", std::ios_base::app);

    const double LO_REF = COMPARE_LEN * 0.05 * SCAN_TIMES;
    ull ans_l = 0, ans_r = 0;
    string ans_stat = "NONE";
    for (auto i = 0; i < ref_cover.size(); ++i) {
        if (ref_cover[i] < LO_REF) {
            ans_r = i;
            if (ans_stat == "NONE") {
                ans_stat = "DEL";
                ans_l = i;
            }
        } else {
            if (i - ans_r > 100) {
                if (ans_r - ans_l > 100 && ans_r - ans_l < 1100 && ans_l > 0) {
                    bool valid = false;
                    double score = 0.0;
                    for (auto j = ans_l; j <= ans_r; ++j) {
                        score += ref_cover[j];
                    }
                    if (ans_stat == "DEL" && score / (ans_r - ans_l) < LO_REF) {
                        valid = true;
                    }
                    if (valid) {
                        if (ref_interval.second - ref_interval.first < ans_r - ans_l) {
                            ref_interval = make_pair(ans_l, ans_r);
                            ref_score = score;
                        }
                    }
                }
                ans_l = ans_r = 0;
                ans_stat = "NONE";
            }
        }
    }

    const double LO_QUERY = COMPARE_LEN * 0.4 * SCAN_TIMES;
    ans_l = 0, ans_r = 0;
    ans_stat = "NONE";
    for (auto i = 0; i < query_cover.size(); ++i) {
        if (query_cover[i] < LO_QUERY) {
            ans_r = i;
            if (ans_stat == "NONE") {
                ans_stat = "DEL";
                ans_l = i;
            }
        } else {
            if (i - ans_r > 100) {
                if (ans_r - ans_l > 100 && ans_r - ans_l < 1100 && ans_l > 0) {
                    bool valid = false;
                    double score = 0.0;
                    for (auto j = ans_l; j <= ans_r; ++j) {
                        score += query_cover[j];
                    }
                    if (ans_stat == "DEL" && score / (ans_r - ans_l) < LO_QUERY) {
                        valid = true;
                    }
                    if (valid) {
                        if (query_interval.second - query_interval.first < ans_r - ans_l) {
                            query_interval = make_pair(ans_l, ans_r);
                            query_score = score;
                        }
                    }
                }
                ans_l = ans_r = 0;
                ans_stat = "NONE";
            }
        }
    }

    auto ref_diff = ref_interval.second - ref_interval.first;
    auto query_diff = query_interval.second - query_interval.first;
    auto abs_len_diff = ref_diff > query_diff ? ref_diff - query_diff : query_diff - ref_diff;
    auto abs_pos_diff = ref_interval.first > query_interval.first ? ref_interval.first - query_interval.first : query_interval.first - ref_interval.first;

    if (ref_diff > 0 && query_diff > 0 && abs_len_diff < min(ref_diff, query_diff) * 0.3 && abs_pos_diff < min(ref_diff, query_diff) * 0.4) {
//        for (auto k = query_interval.first + alignment.pos; k <= query_interval.second + alignment.pos; ++k) {
//            suspects[k].votes["INV"]++;
//        }
        ans_l = query_interval.first;
        ans_r = query_interval.second;
        ofs << "segment length: " << query_len << "\tscore: " << alignment.score << endl;
        ofs << "INV " << ans_l << " " << ans_r << " " << query_score / (ans_r - ans_l);
        ofs << "\t" << ans_l + alignment.pos << " " << ans_r + alignment.pos << endl;
        ofs << endl;
    } else if (ref_diff > 0 && query_diff == 0) {
//        for (auto k = ref_interval.first + alignment.pos; k <= ref_interval.second + alignment.pos; ++k) {
//            suspects[k].votes["DEL"]++;
//        }
        ans_l = ref_interval.first;
        ans_r = ref_interval.second;
        ofs << "segment length: " << query_len << "\tscore: " << alignment.score << endl;
        ofs << "DEL " << ans_l << " " << ans_r << " " << ref_score / (ans_r - ans_l);
        ofs << "\t" << ans_l + alignment.pos << " " << ans_r + alignment.pos << endl;
        ofs << endl;
    } else if (ref_diff == 0 && query_diff > 0) {
//        for (auto k = query_interval.first + alignment.pos; k <= query_interval.second + alignment.pos; ++k) {
//            suspects[k].votes["INS"]++;
//        }
        ans_l = query_interval.first;
        ans_r = query_interval.second;
        ofs << "segment length: " << query_len << "\tscore: " << alignment.score << endl;
        ofs << "INS " << ans_l << " " << ans_r << " " << query_score / (ans_r - ans_l);
        ofs << "\t" << ans_l + alignment.pos << " " << ans_r + alignment.pos << endl;
        ofs << endl;
    }


}

void compare(const string & ref, const string & query, const align_res_t& alignment, vector<ans_t> & suspects)
{
    static const int HASH_LEN = 5, COMPARE_LEN = 50, SCAN_TIMES = 2;
    const ull BASE_LEN = query.length() * 2;

    string base = ref.substr(alignment.pos, BASE_LEN);

    unordered_multimap<ull, ull> hashmap;
    make_hashmap(base, hashmap, HASH_LEN);

    vector<double> ref_cover(BASE_LEN + COMPARE_LEN, 0.0);
    vector<double> query_cover(query.length() + COMPARE_LEN, 0.0);
    for (auto i = 0; i < query.length(); i += COMPARE_LEN / SCAN_TIMES) {
        align_res_t res = align("ref", "query", hashmap, query.substr(i, COMPARE_LEN), HASH_LEN, BASE_LEN);
        for (auto j = res.l1; j < res.r1 + HASH_LEN; ++j) {
            ref_cover[j] += res.score;
        }
        for (auto j = i; j < i + COMPARE_LEN; ++j) {
            query_cover[j] += res.score;
        }
    }

    find_dup(ref_cover, alignment, COMPARE_LEN, SCAN_TIMES, suspects);

    vector<align_res_t> misses;
    identify_misses(ref_cover, query_cover, alignment, COMPARE_LEN, SCAN_TIMES, query.length(), misses);

}

void find_answers(string & ref, vector<string> & concats, string const & id, vector<string> & answers)
{
    auto ref_len = ref.length();
    unordered_multimap<ull, ull> hashmap;
    make_hashmap(ref, hashmap, LEN_D);

    double score_sum = 0.0;
    ull query_len_sum = 0;
    vector<ans_t> suspects(ref_len + 1000);
    for (auto i = 0; i < concats.size(); ++i) {
        cout << "find_answers: " << i+1 << "/" << concats.size() << "\r";
        cout.flush();

        align_res_t align_res = align("ref", "query", hashmap, concats[i], LEN_D, ref_len);
        compare(ref, concats[i], align_res, suspects);

        /* DEBUG */
        score_sum += align_res.score;
        query_len_sum += concats[i].length();
    }
    cout << endl;

    /* DEBUG */
    cout << "Average matching accuracy: " << score_sum / query_len_sum << endl;
    cout << "Average matching score: " << score_sum / concats.size() << endl;


    ofstream ofs;
    ofs.open("../logs/answers.txt", std::ios_base::app);
    string stat;
    auto fst = 0, lst = 0;
    for (auto i = 0; i < ref_len; ++i) {
        string maj = suspects[i].majority();
        if (stat.empty()) {
            if (suspects[i].votes[maj] >= 2 || maj == "INV" && suspects[i].votes[maj] >= 1) {
                stat = maj;
                fst = lst = i;
            }
        } else {
            if (maj.empty() || maj != stat) {
                if (i - lst > 200) {
                    if (lst - fst > 50 && lst - fst < 1500) {
                        answers.push_back(stat + " " + id + " " + to_string(fst) + " " + to_string(lst));
//                        ofs << stat << " " << id << " " << fst << " " << lst << endl;
                    }
                    stat.clear();
                    fst = lst = 0;
                }
            } else if (maj == stat) {
                lst = i;
            }
        }
    }
    ofs.close();

}

