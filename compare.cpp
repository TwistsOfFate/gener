#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <tuple>
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

typedef unsigned long long ull;

void compare(const string & ref, const string & query, const align_res_t& alignment, vector<string> & answers)
{
    static const int HASH_LEN = 5, COMPARE_LEN = 50;

    ull l1 = alignment.l1, l2 = alignment.l2;
    string base = ref.substr(l1, query.length() * 3 / 2);
    auto base_len = base.length();

    unordered_multimap<ull, ull> hashmap;
    make_hashmap(base, hashmap, HASH_LEN);

    for (auto i = l2; i < query.length(); i += COMPARE_LEN / 2) {
        align_res_t res = align("ref", "query", hashmap, query.substr(i, COMPARE_LEN), HASH_LEN, base_len);
//        cout << "compare: [" << l1 + res.l1 << "," << l1 + res.r1 << "]:\t" << res.score << "/" << COMPARE_LEN << endl;
    }


}

vector<string> find_answers(string & ref, vector<string> & concats)
{
    vector<string> answers;
    auto ref_len = ref.length();
    unordered_multimap<ull, ull> hashmap;
    make_hashmap(ref, hashmap, LEN_D);

    double score_sum = 0.0;
    ull query_len_sum = 0;
    for (auto i = 0; i < concats.size(); ++i) {
        cout << "find_answers: initial match " << i+1 << "/" << concats.size() << "\r";
        cout.flush();

        align_res_t align_res = align("ref", "query", hashmap, concats[i], LEN_D, ref_len);
        compare(ref, concats[i], align_res, answers);
        cout << endl << endl;

//        if (i == 10) while (1);

        /* DEBUG */
        score_sum += align_res.score;
        query_len_sum += concats[i].length();
    }
    cout << endl;

    cout << "Average matching accuracy: " << score_sum / query_len_sum << endl;
    cout << "Average matching score: " << score_sum / concats.size() << endl;


    return answers;
}

