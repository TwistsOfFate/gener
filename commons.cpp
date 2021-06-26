#include <string>
#include <unordered_map>
#include <iostream>

using std::string;          using std::unordered_map;
using std::cout;            using std::endl;

unsigned long long get_mask_from_hash_len(int hash_len)
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

string inverse(const string & src)
{
    string res;
    static unordered_map<char, string> comp = {{'A', "T"}, {'T', "A"}, {'C', "G"}, {'G', "C"}};

    for (int i = src.size() - 1; i >= 0; --i) {
        res.append(comp[src[i]]);
    }
    return res;
}

