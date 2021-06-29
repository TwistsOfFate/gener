#ifndef GENER_DEFINES_H
#define GENER_DEFINES_H

#include <string>
#include <iostream>
#include <unordered_map>

#define REF_FASTA "../data/ref.fasta"
#define SV_FASTA "../data/sv.fasta"
#define LONG_FASTA "../data/long.fasta"
#define SVS_BED "../sv.bed"
#define IDENTICAL_LEN 50
#define MAX_SV_LEN 1005
#define MIN_SV_LEN 50
#define SEED_LEN_A 12
#define LEN_B 100
#define SEED_LEN_B 10
#define LEN_G 200
#define LEN_D 12

struct id_pos_t {
    std::string id;
    unsigned long long pos;

    id_pos_t(std::string id, unsigned long long pos): id(std::move(id)), pos(pos) {}

    void print(std::ostream & os = std::cout) const
    {
        os << id << "\t" << pos;
    }
};

struct align_res_t {
    std::string id1, id2;
    unsigned long long len1, len2;
    unsigned long long l1, r1, l2, r2;
    unsigned long long pos;
    double score;

    explicit align_res_t(const std::string& s="", const std::string& t="", unsigned long long x=0, unsigned long long y=0,
                unsigned long long a=0, unsigned long long b=0, unsigned long long c=0, unsigned long long d=0, double e=0)
    {
        this->id1 = s; this->id2 = t;
        this->len1 = x; this->len2 = y;
        this->l1 = a; this->r1 = b; this->l2 = c; this->r2 = d;
        this->pos = l1 - l2;
        this->score = e;
    }

    void update(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long d, double e)
    {
        this->l1 = a;
        this->r1 = b;
        this->l2 = c;
        this->r2 = d;
        this->pos = l1 - l2;
        this->score = e;
    }

    void clear()
    {
        this->id1 = ""; this->id2 = "";
        this->len1 = 0; this->len2 = 0;
        this->l1 = 0; this->r1 = 0; this->l2 = 0; this->r2 = 0;
        this->pos = 0;
        this->score = 0;
    }

    bool operator < (const align_res_t & rhs) const
    {
        return pos < rhs.pos;
    }

    unsigned long long extend_len() const
    {
        return r1 + len2 - r2;
    }
};

struct ans_t {
    std::unordered_map<std::string, unsigned long long> votes;
    ans_t()
    {
        votes["INS"] = 0;
        votes["DEL"] = 0;
        votes["DUP"] = 0;
        votes["INV"] = 0;
        votes["TRA"] = 0;
    }

    std::string majority()
    {
        std::string maj = "";
        unsigned long long mx = 0;
        for (const auto& x: votes) {
            if (x.second > mx) {
                maj = x.first;
                mx = x.second;
            }
        }
        if (maj == "DEL" && votes["INV"]) {
             maj = "INV";
        }
        return maj;
    }
};

#endif //GENER_DEFINES_H
