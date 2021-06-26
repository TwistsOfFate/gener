#ifndef GENER_DEFINES_H
#define GENER_DEFINES_H

#include <string>
#include <iostream>

#define REF_FASTA "../data/ref.fasta"
#define SV_FASTA "../data/sv.fasta"
#define LONG_FASTA "../data/long.fasta"
#define SVS_BED "../sv.bed"
#define IDENTICAL_LEN 50
#define MAX_SV_LEN 1005
#define MIN_SV_LEN 50
#define LEN_A 12
#define LEN_B 100
#define LEN_C 10
#define LEN_G 100

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
    double score;

    explicit align_res_t(const std::string& s="", const std::string& t="", unsigned long long x=0, unsigned long long y=0,
                unsigned long long a=0, unsigned long long b=0, unsigned long long c=0, unsigned long long d=0, double e=0)
    {
        this->id1 = s; this->id2 = t;
        this->len1 = x; this->len2 = y;
        this->l1 = a; this->r1 = b; this->l2 = c; this->r2 = d;
        this->score = e;
    }

    void update(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long d, double e)
    {
        this->l1 = a;
        this->r1 = b;
        this->l2 = c;
        this->r2 = d;
        this->score = e;
    }

    void clear()
    {
        this->id1 = ""; this->id2 = "";
        this->len1 = 0; this->len2 = 0;
        this->l1 = 0; this->r1 = 0; this->l2 = 0; this->r2 = 0;
        this->score = 0;
    }

    bool operator < (const align_res_t & rhs) const
    {
        return l1 < rhs.l1;
    }

    unsigned long long extend_len() const
    {
        return r1 + len2 - r2;
    }
};

#endif //GENER_DEFINES_H
