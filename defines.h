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
#define LEN_B 200

struct id_pos_t {
    std::string id;
    unsigned long long pos;

    id_pos_t(std::string id, unsigned long long pos): id(std::move(id)), pos(pos) {}

    void print(std::ostream & os = std::cout) const
    {
        os << id << "\t" << pos;
    }
};

#endif //GENER_DEFINES_H
