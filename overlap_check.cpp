#include <iostream>
#include <vector>

using namespace std;


int main()
{
    vector<pair<unsigned long long, unsigned long long> > intervals;

    cout << "input intervals: " << endl;

    unsigned long long l, r;

    while (cin >> l >> r) {
        if (l == 0) break;
        intervals.push_back({l, r});
    }

    cout << "input finished" << endl;

    while (cin >> l >> r) {
        unsigned long long match = 0;
        for (auto p: intervals) {
            auto match_len = 100;
            if (p.first <= l && r <= p.second) {
                match = r - l;
            } else if (r <= p.first || p.second <= l) {
                match = 0;
            } else if (p.first <= r && r <= p.second) {
                match = r - p.first;
            } else if (p.first <= l && l <= p.second) {
                match = p.second - l;
            }
            if (match > match_len) {
                cout << "matched to " << p.first << " " << p.second << ": " << match << endl;
                break;
            }
        }
    }
}

