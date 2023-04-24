#include <bits/stdc++.h>
using namespace std;
#define MAXN 15
#define MAXBUF 1000005
#define MAXDB 100005
#define LL long long
#define double long double
const int inf = 1e9;

double f1(double x) { return 3 * x * x * x; }

double fT(int x) {
    // return exp(-x);
    return (double)1 / (double)pow(x, 0.7);
}

int n, sum_buf_size, num_request;

int SLA_act[MAXN];

double SLA_rate[MAXN];

int prior[MAXN], DB_size[MAXN];

struct BufSize {
    int min_buf_size, base_buf_size, max_buf_size;
};

BufSize buf_size_config[MAXN];

struct UserBuf {
    int us, buf_id;
};
UserBuf InCache[MAXBUF];

int where_user_buf[MAXN][MAXDB];
int num_buf_in[MAXN];

int las_ins_pos;

int ad_Time[MAXN][MAXDB];

list<pair<int, int>>::iterator mp[MAXBUF];

struct SetNode {
    int rat, where, time;
};
bool operator<(SetNode x, SetNode y) {
    if (x.rat != y.rat) return x.rat < y.rat;
    return x.where < y.where;
}

set<SetNode>::iterator it[MAXBUF];
struct LFRUCache {
    int pri_size;
    multiset<SetNode> pri;
    list<pair<int, int>> lru_list;
    pair<int, int> getLeastUsed() {
        if (lru_list.empty() == 0) {
            return lru_list.back();
        } else if (pri.empty()) {
            return make_pair(-1, -1);
        }
        return make_pair(pri.begin()->where, pri.begin()->time);
        // if (lru_list.empty()) {
        //     // cerr << "The lru_list is empty!" << endl;
        //     return make_pair(-1, -1);
        // }
        // return lru_list.back();
    }
    void add(int where, int request_id) {
        // cerr << "in " << ad_Time[InCache[where].us][InCache[where].buf_id]
        //      << " " << where << " " << request_id << endl;
        SetNode tmp =
            SetNode {ad_Time[InCache[where].us][InCache[where].buf_id], where,
                     request_id};
        S.insert(tmp);
        it[where] = S.lower_bound(tmp);
        // mp[where] = lru_list.begin();
        return;
    }
    void delPoint(int where) {
        set<SetNode>::iterator iit;
        if (where == 0) {
            iit = S.begin();
            // where = lru_list.back().first;
        } else {
            iit = it[where];
        }

        S.erase(iit);
        // lru_list.erase(mp[where]);
        // mp[where] = NULL;
    }
};
LFRUCache lru_one[MAXN];

struct LRUBase {
    int capacity, SLA;
    unordered_map<int, list<int>::iterator> cache;
    list<int> lru_list;
    LRUBase() {
        capacity = SLA = 0;
        cache.clear();
        lru_list.clear();
    }
    void add(int key) {
        auto it = cache.find(key);
        if (it != cache.end()) {
            lru_list.erase(it->second);
            lru_list.push_front(key);
            cache[key] = lru_list.begin();
        } else {
            ++SLA;
            if (cache.size() == capacity) {
                cache.erase(lru_list.back());
                lru_list.pop_back();
            }
            lru_list.push_front(key);
            cache[key] = lru_list.begin();
        }
    }
};

LRUBase lru_base[MAXN];

int cold_area[MAXN];

// struct LRUBase{

// };

void InitData() {
    cin >> n >> sum_buf_size >> num_request;
    for (int i = 1; i <= n; ++i) {
        cin >> prior[i];
    }
    for (int i = 1; i <= n; ++i) {
        cin >> DB_size[i];
    }
    for (int i = 1; i <= n; ++i) {
        cin >> buf_size_config[i].min_buf_size >>
            buf_size_config[i].base_buf_size >> buf_size_config[i].max_buf_size;
        lru_base[i].capacity = buf_size_config[i].base_buf_size;
    }
    cold_area[0] = 0;
    for (int i = 1; i <= n; ++i) {
        cold_area[i] = cold_area[i - 1] + buf_size_config[i].min_buf_size;
    }
}

// void PUT(int id) {
//     for (auto x : lru_one[id].lru_list) {
//         cerr << x.first << " " << x.second << endl;
//     }
//     return;
// }

void DelBuf(int where) {
    if (InCache[where].us == 0) {
        return;
    }
    int us = InCache[where].us;
    where_user_buf[us][InCache[where].buf_id] = 0;
    num_buf_in[us]--;
    // if (us == 1) {
    //     cerr << "del"
    //          << " " << us << " " << num_buf_in[1] << endl;
    //     PUT(1);
    // }
    lru_one[us].delPoint(where);

    InCache[where].us = 0;
    return;
}

void InsBuf(int request_id, int where, UserBuf req) {
    if (InCache[where].us > 0) {
        cerr << "Can not Ins because Cache is not NULL" << endl;
        return;
    }

    InCache[where] = req;
    where_user_buf[req.us][req.buf_id] = where;
    num_buf_in[req.us]++;
    lru_one[req.us].add(where, request_id);
    // if (req.us == 1) {
    //     cerr << request_id << " " << req.us << " " << num_buf_in[1] << endl;
    //     PUT(1);
    // }
    return;
}

bool CanReplaceBuf(int us, int req_us) {
    if (us == req_us) {
        if (num_buf_in[us] >= buf_size_config[us].min_buf_size &&
            num_buf_in[us] <= buf_size_config[us].max_buf_size) {
            return 1;
        } else {
            return 0;
        }
    }

    if (num_buf_in[req_us] >= buf_size_config[req_us].max_buf_size) {
        return 0;
    }
    if (us > 0 && num_buf_in[us] <= buf_size_config[us].min_buf_size) {
        return 0;
    }

    return 1;
}

int CalcWhereToIns_Traver(int request_id, UserBuf req) {
    if (where_user_buf[req.us][req.buf_id] > 0) {
        int where = where_user_buf[req.us][req.buf_id];
        return where;
    }

    while (1) {
        ++las_ins_pos;
        if (las_ins_pos > sum_buf_size || las_ins_pos < 1) las_ins_pos = 1;
        if (CanReplaceBuf(InCache[las_ins_pos].us, req.us)) {
            return las_ins_pos;
        }
    }
    return 0;
}

int CalcWhereToIns_LRUAll(int request_id, UserBuf req) {
    if (where_user_buf[req.us][req.buf_id] > 0) {
        int where = where_user_buf[req.us][req.buf_id];
        return where;
    }
    ++las_ins_pos;
    if (las_ins_pos <= sum_buf_size) {
        if (num_buf_in[req.us] <= buf_size_config[req.us].base_buf_size) {
            int where = las_ins_pos;
            return where;
        } else {
            las_ins_pos--;
        }
    }
    // todo:
    int Least_time = inf;
    int id = 0;
    for (int i = 1; i <= n; ++i) {
        if (num_buf_in[i] <= buf_size_config[i].base_buf_size) {
            continue;
        }
        pair<int, int> tmp = lru_one[i].getLeastUsed();
        if (tmp.first < 0) {
            continue;
        }
        if (CanReplaceBuf(i, req.us)) {
            if (Least_time > tmp.second) {
                Least_time = tmp.second;
                id = i;
            }
        }
    }
    if (id == 0) {
        for (int i = 1; i <= n; ++i) {
            pair<int, int> tmp = lru_one[i].getLeastUsed();
            if (tmp.first < 0) {
                continue;
            }
            if (CanReplaceBuf(i, req.us)) {
                if (Least_time > tmp.second) {
                    Least_time = tmp.second;
                    id = i;
                }
            }
        }
    }

    if (num_buf_in[req.us] >= buf_size_config[req.us].base_buf_size) {
        id = req.us;
    }
    pair<int, int> tmp = lru_one[id].getLeastUsed();
    return tmp.first;
}

int CalcWhereToIns_CheckSLA(int request_id, UserBuf req) {
    if (where_user_buf[req.us][req.buf_id] > 0) {
        int where = where_user_buf[req.us][req.buf_id];
        return where;
    }
    ++SLA_act[req.us];
    ++las_ins_pos;
    if (las_ins_pos <= sum_buf_size) {
        if (num_buf_in[req.us] <= buf_size_config[req.us].base_buf_size) {
            int where = las_ins_pos;
            return where;
        } else {
            las_ins_pos--;
        }
    }
    // todo:
    int Least_time = inf;
    int id = 0;
    for (int i = 1; i <= n; ++i) {
        if (num_buf_in[i] <= buf_size_config[i].base_buf_size) {
            continue;
        }
        pair<int, int> tmp = lru_one[i].getLeastUsed();
        if (tmp.first < 0) {
            continue;
        }
        if (CanReplaceBuf(i, req.us)) {
            if (Least_time > tmp.second) {
                Least_time = tmp.second;
                id = i;
            }
        }
    }
    if (id == 0) {
        for (int i = 1; i <= n; ++i) {
            pair<int, int> tmp = lru_one[i].getLeastUsed();
            if (tmp.first < 0) {
                continue;
            }
            if (CanReplaceBuf(i, req.us)) {
                if (Least_time > tmp.second) {
                    Least_time = tmp.second;
                    id = i;
                }
            }
        }
    }

    if (num_buf_in[req.us] >= buf_size_config[req.us].base_buf_size) {
        id = req.us;
    }
    pair<int, int> tmp = lru_one[id].getLeastUsed();
    return tmp.first;
}

int CalcWhereToIns_MinLoss(int request_id, UserBuf req) {
    if (where_user_buf[req.us][req.buf_id] > 0) {
        int where = where_user_buf[req.us][req.buf_id];
        return where;
    }
    ++SLA_act[req.us];
    ++las_ins_pos;
    if (las_ins_pos <= sum_buf_size) {
        if (num_buf_in[req.us] < buf_size_config[req.us].max_buf_size) {
            int where = las_ins_pos;
            return where;
        } else {
            las_ins_pos--;
        }
    }
    // todo:
    int Least_time = inf;
    double MinLoss = inf;
    LL UU = inf, DD = 0;
    int id = 0;
    // for (int i = 1; i <= n; ++i) {
    //     if (num_buf_in[i] <= buf_size_config[i].base_buf_size) {
    //         continue;
    //     }
    //     pair<int, int> tmp = lru_one[i].getLeastUsed();
    //     if (tmp.first < 0) {
    //         continue;
    //     }
    //     if (CanReplaceBuf(i, req.us)) {
    //         if (Least_time > tmp.second) {
    //             Least_time = tmp.second;
    //             id = i;
    //         }
    //     }
    // }
    // if (id > 0) {
    // cerr << "ok" << endl;
    // }
    // cerr << "??? " << id << endl;
    if (id == 0) {
        // cerr << "ok" << endl;
        for (int i = 1; i <= n; ++i) {
            pair<int, int> tmp = lru_one[i].getLeastUsed();
            // cerr << "? " << i << " " << tmp.first << " " << num_buf_in[i]
            //      << endl;
            if (tmp.first < 0) {
                continue;
            }
            // cerr << i << " " << req.us << endl;
            if (CanReplaceBuf(i, req.us)) {
                // cerr << lru_base[i].SLA << endl;
                LL sla_if =
                    (max(SLA_act[i] + 1, lru_base[i].SLA) - lru_base[i].SLA);
                LL sla_now =
                    (max(SLA_act[i], lru_base[i].SLA) - lru_base[i].SLA);
                LL diff = sla_if * sla_if - sla_now * sla_now;
                LL U = diff * (LL)3 * (LL)prior[i];
                LL D = (LL)lru_base[i].SLA * (LL)lru_base[i].SLA;
                double diffLoss =
                    (double)((LL)diff * (LL)3 * (LL)prior[i]) /
                    (double)((LL)lru_base[i].SLA * (LL)lru_base[i].SLA);

                double rate_if = (double)(max(SLA_act[i] + 1, lru_base[i].SLA) -
                                          lru_base[i].SLA) /
                                 (double)lru_base[i].SLA;
                double nowLoss = f1(rate_if) * (double)prior[i],
                       lasLoss = f1(SLA_rate[i]) * (double)prior[i],
                       diffLoss2 = nowLoss - lasLoss;
                // if (fabs(diffLoss - diffLoss2) > 1e-2) {
                //     cerr << "? " << request_id << " "
                //          << fabs(diffLoss - diffLoss2) << endl;
                // }

                if (MinLoss > diffLoss * (double)fT(request_id - tmp.second)) {
                    MinLoss = diffLoss * (double)fT(request_id - tmp.second);
                    id = i;
                }
                // if (UU * D > U * DD) {
                //     UU = U;
                //     DD = D;
                //     // MinLoss = diffLoss;
                //     id = i;
                // }
            }
        }
    }
    // cerr << id << " " << num_buf_in[req.us] << " "
    //      << buf_size_config[req.us].min_buf_size << " "
    //      << buf_size_config[req.us].max_buf_size << " "
    //      << CanReplaceBuf(req.us, req.us) << endl;
    // if (num_buf_in[req.us] >= buf_size_config[req.us].base_buf_size) {
    //     id = req.us;
    // }

    pair<int, int> tmp = lru_one[id].getLeastUsed();
    return tmp.first;
}

void ReplaceBuf(int request_id, int where, UserBuf req) {
    // if (!CanReplaceBuf(where, req.us)) {
    //     cerr << "Can not Replace Buf" << endl;
    //     return;
    // }
    DelBuf(where);
    InsBuf(request_id, where, req);
    return;
}

void SolveRequest(int request_id, UserBuf req) {
    int where_ins_buf = CalcWhereToIns_MinLoss(request_id, req);
    ReplaceBuf(request_id, where_ins_buf, req);
    cout << where_ins_buf << endl;
}

void MAIN() {
    UserBuf req;

    // Init
    InitData();
    for (int i = 1; i <= n; ++i) {
        lru_one[i].pri_size = buf_size_config[i].base_buf_size / 3 + 1;
    }
    // solve requests
    for (int request_id = 1; request_id <= num_request; ++request_id) {
        cin >> req.us >> req.buf_id;
        ++ad_Time[req.us][req.buf_id];
        // cerr << request_id << endl;
        // cerr << req.us << " " << req.buf_id << endl;
        lru_base[req.us].add(req.buf_id);
        // cerr << "?" << endl;
        SolveRequest(request_id, req);
        // cerr << "??" << endl;
        SLA_rate[req.us] = (double)(max(SLA_act[req.us], lru_base[req.us].SLA) -
                                    lru_base[req.us].SLA) /
                           (double)lru_base[req.us].SLA;
        // cerr << "now " << request_id << endl;
        // for (int j = 1; j <= n; ++j) {
        //     cerr << "id " << j << " "
        //          << "SLA  base & act " << lru_base[j].SLA << " " <<
        //          SLA_act[j]
        //          << endl;
        // }
        // for (int j = 1; j <= n; ++j) {
        //     cerr << "rate " << SLA_rate[j] << " " << prior[j] << " "
        //          << num_buf_in[j] << " " << buf_size_config[j].min_buf_size
        //          << endl;
        // }
    }
    for (int j = 1; j <= n; ++j) {
        cerr << "id " << j << " "
             << "SLA  base & act " << lru_base[j].SLA << " " << SLA_act[j]
             << endl;
    }
    for (int i = 1; i <= n; ++i) {
        SLA_rate[i] =
            (double)(max(SLA_act[i], lru_base[i].SLA) - lru_base[i].SLA) /
            (double)lru_base[i].SLA;
    }
    for (int j = 1; j <= n; ++j) {
        cerr << "rate " << SLA_rate[j] << " " << prior[j] << " "
             << num_buf_in[j] << " " << buf_size_config[j].min_buf_size << endl;
    }
    double cost = 0;
    for (int i = 1; i <= n; ++i) {
        cost += f1(SLA_rate[i]) * (double)prior[i];
    }
    cerr << "cost " << cost << endl;
    return;
}

int main() {
    // freopen("data/8.in", "r", stdin);
    // freopen("data/8.out", "w", stdout);
    srand(time(0));
    std::ios::sync_with_stdio(false);
    int ttt = 1;
    // cin >> ttt;
    for (int i = 1; i <= ttt; ++i) {
        MAIN();
    }
    return 0;
}