#include <bits/stdc++.h>
using namespace std;
#define MAXN 15
#define MAXBUF 1000005
#define MAXDB 100005
#define LL long long
#define double long double
const int inf = 1e9;

double f1(double x) { return pow(x, 2.9); }
double f2(double x) { return 3 * x * x; }
double fT(int x) {
    // return exp(-x);
    return (double)1 / (double)pow(x, 0.24);
}
double fR(int x) { return pow(x, 0.18); }

int n, sum_buf_size, num_request;

int SLA_act[MAXN];
int num_ask[MAXN];

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

pair<int, int> TM[15];

int where_user_buf[MAXN][MAXDB];
int num_buf_in[MAXN];

int las_ins_pos;

int ad_Time[MAXN][MAXDB];

struct SetNode {
    int rat, where, time;
};
bool operator<(SetNode x, SetNode y) {
    if (x.rat != y.rat) return x.rat < y.rat;
    if (x.time != y.time) return x.time < y.time;
    return x.where < y.where;
}

struct SetNode2 {
    int rat, where, time;
};
bool operator<(SetNode2 x, SetNode2 y) {
    if (x.time != y.time) return x.time < y.time;
    // if (x.rat != y.rat) return x.rat < y.rat;
    return x.where < y.where;
}

set<SetNode2>::iterator mp[MAXBUF];
set<SetNode>::iterator it[MAXBUF];
bool is_pri[MAXBUF];
struct LFRUCache {
    double pri_size;
    int lru_op;
    multiset<SetNode> pri;
    multiset<SetNode2> lru;

    void resize() {
        if (pri.size() + lru.size() == 0) return;
        // cerr << "? " << pri.size() << " " << lru.size() << " " << pri_size
        //      << endl;
        if (pri_size == 0) {
            int t = min((int)pri.size(), 50);
            t = pri.size();
            while (t--) {
                // cerr << "??? " << pri.size() - 1 << " " << (pri.size() +
                // lru.size())
                //      << " " << pri_size << endl;
                SetNode tmp = *pri.begin();
                // del pri.begin
                pri.erase(it[tmp.where]);

                is_pri[tmp.where] = 0;

                // put begin to lrulist
                SetNode2 tmp2 = SetNode2 {tmp.rat, tmp.where, tmp.time};
                lru.insert(tmp2);
                mp[tmp2.where] = lru.lower_bound(tmp2);
            }
        }
        if (pri_size == 1) {
            int t = min((int)lru.size(), 50);
            while (t--) {
                // cerr << "??? " << pri.size() - 1 << " " << (pri.size() +
                // lru.size())
                //      << " " << pri_size << endl;
                multiset<SetNode2>::iterator iit;
                if (lru_op == 0) {
                    iit = lru.begin();
                } else {
                    iit = lru.end();
                    iit--;
                }
                // iit--;
                SetNode2 tmp = *iit;

                // del pri.begin
                lru.erase(mp[tmp.where]);

                is_pri[tmp.where] = 1;

                // put begin to lrulist
                SetNode tmp2 = SetNode {tmp.rat, tmp.where, tmp.time};
                pri.insert(tmp2);
                it[tmp2.where] = pri.lower_bound(tmp2);
            }
        }

        // cerr << "??" << endl;
    }

    pair<int, int> getLeastUsed() {
        if (lru.empty() == 0) {
            if (lru_op > 0) {
                return make_pair(lru.begin()->where, lru.begin()->time);
            } else {
                return make_pair(lru.rbegin()->where, lru.rbegin()->time);
            }
        } else if (pri.empty()) {
            return make_pair(-1, -1);
        } else {
            if (lru_op > 0) {
                return make_pair(pri.begin()->where, pri.begin()->time);
            } else {
                return make_pair(pri.rbegin()->where, pri.rbegin()->time);
            }
        }
        // if (lru_list.empty()) {
        //     // cerr << "The lru_list is empty!" << endl;
        //     return make_pair(-1, -1);
        // }
        // return lru_list.back();
    }
    void add(int where, int request_id) {
        // cerr << "in " << ad_Time[InCache[where].us][InCache[where].buf_id]
        //      << " " << where << " " << request_id << endl;
        is_pri[where] = 1;
        SetNode tmp =
            SetNode {ad_Time[InCache[where].us][InCache[where].buf_id], where,
                     request_id};
        pri.insert(tmp);
        it[where] = pri.lower_bound(tmp);
        resize();

        return;
    }
    void delPoint(int where) {
        set<SetNode>::iterator iit;
        set<SetNode2>::iterator iit2;
        if (where == 0) {
            if (lru.empty() == 0) {
                if (lru_op > 0) {
                    iit2 = lru.begin();
                    lru.erase(iit2);
                    resize();
                    return;
                } else {
                    iit2 = lru.end();
                    iit2--;
                    lru.erase(iit2);
                    resize();
                    return;
                }
            } else {
                if (lru_op > 0) {
                    iit = pri.begin();
                    pri.erase(iit);
                    resize();
                    return;
                } else {
                    iit = pri.end();
                    iit--;
                    pri.erase(iit);
                    resize();
                    return;
                }
            }

            // where = lru_list.back().first;
        } else {
            if (is_pri[where]) {
                iit = it[where];
                pri.erase(iit);
                resize();
                return;
            } else {
                iit2 = mp[where];
                lru.erase(iit2);
                resize();
            }
        }

        // S.erase(iit);
        // lru_list.erase(mp[where]);
        // mp[where] = NULL;
    }
};
LFRUCache lru_one[MAXN];

vector<pair<int, int>> jilu[MAXN];

struct Simu {
    int be_lru, be_lfu, be_mru, re_lru, re_lfu, re_mru, id;
    int rat_pri[MAXDB], las_pri[MAXDB];
    int rat_lru[MAXDB], las_lru[MAXDB];
    int rat_mru[MAXDB], las_mru[MAXDB];
    multiset<SetNode> pri, mru;
    multiset<SetNode2> lru;
    Simu() {
        be_lfu = be_lru = be_mru = 0;
        re_lru = re_lfu = re_mru = 0;
        for (int i = 1; i <= 100000; ++i) {
            rat_lru[i] = rat_pri[i] = rat_mru[i] = 0;
            las_lru[i] = las_pri[i] = las_mru[i] = -1;
        }
        pri.clear();
        lru.clear();
        mru.clear();
    }
    int ask_mru() {
        set<SetNode>::iterator iit;
        for (int i = be_mru; i < jilu[id].size(); ++i) {
            int lim = jilu[id][i].first;
            int x = jilu[id][i].second;
            while (mru.size() > lim) {
                iit = mru.end();
                iit--;
                mru.erase(iit);
            }
            rat_mru[x]++;
            SetNode tmp, tmp2;
            if (las_mru[x] >= 0) {
                tmp = SetNode {rat_mru[x] - 1, x, las_mru[x]};
                iit = mru.lower_bound(tmp);
                if (iit != mru.end()) {
                    tmp2 = *iit;
                    if (tmp2.rat == tmp.rat && tmp2.time == tmp.time &&
                        tmp2.where == tmp.where) {
                        mru.erase(iit);
                        --re_mru;
                    }
                }
            }
            ++re_mru;
            tmp = SetNode {rat_mru[x], x, i};
            las_mru[x] = i;
            mru.insert(tmp);
        }
        be_mru = jilu[id].size();
        return re_mru;
    }
    int ask_lru() {
        set<SetNode2>::iterator iit;
        for (int i = be_lru; i < jilu[id].size(); ++i) {
            int lim = jilu[id][i].first;
            int x = jilu[id][i].second;
            while (lru.size() > lim) {
                iit = lru.begin();
                lru.erase(iit);
            }
            rat_lru[x]++;
            SetNode2 tmp, tmp2;
            if (las_lru[x] >= 0) {
                tmp = SetNode2 {rat_lru[x] - 1, x, las_lru[x]};
                iit = lru.lower_bound(tmp);
                if (iit != lru.end()) {
                    tmp2 = *iit;
                    if (tmp2.rat == tmp.rat && tmp2.time == tmp.time &&
                        tmp2.where == tmp.where) {
                        lru.erase(iit);
                        --re_lru;
                    }
                }
            }
            ++re_lru;
            tmp = SetNode2 {rat_lru[x], x, i};
            las_lru[x] = i;
            lru.insert(tmp);
        }
        be_lru = jilu[id].size();
        return re_lru;
    }
    int ask_lfu() {
        set<SetNode>::iterator iit;
        for (int i = be_lfu; i < jilu[id].size(); ++i) {
            int lim = jilu[id][i].first;
            int x = jilu[id][i].second;
            while (pri.size() > lim) {
                iit = pri.begin();
                pri.erase(iit);
            }
            rat_pri[x]++;
            SetNode tmp, tmp2;
            if (las_pri[x] >= 0) {
                tmp = SetNode {rat_pri[x] - 1, x, las_pri[x]};
                iit = pri.lower_bound(tmp);
                if (iit != pri.end()) {
                    tmp2 = *iit;
                    if (tmp2.rat == tmp.rat && tmp2.time == tmp.time &&
                        tmp2.where == tmp.where) {
                        pri.erase(iit);
                        --re_lfu;
                    }
                }
            }
            ++re_lfu;
            tmp = SetNode {rat_pri[x], x, i};
            las_pri[x] = i;
            pri.insert(tmp);
        }
        be_lfu = jilu[id].size();
        return re_lfu;
    }
} simu[MAXN];

struct LRUBase {
    int capacity, SLA;

    unordered_map<int, list<pair<int, int>>::iterator> cache;
    list<pair<int, int>> lru_list;
    LRUBase() {
        capacity = SLA = 0;
        cache.clear();
        lru_list.clear();
    }
    void add(int key, int time) {
        auto it = cache.find(key);
        if (it != cache.end()) {
            lru_list.erase(it->second);
            lru_list.push_front(make_pair(key, time));
            cache[key] = lru_list.begin();

        } else {
            ++SLA;
            if (cache.size() == capacity) {
                cache.erase(lru_list.back().first);
                lru_list.pop_back();
            }
            lru_list.push_front(make_pair(key, time));
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
    TM[us] = lru_one[us].getLeastUsed();
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

    TM[req.us] = lru_one[req.us].getLeastUsed();

    // if (req.us == 1) {
    //     cerr << request_id << " " << req.us << " " << num_buf_in[1] <<
    //     endl; PUT(1);
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
    for (int i = 1; i <= n; ++i) {
        if (SLA_act[i] >= lru_base[i].SLA) {
            continue;
        }
        if (num_buf_in[i] <= buf_size_config[i].base_buf_size) {
            continue;
        }
        // pair<int, int> tmp = lru_one[i].getLeastUsed();
        pair<int, int> tmp = TM[i];
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
            if (num_buf_in[i] <= buf_size_config[i].base_buf_size) {
                continue;
            }
            // pair<int, int> tmp = lru_one[i].getLeastUsed();
            pair<int, int> tmp = TM[i];
            // cerr << "? " << i << " " << tmp.first << " " << num_buf_in[i]
            //      << endl;
            if (tmp.first < 0) {
                continue;
            }
            // cerr << i << " " << req.us << endl;
            if (CanReplaceBuf(i, req.us)) {
                // cerr << lru_base[i].SLA << endl;

                double rate_if = (double)(max(SLA_act[i] + 1, lru_base[i].SLA) -
                                          lru_base[i].SLA) /
                                 (double)lru_base[i].SLA;
                double nowLoss = f1(rate_if) * (double)prior[i],
                       lasLoss = f1(SLA_rate[i]) * (double)prior[i],
                       diffLoss2 = nowLoss - lasLoss;
                double vaa = diffLoss2 * (double)fT(request_id - tmp.second) *
                             fR(ad_Time[InCache[tmp.first].us]
                                       [InCache[tmp.first].buf_id]);
                if (MinLoss > vaa) {
                    MinLoss = vaa;
                    id = i;
                }
            }
        }
    }
    // if (id > 0) {
    // cerr << "ok" << endl;
    // }
    // cerr << "??? " << id << endl;
    if (id == 0) {
        // cerr << "ok" << endl;
        for (int i = 1; i <= n; ++i) {
            // pair<int, int> tmp = lru_one[i].getLeastUsed();
            pair<int, int> tmp = TM[i];
            // cerr << "? " << i << " " << tmp.first << " " << num_buf_in[i]
            //      << endl;
            if (tmp.first < 0) {
                continue;
            }
            // cerr << i << " " << req.us << endl;
            if (CanReplaceBuf(i, req.us)) {
                // cerr << lru_base[i].SLA << endl;

                double rate_if = (double)(max(SLA_act[i] + 1, lru_base[i].SLA) -
                                          lru_base[i].SLA) /
                                 (double)lru_base[i].SLA;
                double nowLoss = f1(rate_if) * (double)prior[i],
                       lasLoss = f1(SLA_rate[i]) * (double)prior[i],
                       diffLoss2 = nowLoss - lasLoss;
                double vaa = diffLoss2 * (double)fT(request_id - tmp.second) *
                             fR(ad_Time[InCache[tmp.first].us]
                                       [InCache[tmp.first].buf_id]);
                if (MinLoss > vaa) {
                    MinLoss = vaa;
                    id = i;
                }
            }
        }
    }
    // cerr << id << " " << num_buf_in[req.us] << " "
    //      << buf_size_config[req.us].min_buf_size << " "
    //      << buf_size_config[req.us].max_buf_size << " "
    //      << CanReplaceBuf(req.us, req.us) << endl;
    if (num_buf_in[req.us] >= buf_size_config[req.us].base_buf_size) {
        id = req.us;
    }

    pair<int, int> tmp_ = lru_one[id].getLeastUsed();
    return tmp_.first;
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
        lru_one[i].pri_size = 0;
        lru_one[i].lru_op = 1;
        // lru_one[i].pri_size = 1;
    }
    // solve requests
    int kuai = sqrt(num_request) + 1;
    // int mo[MAXN];
    for (int i = 1; i <= n; ++i) {
        simu[i].id = i;
    }

    for (int i = 1; i <= n; ++i) {
        TM[i] = lru_one[i].getLeastUsed();
    }

    for (int request_id = 1; request_id <= num_request; ++request_id) {
        cin >> req.us >> req.buf_id;
        ++ad_Time[req.us][req.buf_id];
        ++num_ask[req.us];

        jilu[req.us].push_back(
            make_pair(lru_one[req.us].pri.size() + lru_one[req.us].lru.size(),
                      req.buf_id));

        if (1) {
            int res_lru = simu[req.us].ask_lru();
            int res_lfu = simu[req.us].ask_lfu();
            int res_mru = simu[req.us].ask_mru();
            if (res_mru < res_lru && res_mru < res_lfu) {
                lru_one[req.us].pri_size = 1;
                lru_one[req.us].lru_op = 0;
            } else if (res_lru > res_lfu) {
                lru_one[req.us].pri_size = 1;
                lru_one[req.us].lru_op = 1;
            } else {
                lru_one[req.us].pri_size = 0;
                lru_one[req.us].lru_op = 1;
            }
        }
        lru_base[req.us].add(req.buf_id, num_ask[req.us]);
        SolveRequest(request_id, req);
        SLA_rate[req.us] = (double)(max(SLA_act[req.us], lru_base[req.us].SLA) -
                                    lru_base[req.us].SLA) /
                           (double)lru_base[req.us].SLA;
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
        cost += f2(SLA_rate[i]) * (double)prior[i];
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