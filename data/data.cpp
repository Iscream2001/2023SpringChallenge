#include <bits/stdc++.h>
using namespace std;
#define MAXN 15
#define MAXBUF 1000005
#define MAXDB 100005
int n, Q, M, DBlim;
int prior[MAXN];
int DB[MAXN];

int DBconfig[MAXN][4];

int cnt[MAXN];
int las_id[MAXN][MAXBUF];

void MAIN() {
    n = 10;
    Q = 100;
    M = 100000;

    cout << n << " " << Q << " " << M << endl;

    for (int i = 1; i <= n; ++i) {
        prior[i] = i;
    }
    random_shuffle(prior + 1, prior + 1 + n);
    for (int i = 1; i <= n; ++i) {
        cout << prior[i] << " ";
    }
    cout << endl;

    DBlim = 50;

    for (int i = 1; i <= n; ++i) {
        DB[i] = rand() % DBlim + 1;
        cout << DB[i] << " ";
    }
    cout << endl;

    for (int i = 1; i <= n; ++i) {
        DBconfig[i][0] = rand() % DB[i] + 1;
        DBconfig[i][1] = rand() % DB[i] + 1;
        DBconfig[i][2] = rand() % DB[i] + 1;
        sort(DBconfig[i], DBconfig[i] + 3);
        int mid = (DBconfig[i][0] + DBconfig[i][2]) / 2;
        DBconfig[i][1] = max(DBconfig[i][1], mid);
        cout << DBconfig[i][0] << " " << DBconfig[i][1] << " " << DBconfig[i][2]
             << endl;
    }

    for (int i = 1; i <= M; ++i) {
        int u = rand() % (n * 2) + 1;
        u = min(u, n);
        int id = rand() % DB[u] + 1;
        if (cnt[u] > DBconfig[u][1]) {
            int op = rand() % 100;
            if (op > 30) {
                id = las_id[u][cnt[u] - rand() % (DBconfig[u][1] / 2 + 1)];
            }
        }
        cout << u << " " << id << endl;
        ++cnt[u];
        las_id[u][cnt[u]] = id;
    }
}
int main() {
    // freopen("1.in", "r", stdin);
    // freopen("7.in", "w", stdout);
    srand(time(0));
    std::ios::sync_with_stdio(false);
    int ttt = 1;
    // cin >> ttt;
    for (int i = 1; i <= ttt; ++i) {
        MAIN();
    }
    return 0;
}