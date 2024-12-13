#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <numeric>
#include <cmath>
#include <random>
#include <map>

using namespace std;


mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
// similar results hold regardless of seed used

// standard linear programming solver
#define rep(i, from, to) for (int(i) = from; (i) < (to); ++(i))
#define all(x) x.begin(), (x).end()
#define sz(x) (int) (x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef long double T; // long double, Rational, double + mod<P>...
typedef vector<T> vd;
typedef vector<vd> vvd;

/**
 * Author: Stanford
 * Source: Stanford Notebook
 * License: MIT
 * Description: Solves a general linear maximization problem: maximize $c^T x$ subject to $Ax \le b$, $x \ge 0$.
 * Returns -inf if there is no solution, inf if there are arbitrarily good solutions, or the maximum value of $c^T x$ otherwise.
 * The input vector is set to an optimal  (or in the unbounded case, an arbitrary solution fulfilling the constraints).
 * Numerical stability is not guaranteed. For better performance, define variables such that $x = 0$ is viable.
 * Usage:
 * vvd A = {{1,-1}, {-1,1}, {-1,-2}};
 * vd b = {1,1,-4}, c = {-1,-1}, x;
 * T val = LPSolver(A, b, c).solve(x);
 * Time: O(NM * \#pivots), where a pivot may be e.g. an edge relaxation. O(2^n) in the general case.
 */

const T eps = 1e-8, inf = 1 / .0;
#define MP make_pair
#define ltj(X) if(s == -1 || MP(X[j],N[j]) < MP(X[s],N[s])) s=j

struct LPSolver {
    int m, n;
    vi N, B;
    vvd D;

    LPSolver(const vvd &A, const vd &b, const vd &c) :
            m(sz(b)), n(sz(c)), N(n + 1), B(m), D(m + 2, vd(n + 2)) {
        rep(i, 0, m) rep(j, 0, n) D[i][j] = A[i][j];
        rep(i, 0, m) {
            B[i] = n + i;
            D[i][n] = -1;
            D[i][n + 1] = b[i];
        }
        rep(j, 0, n) {
            N[j] = j;
            D[m][j] = -c[j];
        }
        N[n] = -1;
        D[m + 1][n] = 1;
    }

    void pivot(int r, int s) {
        T *a = D[r].data(), inv = 1 / a[s];
        rep(i, 0, m + 2) if (i != r && abs(D[i][s]) > eps) {
                T *b = D[i].data(), inv2 = b[s] * inv;
                rep(j, 0, n + 2) b[j] -= a[j] * inv2;
                b[s] = a[s] * inv2;
            }
        rep(j, 0, n + 2) if (j != s) D[r][j] *= inv;
        rep(i, 0, m + 2) if (i != r) D[i][s] *= -inv;
        D[r][s] = inv;
        swap(B[r], N[s]);
    }

    bool simplex(int phase) {
        int x = m + phase - 1;
        for (;;) {
            int s = -1;
            rep(j, 0, n + 1) if (N[j] != -phase) ltj(D[x]);
            if (D[x][s] >= -eps) return true;
            int r = -1;
            rep(i, 0, m) {
                if (D[i][s] <= eps) continue;
                if (r == -1 || MP(D[i][n + 1] / D[i][s], B[i])
                               < MP(D[r][n + 1] / D[r][s], B[r]))
                    r = i;
            }
            if (r == -1) return false;
            pivot(r, s);
        }
    }

    T solve(vd &x) {
        int r = 0;
        rep(i, 1, m) if (D[i][n + 1] < D[r][n + 1]) r = i;
        if (D[r][n + 1] < -eps) {
            pivot(r, n);
            if (!simplex(2) || D[m + 1][n + 1] < -eps) return -inf;
            rep(i, 0, m) if (B[i] == -1) {
                    int s = 0;
                    rep(j, 1, n + 1) ltj(D[i]);
                    pivot(i, s);
                }
        }
        bool ok = simplex(1);
        x = vd(n);
        rep(i, 0, m) if (B[i] < n) x[B[i]] = D[i][n + 1];
        return ok ? D[m][n + 1] : inf;
    }
};

// instance of fair division problem
struct instance {
    int n{}, m{};
    vector<vector<long double>> v;
    vector<long double> tot_u, ccs, efs, prop, ccs_guarantee, efs_guarantee, prop_guarantee;
    long double welfare = -1, ccs_theta = -1, ccs_sum = -1,
            efs_theta = -1, efs_sum = -1, prop_theta = -1;
    // should probably refactor to make a struct for each share

    friend ostream &operator<<(ostream &os, const instance &instance) {
        os << "n: " << instance.n << " m: " << instance.m << " welfare: " << instance.welfare << " ccs_sum: "
           << instance.ccs_sum << " ccs_theta: " << instance.ccs_theta << " efs_sum: "
           << instance.efs_sum << " efs_theta: " << instance.efs_theta << " prop_theta: " << instance.prop_theta;
        return os;
    }

    instance(int n, int m) {
        assert(n > 0 && m > 0);
        this->n = n;
        this->m = m;
        v = vector<vector<long double>>(n, vector<long double>(m));
        ccs.resize(n, -1);
        efs.resize(n, -1);
        tot_u.resize(n, -1);
        prop.resize(n, -1);
        ccs_guarantee.resize(n, -1);
        efs_guarantee.resize(n, -1);
        prop_guarantee.resize(n, -1);
    }

    instance(int n, int m, const vector<vector<long double>> &v) : instance(n, m) {
        this->v = v;
        for (int i = 0; i < n; ++i) {
            tot_u[i] = 0;
            for (int k = 0; k < m; ++k) {
                tot_u[i] += v[i][k];
            }
            prop[i] = tot_u[i] / n;
        }
        welfare = 0;
        for (int k = 0; k < m; ++k) {
            long double cur_mx = v[0][k];
            for (int i = 0; i < n; ++i) {
                cur_mx = max(cur_mx, v[i][k]);
            }
            welfare += cur_mx;
        }
    }

    void calc_ccs() {
        for (int i = 0; i < n; ++i) {
            // maximize $c^T x$ subject to $Ax \le b$, $x \ge 0$.
            vector<long double> c(m), b(n + m - 1), x;
            vector<vector<long double>> a(n + m - 1, vector<long double>(m, 0));

            // objective: maximize utility of bundle
            for (int j = 0; j < m; ++j) {
                c[j] = v[i][j];
            }
            // constraint: bundle satisfies proportionality for every other agent
            int ind = 0;
            for (int j = 0; j < n - 1; ++j, ++ind) {
                if (ind == i) {
                    ++ind;
                }
                for (int k = 0; k < m; ++k) {
                    a[j][k] = v[ind][k];
                    b[j] = tot_u[ind] / n;
                }
            }
            // constraint: at most 1 quantity of each item is allocated
            for (int j = n - 1; j < n + m - 1; ++j) {
                a[j][j - (n - 1)] = 1;
                b[j] = 1;
            }
            ccs[i] = LPSolver(a, b, c).solve(x);
        }
        ccs_sum = 0;
        for (int i = 0; i < n; ++i) {
            ccs_sum += ccs[i];
        }
    }

    long double
    calc_ratio(vector<long double> &shares) { // calculate approximation ratio \theta given a vector of share values
        // maximize $c^T x$ subject to $Ax \le b$, $x \ge 0$.
        vector<long double> c(n * m + 1), b(n + m), x; // x_{ik} = x[i * m + k], \theta = x[n * m]
        vector<vector<long double>> a(n + m, vector<long double>(n * m + 1, 0));

        // objective: maximize \theta
        c[n * m] = 1;
        // constraint: value of agent i's allocation is at least \theta \cdot \share_i
        for (int i = 0; i < n; ++i) {
            b[i] = 0;
            a[i][n * m] = shares[i];
            for (int k = 0; k < m; ++k) {
                a[i][i * m + k] = -v[i][k];
            }
        }
        // constraint: at most unit supply of each item is allocated
        for (int k = 0; k < m; ++k) {
            b[n + k] = 1;
            for (int i = 0; i < n; ++i) {
                a[n + k][i * m + k] = 1;
            }
        }
        return LPSolver(a, b, c).solve(x);
    }

    void calc_ccs_ratio() {
        if (ccs[0] == -1) {
            calc_ccs();
        }
        ccs_theta = calc_ratio(ccs);
        for (int i = 0; i < n; ++i) {
            ccs_guarantee[i] = ccs[i] * ccs_theta;
        }
    }

    void calc_efs() {
        for (int i = 0; i < n; ++i) {
            // maximize $c^T x$ subject to $Ax \le b$, $x \ge 0$.
            vector<long double> c(n * m, 0), b(n + m - 1, 0), x; // x_{ik} = x[i * m + k]
            vector<vector<long double>> a(n + m - 1, vector<long double>(n * m, 0));

            // objective: maximize utility of bundle
            for (int k = 0; k < m; ++k) {
                c[i * m + k] = v[i][k];
            }
            // constraint: for every agent j, they value their bundle weakly more than agent i's bundle
            int ind = 0;
            for (int j = 0; j < n; ++j) {
                if (j == i) {
                    continue;
                }
                for (int k = 0; k < m; ++k) {
                    a[ind][i * m + k] = v[j][k];
                    a[ind][j * m + k] = -v[j][k];
                }
                b[ind] = 0;
                ++ind;
            }
            // constraint: at most 1 quantity of each item is allocated
            for (int k = 0; k < m; ++k) {
                for (int j = 0; j < n; ++j) {
                    a[k + n - 1][j * m + k] = 1;
                }
                b[k + n - 1] = 1;
            }
            efs[i] = LPSolver(a, b, c).solve(x);
        }
        efs_sum = 0;
        for (int i = 0; i < n; ++i) {
            efs_sum += efs[i];
        }
    }

    void calc_efs_ratio() {
        if (efs[0] == -1) {
            calc_efs();
        }
        efs_theta = calc_ratio(efs);
        for (int i = 0; i < n; ++i) {
            efs_guarantee[i] = efs[i] * efs_theta;
        }
    }

    void calc_prop_ratio() {
        prop_theta = calc_ratio(prop);
        for (int i = 0; i < n; ++i) {
            prop_guarantee[i] = prop[i] * prop_theta;
        }
    }

    vector<long double> calc_p_efs_trial(long double delta) {
        vector<long double> p_efs(n);

        for (int i = 0; i < n; ++i) {
            vector<int> tot_agents(n);
            for (int j = 0; j < n; ++j) {
                tot_agents[j] = j;
            }
            tot_agents.erase(find(tot_agents.begin(), tot_agents.end(), i));
            shuffle(tot_agents.begin(), tot_agents.end(), rng);
            // randomly shuffle agents and take w to be the first floor((n-1)/delta)
            vector<int> w(tot_agents.begin(), tot_agents.begin() + min(n - 1, (int) floor((n - 1) / delta)));

            map<int, int> agent_ind;
            int cur_ind = 0;
            for (int j = 0; j < n; ++j) {
                if (find(w.begin(), w.end(), j) == w.end()) {
                    agent_ind[j] = cur_ind++;
                }
            }
            // maximize $c^T x$ subject to $Ax \le b$, $x \ge 0$.
            vector<long double> c(agent_ind.size() * m, 0), b(agent_ind.size() - 1 + m,
                                                              0), x; // x_{ik} = x[agent_ind[i] * m + k]
            vector<vector<long double>> a(agent_ind.size() - 1 + m, vector<long double>(agent_ind.size() * m, 0));
            // objective: maximize utility of bundle
            for (int k = 0; k < m; ++k) {
                c[agent_ind[i] * m + k] = v[i][k];
            }
            // constraint: for every agent j \not \in Z, they value their bundle weakly more than agent i's bundle
            int cons_ind = 0;
            for (int j = 0; j < n; ++j) {
                if (agent_ind.find(j) == agent_ind.end() || j == i) {
                    continue;
                }
                for (int k = 0; k < m; ++k) {
                    a[cons_ind][agent_ind[i] * m + k] = v[j][k];
                    a[cons_ind][agent_ind[j] * m + k] = -v[j][k];
                }
                b[cons_ind] = 0;
                ++cons_ind;
            }
            // constraint: at most 1 quantity of each item is allocated
            for (int k = 0; k < m; ++k) {
                for (int j = 0; j < n; ++j) {
                    if (agent_ind.find(j) != agent_ind.end() && j != i) {
                        a[cons_ind][agent_ind[j] * m + k] = 1;
                    }
                }
                a[cons_ind][agent_ind[i] * m + k] = (int) w.size() + 1;
                b[cons_ind] = 1;
                ++cons_ind;
            }
            p_efs[i] = LPSolver(a, b, c).solve(x);
        }
        return p_efs;
    }

    vector<long double> calc_p_efs(long double delta) {
        // estimates \pEFS value by random sampling 20 W_i
        int trials = 20;
        vector<long double> res(n);
        for (int i = 0; i < trials; ++i) {
            vector<long double> cur_p_efs = calc_p_efs_trial(delta);
            for (int j = 0; j < n; ++j) {
                res[j] += cur_p_efs[j];
            }
        }
        for (auto &i: res) {
            i /= trials;
        }
        return res;
    }

    void calc_all() {
        calc_ccs();
        calc_ccs_ratio();
        calc_efs();
        calc_efs_ratio();
        calc_prop_ratio();
    }
};

instance create_instance(const string &input_file) { // creates fair division from file input data
    // file format:
    // n m
    // n x m matrix of utilities
    ifstream in(input_file);
    int n, m;
    in >> n >> m;
    vector<vector<long double>> v(n, vector<long double>(m));
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < m; ++k) {
            in >> v[j][k];
        }
    }
    instance cur(n, m, v);
    return cur;
}

long double average(const vector<long double> &v) {
    return accumulate(v.begin(), v.end(), (long double) 0) / (long double) v.size();
}

long double stdev(const vector<long double> &v) {
    long double m = average(v);

    long double accum = 0.0;
    std::for_each(std::begin(v), std::end(v), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    long double stdev = sqrt(accum / (long double) (v.size() - 1));
    return stdev;
}

// runs all test cases in a directory and computes ccs and efs data
void run_all(const string &dir) {
    long double mx_ccs_sum_ratio = 0;
    long double mx_efs_sum_ratio = 0;
    long double mn_ccs_prop_ratio = 1e18;
    long double mn_efs_prop_ratio = 1e18;

    vector<long double> ccs_sum_ratios, efs_sum_ratios, ccs_prop_ratios, efs_prop_ratios, ccs_theta, efs_theta, prop_theta;

    int test_cases = 0;
    for (auto &entry: filesystem::directory_iterator(dir)) {
        if (entry.path().filename() == "output.txt") {
            continue;
        }
        ++test_cases;
        instance cur = create_instance(entry.path());
        cout << "Case " << test_cases << ": " << entry.path() << endl;
        cur.calc_all();
        cout << cur << endl;
        if (cur.welfare == 0) {
            continue;
        }
        long double ccs_sum_ratio = cur.ccs_sum / cur.welfare;
        long double efs_sum_ratio = cur.efs_sum / cur.welfare;
        mx_ccs_sum_ratio = max(mx_ccs_sum_ratio, ccs_sum_ratio);
        ccs_sum_ratios.push_back(ccs_sum_ratio);
        mx_efs_sum_ratio = max(mx_efs_sum_ratio, efs_sum_ratio);
        efs_sum_ratios.push_back(efs_sum_ratio);

        ccs_theta.push_back(cur.ccs_theta);
        efs_theta.push_back(cur.efs_theta);
        prop_theta.push_back(cur.prop_theta);

        for (int i = 0; i < cur.n; ++i) {
            cout << cur.prop[i] << " " << cur.ccs_guarantee[i] << " " << cur.efs_guarantee[i] << endl;
            if (cur.prop[i] == 0) {
                continue;
            }
            long double ccs_prop_ratio = cur.ccs_guarantee[i] / cur.prop[i];
            long double efs_prop_ratio = cur.efs_guarantee[i] / cur.prop[i];
            ccs_prop_ratios.push_back(ccs_prop_ratio);
            efs_prop_ratios.push_back(efs_prop_ratio);
            mn_ccs_prop_ratio = min(mn_ccs_prop_ratio, ccs_prop_ratio);
            mn_efs_prop_ratio = min(mn_efs_prop_ratio, efs_prop_ratio);
        }
        long double avg_ccs_sum_ratio = average(ccs_sum_ratios);
        long double avg_efs_sum_ratio = average(efs_sum_ratios);
        long double avg_ccs_prop_ratio = average(ccs_prop_ratios);
        long double avg_efs_prop_ratio = average(efs_prop_ratios);
        cout << "Maximum C(I) / SW(I) ratio: " << mx_ccs_sum_ratio << "\n";
        cout << "Average C(I) / SW(I) ratio: " << avg_ccs_sum_ratio << "\n";
        cout << "Stdev C(I) / SW(I) ratio: " << stdev(ccs_sum_ratios) << "\n";
        cout << "Maximum EF(I) / SW(I) ratio: " << mx_efs_sum_ratio << "\n";
        cout << "Average EF(I) / SW(I) ratio: " << avg_efs_sum_ratio << "\n";
        cout << "Stdev EF(I) / SW(I) ratio: " << stdev(efs_sum_ratios) << "\n";
        cout << "Minimum CCS allocation vs PROP allocation ratio: " << mn_ccs_prop_ratio << "\n";
        cout << "Average CCS allocation vs PROP allocation ratio: " << avg_ccs_prop_ratio << "\n";
        cout << "Stdev CCS allocation vs PROP allocation ratio: " << stdev(ccs_prop_ratios) << "\n";
        cout << "Minimum EFS allocation vs PROP allocation ratio: " << mn_efs_prop_ratio << "\n";
        cout << "Average EFS allocation vs PROP allocation ratio: " << avg_efs_prop_ratio << "\n";
        cout << "Stdev EFS allocation vs PROP allocation ratio: " << stdev(efs_prop_ratios) << "\n\n";
        cout << "Minimum PROP theta: " << *min_element(prop_theta.begin(), prop_theta.end()) << "\n";
        cout << "Average PROP theta: " << average(prop_theta) << "\n";
        cout << "Stdev PROP theta: " << stdev(prop_theta) << "\n";
        cout << "Minimum CCS theta: " << *min_element(ccs_theta.begin(), ccs_theta.end()) << "\n";
        cout << "Average CCS theta: " << average(ccs_theta) << "\n";
        cout << "Stdev CCS theta: " << stdev(ccs_theta) << "\n";
        cout << "Minimum EFS theta: " << *min_element(efs_theta.begin(), efs_theta.end()) << "\n";
        cout << "Average EFS theta: " << average(efs_theta) << "\n";
        cout << "Stdev EFS theta: " << stdev(efs_theta) << "\n";
    }
    cout << "\nprop_theta values:";
    for (auto &i: prop_theta) {
        cout << " " << i;
    }
    cout << "\n";
    cout << "ccs_theta values:";
    for (auto &i: ccs_theta) {
        cout << " " << i;
    }
    cout << "\n";
    cout << "efs_theta values:";
    for (auto &i: efs_theta) {
        cout << " " << i;
    }
    cout << "\n";
}

// runs all test cases in a directory and computes \pEFS values
void run_all2(const string &dir) {
    int test_cases = 0;
    int pre = -1, n = 20; // change depending on value of n
    vector<long double> delta_values;
    for (long double delta = 1; pre != 0; delta += 0.01) {
        int cur = floor((n - 1) / delta);
        if (cur != pre) {
            delta_values.push_back(delta);
            pre = cur;
        }
    }
    map<long double, vector<long double>> p_efs_theta;
    for (auto &entry: filesystem::directory_iterator(dir)) {
        if (entry.path().filename() == "output.txt") {
            continue;
        }
        ++test_cases;
        cout << "Test case: " << test_cases << endl;
        instance cur = create_instance(entry.path());
        for (auto &i: delta_values) {
            vector<long double> p_efs = cur.calc_p_efs(i);
            p_efs_theta[i].push_back(cur.calc_ratio(p_efs));
            cout << i << ": " << average(p_efs_theta[i]) << endl;
        }
        cout << "\n\n";


        if (test_cases % 50 == 0) {
            cout << "p_efs_theta values:\n{";
            for (auto &i: p_efs_theta) {
                cout << i.first << ": [";
                for (int j = 0; j < int(i.second.size()); ++j) {
                    cout << i.second[j];
                    if (j != (int) i.second.size() - 1) {
                        cout << ", ";
                    }
                }
                cout << "],\n";
            }
            cout << "}\n";
        }
    }
    cout << "p_efs_theta values:\n{";
    for (auto &i: p_efs_theta) {
        cout << i.first << ": [";
        for (int j = 0; j < int(i.second.size()); ++j) {
            cout << i.second[j];
            if (j != (int) i.second.size() - 1) {
                cout << ", ";
            }
        }
        cout << "],\n";
    }
    cout << "}\n";

    cout << "p_efs_theta average values:\n";
    for (auto &i: p_efs_theta) {
        cout << i.first << ": " << average(i.second) << "," << endl;
    }
    cout << "\n";
}

int main() {
//    run_all("data/n=25,m=75,uniform,sum=1000"); // specify test cases directory
    run_all2("data/n=20,m=20,yahoo");
    return 0;
}
