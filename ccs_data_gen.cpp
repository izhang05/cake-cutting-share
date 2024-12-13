#include <bits/stdc++.h>

using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

long long rnd(long long l, long long r) {
    return uniform_int_distribution<long long>(l, r)(rng);
}

long double rnd(long double l, long double r) {
    return uniform_real_distribution<long double>(l, r)(rng);
}

// Generate data for Bernoulli model
int main() {
    int data_files = 200;
    bernoulli_distribution b(.5);

    for (int i = 0; i < data_files; ++i) {
        string output_file = "data/n=25,m=75,binary,p=.5/in" + to_string(i + 1) + ".txt";
        if (filesystem::exists(output_file)) {
            cout << "File \"" << output_file << "\" exists. Press \"y\" to continue.\n";
            string tmp;
            cin >> tmp;
            if (tmp != "y") {
                exit(0);
            }
        }
        ofstream out(output_file);
        int n = 25, m = 75;
        out << n << " " << m << "\n";
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < m; ++k) {
                out << b(rng) << " \n"[k == m - 1];
            }
        }
        out.close();
    }
}

// Generate data for intrinsic value model
//int main() {
//    int data_files = 200;
//    for (int i = 0; i < data_files; ++i) {
//        string output_file = "data/n=25,m=75,intrinsic_value,d=.3/in" + to_string(i + 1) + ".txt";
//        if (filesystem::exists(output_file)) {
//            cout << "File \"" << output_file << "\" exists. Press \"y\" to continue.\n";
//            string tmp;
//            cin >> tmp;
//            if (tmp != "y") {
//                exit(0);
//            }
//        }
//        ofstream out(output_file);
//        int n = 25, m = 75;
//        out << n << " " << m << "\n";
//        vector<long double> v(m); // intrinsic value
//        for (auto &j: v) {
//            j = rnd((long double) 0, 1);
//        }
//
//        for (int j = 0; j < n; ++j) {
//            for (int k = 0; k < m; ++k) {
//                out << v[k] + rnd((long double) 0, 0.3) << " \n"[k == m - 1];
//            }
//        }
//        out.close();
//    }
//}

// Generate data for uniform distribution model
//
//int main() {
//    int data_files = 200;
//    for (int i = 0; i < data_files; ++i) {
//        string output_file = "data/n=25,m=75,uniform,sum=1000/in" + to_string(i + 1) + ".txt";
//        if (filesystem::exists(output_file)) {
//            cout << "File \"" << output_file << "\" exists. Press \"y\" to continue.\n";
//            string tmp;
//            cin >> tmp;
//            if (tmp != "y") {
//                exit(0);
//            }
//        }
//        ofstream out(output_file);
//        int n = 25, m = 75;
//        out << n << " " << m << "\n";
//        int sum = 1000;
//        vector<int> cur(m + sum - 1);
//        for (int j = 0; j < m + sum - 1; ++j) {
//            cur[j] = j;
//        }
//        for (int j = 0; j < n; ++j) {
//            shuffle(cur.begin(), cur.end(), rng);
//            vector<int> cur_nums(cur.begin(), cur.begin() + m - 1);
//            sort(cur_nums.begin(), cur_nums.end());
//            vector<int> v(m);
//            v[0] = cur_nums[0] - 1;
//            for (int k = 1; k < m - 1; ++k) {
//                v[k] = cur_nums[k] - cur_nums[k - 1] - 1;
//            }
//            v[m - 1] = sum + m - 1 - cur_nums[m - 2];
//            for (int k = 0; k < m; ++k) {
//                out << v[k] << " \n"[k == m - 1];
//            }
//        }
//        out.close();
//    }
//}

