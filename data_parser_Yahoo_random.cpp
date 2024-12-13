#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <random>
#include <chrono>

using namespace std;

// n = agents (advertisers), m = items (advertisements)
const int mxn = 10475, mxm = 1000;

long double agent_vals[mxn][mxm];
//set id on each item and agent for easy access
map<int, int> mp_agent, mp_item;
int id_agent = 0, id_item = 0;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

long long rnd(long long l, long long r) {
    return uniform_int_distribution<long long>(l, r)(rng);
}

long double rnd(long double l, long double r) {
    return uniform_real_distribution<long double>(l, r)(rng);
}

//records the bid of each advertiser on the advertisements
map<pair<int, int>, long double> bid;

map<int, map<int, long double>> agent_bids; // agent_bids[i][j] is agent i's bid for item j
map<int, map<int, long double>> item_bids; // item_bids[i][j] is agent j's bid for item i

void read_data() {
    ifstream fin("ydata-ysm-advertiser-bids-v1_0.txt");
    string useless;

    int cnt = 0;

    //read data
    while (fin >> useless) {
        fin >> useless;
        int item, agent;
        long double price;
        fin >> item >> agent >> price;
        fin >> useless;
        if (price == -1) {
            continue;
        }
        if (!mp_agent.count(agent)) {
            mp_agent[agent] = id_agent++;
        }
        if (!mp_item.count(item)) {
            mp_item[item] = id_item++;
        }
        bid[make_pair(mp_agent[agent], mp_item[item])] = price;
        agent_bids[mp_agent[agent]][mp_item[item]] = price;
        item_bids[mp_item[item]][mp_agent[agent]] = price;
        agent_vals[mp_agent[agent]][mp_item[item]] = price;

        cnt++;
//        if (cnt > 1e7) break;
    }

}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    read_data();
    cout << "done reading data" << endl;

    int data_files = 200;
    for (int file_num = 0; file_num < data_files; ++file_num) {
        string output_file = "data/n=20,m=20,yahoo/in" + to_string(file_num + 1) + ".txt";
        int n = 20, m = 20;
        vector<int> cur_agents(n);
        vector<int> cur_items(m);

        int starting_item = rnd(0ll, id_item - m);
        for (int i = 0; i < m; ++i) {
            cur_items[i] = i + starting_item;
        }

        map<int, map<int, long double>> cur_agent_bids; // item_bids[i][j] is agent i's bid for item j
        for (auto &i: cur_items) {
            for (auto &j: item_bids[i]) {
                cur_agent_bids[j.first][i] = j.second;
            }
        }

        vector<pair<int, int>> items_per_agent; // vector<pair<item_id, # items that the agent likes>>
        for (int i = 0; i < id_agent; ++i) {
            items_per_agent.emplace_back(i, cur_agent_bids[i].size());
        }
        sort(items_per_agent.begin(), items_per_agent.end(), [](auto left, auto right) {
            return left.second == right.second ? left.first < right.first : left.second < right.second;
        }); // sort by second element
        reverse(items_per_agent.begin(), items_per_agent.end());
        // items are sorted by the number of agents that like that item

        // take the most popular agents
        for (int i = 0; i < m; ++i) {
            cur_agents[i] = items_per_agent[i].first;
        }

        vector<vector<long double>> value(n, vector<long double>(m));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                value[i][j] = agent_bids[cur_agents[i]][cur_items[j]];
            }
        }

        ofstream out(output_file);
        out << n << " " << m << "\n";
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                out << value[i][j] << " \n"[j == m - 1];
            }
        }
        out.close();
    }
}
