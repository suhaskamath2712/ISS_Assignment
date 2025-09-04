#include <bits/stdc++.h>
#include "../user_code.h"
using namespace std;

// Utility: parse list of integers from string "[1,2,3]"
vector<int> parse_int_list(const string &s) {
    vector<int> res;
    string num;
    for (char c : s) {
        if (isdigit(c) || c == '-') num.push_back(c);
        else {
            if (!num.empty()) {
                res.push_back(stoi(num));
                num.clear();
            }
        }
    }
    return res;
}

// Utility: parse edges [[u,v,w], ...]
vector<vector<int>> parse_edges(const string &s) {
    vector<vector<int>> res;
    vector<int> cur;
    string num;
    for (char c : s) {
        if (isdigit(c) || c == '-') num.push_back(c);
        else {
            if (!num.empty()) {
                cur.push_back(stoi(num));
                num.clear();
            }
            if (c == ']') {
                if (!cur.empty()) {
                    res.push_back(cur);
                    cur.clear();
                }
            }
        }
    }
    return res;
}

// Extract substring for key
string extract_section(const string &input_str, const string &key) {
    size_t pos = input_str.find(key);
    if (pos == string::npos) return "";
    size_t start = input_str.find('[', pos);
    if (start == string::npos) return "";
    int depth = 0;
    size_t end = start;
    for (size_t i = start; i < input_str.size(); i++) {
        if (input_str[i] == '[') depth++;
        else if (input_str[i] == ']') {
            depth--;
            if (depth == 0) {
                end = i;
                break;
            }
        }
    }
    return input_str.substr(start, end - start + 1);
}

int main() {
    ifstream fin_inputs("q3_inputs.csv");
    ifstream fin_outputs("q3_outputs.csv");
    ofstream fout("q3_timings.csv");

    if (!fin_inputs || !fin_outputs) {
        cerr << "Error: could not open q3_inputs.csv or q3_outputs.csv\n";
        return 1;
    }

    string line_in, line_out;
    getline(fin_inputs, line_in); // skip headers
    getline(fin_outputs, line_out);

    fout << "num_nodes,num_edges,num_metros,time_microseconds,correct\n";

    int case_num = 0;
    while (getline(fin_inputs, line_in) && getline(fin_outputs, line_out)) {
        case_num++;
        if (!line_in.empty() && line_in.front() == '"' && line_in.back() == '"')
            line_in = line_in.substr(1, line_in.size()-2);
        if (!line_out.empty() && line_out.front() == '"' && line_out.back() == '"')
            line_out = line_out.substr(1, line_out.size()-2);

        string edges_str = extract_section(line_in, "edges");
        string metros_str = extract_section(line_in, "metro_cities");

        vector<vector<int>> edges = parse_edges(edges_str);
        vector<int> metro_cities = parse_int_list(metros_str);

        // Expected output
        long long expected = stoll(line_out);

        // Run & time
        auto start = chrono::high_resolution_clock::now();
        long long got = question_three(edges, metro_cities);
        auto end = chrono::high_resolution_clock::now();
        long long elapsed = chrono::duration_cast<chrono::microseconds>(end - start).count();

        int n = 0;
        for (auto &e : edges) n = max(n, max(e[0], e[1]));

        bool correct = (got == expected);
        fout << n << "," << edges.size() << "," << metro_cities.size() << "," << elapsed << "," << (correct ? "1" : "0") << "\n";
    }

    cout << "Finished Q3 testing. Results in q3_timings.csv\n";
    return 0;
}
