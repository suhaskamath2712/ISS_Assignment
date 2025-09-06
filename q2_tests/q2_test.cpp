#include <bits/stdc++.h>
#include "../user_code.h"
#include <sys/resource.h>  // getrusage
#include <unistd.h>        // getpid
using namespace std;

// Helper to get current memory usage in KB (resident set size)
long getMemoryUsageKB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__) && defined(__MACH__)
    return usage.ru_maxrss / 1024; // macOS reports bytes
#else
    return usage.ru_maxrss;        // Linux reports KB
#endif
}

// Utility: safely parse a list from a string like "[1,2,3]"
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

// Utility: parse nested list, e.g. "[[1,2],[3,4]]"
vector<vector<int>> parse_nested_int_list(const string &s) {
    vector<vector<int>> res;
    vector<int> current;
    string num;
    for (char c : s) {
        if (isdigit(c) || c == '-') num.push_back(c);
        else {
            if (!num.empty()) {
                current.push_back(stoi(num));
                num.clear();
            }
            if (c == ']') {
                if (!current.empty()) {
                    res.push_back(current);
                    current.clear();
                }
            }
        }
    }
    return res;
}

// Parse the big input dictionary string into components
void parse_input(const string &input_str,
                 vector<int> &preorder,
                 vector<int> &inorder,
                 vector<vector<int>> &leafParcels,
                 vector<vector<int>> &queries) {
    // crude parsing: search for keys and substrings
    auto get_section = [&](const string &key) {
        size_t pos = input_str.find(key);
        if (pos == string::npos) return string("");
        size_t start = input_str.find('[', pos);
        size_t end = input_str.find(']', start);
        // For nested, we need to balance brackets
        int depth = 0;
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
    };

    string pre_str = get_section("preorder");
    string in_str = get_section("inorder");
    string leaf_str = get_section("leafParcels");
    string query_str = get_section("queries");

    preorder = parse_int_list(pre_str);
    inorder = parse_int_list(in_str);
    leafParcels = parse_nested_int_list(leaf_str);
    queries = parse_nested_int_list(query_str);
}

// Format output vector<int> to string
string format_output(const vector<int> &res) {
    stringstream ss;
    ss << "[";
    for (size_t i = 0; i < res.size(); i++) {
        ss << res[i];
        if (i + 1 < res.size()) ss << ", ";
    }
    ss << "]";
    return ss.str();
}

int main() {
    ifstream fin_inputs("q2_inputs.csv");
    ifstream fin_outputs("q2_outputs.csv");
    ofstream fout("q2_timings.csv");

    if (!fin_inputs || !fin_outputs) {
        cerr << "Error: could not open q2_inputs.csv or q2_outputs.csv\n";
        return 1;
    }

    string line_in, line_out;
    getline(fin_inputs, line_in);  // skip header
    getline(fin_outputs, line_out);

    fout << "num_nodes,num_queries,time_microseconds,memory_kb,correct\n";

    int case_num = 0;
    while (getline(fin_inputs, line_in) && getline(fin_outputs, line_out)) {
        case_num++;

        // Remove quotes if present
        if (!line_in.empty() && line_in.front() == '"' && line_in.back() == '"')
            line_in = line_in.substr(1, line_in.size()-2);
        if (!line_out.empty() && line_out.front() == '"' && line_out.back() == '"')
            line_out = line_out.substr(1, line_out.size()-2);

        // Parse input
        vector<int> preorder, inorder;
        vector<vector<int>> leafParcels, queries;
        parse_input(line_in, preorder, inorder, leafParcels, queries);

        // Expected output
        vector<int> expected = parse_int_list(line_out);

    // Measure memory before, run & time, then measure memory after
    long mem_before = getMemoryUsageKB();

    auto queries_copy = queries; // pass by value
    auto start = chrono::high_resolution_clock::now();
    vector<int> got = question_two(preorder, inorder, leafParcels, queries_copy);
    auto end = chrono::high_resolution_clock::now();

    long mem_after = getMemoryUsageKB();
    long mem_used = mem_after - mem_before;
    if (mem_used < 0) mem_used = 0;

    long long elapsed = chrono::duration_cast<chrono::microseconds>(end - start).count();

        // Compare
        bool correct = (got == expected);

        // Write result
    fout << preorder.size() << "," << queries.size() << "," << elapsed << "," << mem_used << "," << (correct ? "1" : "0") << "\n";
    }

    cout << "Finished Q2 testing. Results in q2_timings.csv\n";
    return 0;
}
