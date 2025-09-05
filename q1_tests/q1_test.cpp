#include <bits/stdc++.h>
#include "../user_code.h"   // your function is here
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

// Helper: parse string like "[[1, 10], [2, 5]]" -> vector<vector<int>>
vector<vector<int>> parse_parcels(const string& s) {
    vector<vector<int>> parcels;
    vector<int> current;
    string num;
    for (char c : s) {
        if (isdigit(c) || c == '-') {
            num.push_back(c);
        } else {
            if (!num.empty()) {
                current.push_back(stoi(num));
                num.clear();
            }
            if (c == ']') {
                if (!current.empty()) {
                    parcels.push_back(current);
                    current.clear();
                }
            }
        }
    }
    return parcels;
}

// Helper: format result into string like "[[id, weight], ...]"
string format_output(const vector<vector<int>>& res) {
    stringstream ss;
    ss << "[";
    for (size_t i = 0; i < res.size(); i++) {
        ss << "[" << res[i][0] << ", " << res[i][1] << "]";
        if (i + 1 < res.size()) ss << ", ";
    }
    ss << "]";
    return ss.str();
}

int main() {
    ifstream fin_inputs("q1_inputs.csv");
    ifstream fin_outputs("q1_outputs.csv");
    ofstream fout("q1_timings.csv");

    if (!fin_inputs || !fin_outputs) {
        cerr << "Error: could not open input/output files.\n";
        return 1;
    }

    string line_in, line_out;
    getline(fin_inputs, line_in);  // skip header
    getline(fin_outputs, line_out);

    fout << "num_parcels,time_microseconds,correct\n";

    int case_num = 0;
    while (getline(fin_inputs, line_in) && getline(fin_outputs, line_out)) {
        case_num++;

        // Input string and expected output string
        string input_str = line_in;
        string expected_str = line_out;

        // Remove quotes if CSV wrapped the field in quotes
        if (input_str.size() > 1 && input_str.front() == '"' && input_str.back() == '"')
            input_str = input_str.substr(1, input_str.size() - 2);
        if (expected_str.size() > 1 && expected_str.front() == '"' && expected_str.back() == '"')
            expected_str = expected_str.substr(1, expected_str.size() - 2);

        // Parse
        vector<vector<int>> parcels = parse_parcels(input_str);
        vector<vector<int>> expected = parse_parcels(expected_str);

        // Run and time
        auto start = chrono::high_resolution_clock::now();
        vector<vector<int>> got = question_one(parcels);
        auto end = chrono::high_resolution_clock::now();

        long long elapsed =
            chrono::duration_cast<chrono::microseconds>(end - start).count();

        // Compare
        bool correct = (got == expected);

        // Write result
        fout << parcels.size() << "," << elapsed << "," << (correct ? "1" : "0") << "\n";

        cout << "Completed test no. " << case_num << endl;
    }

    cout << "Finished testing. Results in q1_timings.csv\n";
    return 0;
}
