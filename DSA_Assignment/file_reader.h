#ifndef FILE_READER_H
#define FILE_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <iterator>
#include <chrono>
using namespace std;

void question1_reader(const string& filename, vector<vector<int>>&data) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;
    bool headerSkipped = false;

    while (getline(file, line)) {
        if (!headerSkipped) {  
            // skip the header (first line)
            headerSkipped = true;
            continue;
        }

        if (line.empty()) continue;  // skip empty lines

        stringstream ss(line);
        string item;
        vector<int> row;

        while (getline(ss, item, ',')) {
            if (!item.empty()) {  // avoid trailing commas
                row.push_back(stoi(item));
            }
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    file.close();

}





void question2_reader(
    const string &file_path,
    vector<int> &preorder,
    vector<int> &inorder,
    vector<vector<int>> &leaf_parcels,
    vector<vector<int>> &queries
) {
    ifstream infile(file_path);
    if (!infile.is_open()) {
        cerr << "Error: Cannot open file " << file_path << endl;
        return;
    }

    string line;
    string section;

    while (getline(infile, line)) {
        if (line.empty()) continue;  // skip blank lines
        stringstream ss(line);

        // Detect section headers
        if (line.find("preorder:") == 0) {
            section = "preorder";
            line = line.substr(9); // remove "preorder: "
        } else if (line.find("inorder:") == 0) {
            section = "inorder";
            line = line.substr(8);
        } else if (line.find("parcel_on_leafs:") == 0) {
            section = "parcel";
            continue;
        } else if (line.find("queries:") == 0) {
            section = "queries";
            continue;
        }

        if (section == "preorder") {
            ss.clear(); ss.str(line);
            int num;
            while (ss >> num) preorder.push_back(num);
        } 
        else if (section == "inorder") {
            ss.clear(); ss.str(line);
            int num;
            while (ss >> num) inorder.push_back(num);
        } 
        else if (section == "parcel") {
            vector<int> parcels;
            int num;
            while (ss >> num) parcels.push_back(num);
            if (!parcels.empty()) leaf_parcels.push_back(parcels);
        } 
        else if (section == "queries") {
            vector<int> q;
            int num;
            while (ss >> num) q.push_back(num);
            if (!q.empty()) queries.push_back(q);
        }
    }
    infile.close();
}




void question3_reader(const string& file_path, vector<vector<int>>& edges, vector<int>& metro_cities) {
    ifstream infile(file_path);
    if (!infile.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    string line;
    bool reading_edges = false, reading_metro = false;

    while (getline(infile, line)) {
        if (line.empty()) continue;  // skip empty lines

        if (line.find("edges:") != string::npos) {
            reading_edges = true;
            reading_metro = false;
            continue;
        }
        if (line.find("metro_cities:") != string::npos) {
            reading_edges = false;
            reading_metro = true;
            continue;
        }

        stringstream ss(line);

        if (reading_edges) {
            vector<int> edge;
            int val;
            while (ss >> val) edge.push_back(val);
            if (!edge.empty()) edges.push_back(edge);
        }
        else if (reading_metro) {
            int city;
            while (ss >> city) metro_cities.push_back(city);
        }
    }
    infile.close();
}




#endif 