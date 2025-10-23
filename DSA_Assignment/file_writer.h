#ifndef FILE_WRITER_H
#define FILE_WRITER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <string>


void question1_writer(const std::string& file_path, const std::vector<std::vector<int>>& output) {
    std::ofstream fout(file_path);  // open file for writing
    if (!fout.is_open()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return;
    }

    fout << "id,weight\n";  // header (optional)

    for (const auto& row : output) {
        for (size_t i = 0; i < row.size(); i++) {
            fout << row[i];
            if (i != row.size() - 1) fout << ","; // add commas between values
        }
        fout << "\n";  // new line for each row
    }

    fout.close();
    std::cout << "Output written to " << file_path << std::endl;
}



// Function to write a vector<int> to file
void question2_writer(const std::string &file_path, const std::vector<int> &vec) {
    std::ofstream fout(file_path);
    if (!fout.is_open()) {
        std::cerr << "Error opening file for writing: " << file_path << std::endl;
        return;
    }

    for (size_t i = 0; i < vec.size(); i++) {
        fout << vec[i];
        if (i != vec.size() - 1) fout << " "; // space separated
    }
    fout << "\n";
    fout.close();
}



// Function to write a single integer to file
void question3_writer(const std::string &file_path, int value) {
    std::ofstream fout(file_path);
    if (!fout.is_open()) {
        std::cerr << "Error opening file for writing: " << file_path << std::endl;
        return;
    }
    fout << value << "\n";
    fout.close();
}

#endif
