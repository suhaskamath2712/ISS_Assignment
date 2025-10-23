#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <chrono>
#include <sstream>
#include "user_code.h"
#include "file_reader.h"
#include "file_writer.h"
using namespace std;
using namespace std::chrono;



/**
 * @brief Do not modify this code.
 *
 */
void question1(string input_file, string output_file)
{
    vector<vector<int>> parcels;
    // read from input file
    question1_reader(input_file, parcels);

    auto start = high_resolution_clock::now();

    vector<vector<int>> output =  question_one(parcels);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    auto computeDuration = duration.count();
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
    std::cout << "Total time taken: " << computeDuration << " microseconds" << endl;

    // save output to output_file
    question1_writer(output_file, output);

}

/**
 * @brief Do not modify this code.
 *
 */
void question2(string input_file, string output_file)
{

    vector<int> preorder, inorder;
    vector<vector<int>> leaf_parcels, queries;
    // read from input file
    question2_reader(input_file, preorder, inorder, leaf_parcels, queries);

    auto start = high_resolution_clock::now();

    vector<int> output =  question_two(preorder, inorder, leaf_parcels, queries);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    auto computeDuration = duration.count();
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
    std::cout << "Total time taken: " << computeDuration << " microseconds" << endl;
   
    // write output to file
    question2_writer(output_file, output);

}


/**
 * @brief Do not modify this code.
 *
 */

void question3(string input_file, string output_file) 
{

    vector<vector<int>> edges;
    vector<int> metro_cities;  
    // read from input file
    question3_reader(input_file, edges, metro_cities);


    auto start = high_resolution_clock::now();

    int output =  question_three(edges, metro_cities);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    auto computeDuration = duration.count();
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
    std::cout << "Total time taken: " << computeDuration << " microseconds" << endl;

    // write output ot file
    question3_writer(output_file, output);

}


//////////////////////////////////////////////////////////////////////////////////
// DO NOT MODIFY THIS SECTION
//////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Do not modify this code. FIXME add comments.
 *
 * @param argc
 * @param argv
 * @return int
 */

int main(int argc, char **argv)
{
    /*
    string question1_input_file = "input1.txt";
    string question1_output_file = "question1_output.txt";
    string question2_input_file = "input2.txt";
    string question2_output_file = "question2_output.txt";
    */

    string question3_input_file = "question3_edgecase.txt";
    string question3_output_file = "question3_edgecase_output.txt";
    
    /*
    string question1_input_file = argv[1];
    string question1_output_file = argv[2];
    string question2_input_file = argv[3];
    string question2_output_file = argv[4];
    string question3_input_file = argv[5];
    string question3_output_file = argv[6];
    */

    /*
    question1(question1_input_file, question1_output_file);

    question2(question2_input_file, question2_output_file);
    */
   
    question3(question3_input_file, question3_output_file);

    return 0;
}
