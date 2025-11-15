#include "sparse_matrix.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <fstream> // Needed for file I/O in parallel
#include <mpi.h>

//Preprocessor macro for ceiling division
#define CEIL_DIV(a, b) (((a) + (b) - 1) / (b))

//Preprocessor macro to add commas for CSV
#define COMMA ","

//Preprocessor macro for separating different runs in output
#define SEPARATOR "==========================="

using namespace std;

//File paths for matrix and vector 
const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";

// Output file stream for writing results (append)
ofstream output_file("/scratch/suhaskamath/MPI/individual_computation_times.csv", ios::app); 

int read_metadata(int& n, int& nnz, vector<int>& full_row_ptr, vector<double>& x)
{

    ifstream file(matrix_file, ios::binary);

    if (!file.is_open())
    {
        cerr << "Root failed to open matrix file." << endl;
        return -1;
    }

    file.read(reinterpret_cast<char*>(&n), sizeof(int));
    file.read(reinterpret_cast<char*>(&nnz), sizeof(int));

    full_row_ptr.resize(n + 1);
    file.read(reinterpret_cast<char*>(full_row_ptr.data()), (n + 1) * sizeof(int));
    file.close();

    if (!load_vector(vector_file, x)) {
        cerr << "Root failed to load vector." << endl;
        return -1;
    }

    return 0;
}

vector<int> get_row_offset_of_each_process(int size, const CSRMatrix& A)
{
    // Compute nnz per row
    vector<int> row_nnz(A.n, 0);
    for (int i = 0; i < A.n; ++i)
        row_nnz[i] = A.row_ptr[i + 1] - A.row_ptr[i];

    long long total_nnz = 0;
    for (int v : row_nnz) total_nnz += v;
    long long target = (size > 0) ? ((total_nnz + size - 1) / size) : total_nnz; // ceil div

    vector<int> starts;
    starts.reserve(size + 1);
    starts.push_back(0);
    long long acc = 0;
    int cur_row = 0;
    for (int r = 1; r < size; ++r) {
        long long goal = (long long)r * target;
        while (cur_row < A.n && acc + row_nnz[cur_row] <= goal) {
            acc += row_nnz[cur_row];
            ++cur_row;
        }
        starts.push_back(cur_row);
    }
    starts.push_back(A.n);

    // Clamp monotonicity and bounds
    for (int i = 1; i < (int)starts.size(); ++i) {
        if (starts[i] < starts[i - 1]) starts[i] = starts[i - 1];
        if (starts[i] > A.n) starts[i] = A.n;
    }
    return starts;
}

int mpi(int argc, char** argv)
{
    double total_time = MPI_Wtime();

    // Initialize MPI
    MPI_Init(&argc, &argv);

    //Rank denotes the process ID
    //Size denotes the total number of processes
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Local CSR chunk for this rank
    CSRMatrix A;
    vector<double> x;

    //n denotes matrix size
    //nnz denotes number of non-zero elements
    int n = 0, nnz = 0;

    // Full row pointer for the entire matrix
    // Contains the offset of each row in the global values/col_idx arrays
    vector<int> full_row_ptr;

    //read matrix metadata and vector on root process
    //No one except for rank 0 will read metadata
    //If the rank of a process is not 0, the first condition will be false
    //Due to short-circuit evaluation, the second condition will not be evaluated
    if (rank == 0 && (read_metadata(n, nnz, full_row_ptr, x) != 0 || n <= 0))
        // If read_metadata returns non-zero (error), abort MPI
        // If matrix size is non-positive, abort MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
    
    //Broadcast matrix metadata and vector to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Resize vectors on non-root processes
    //Vectors on the root process are already correctly sized
    if (rank != 0)
    {
        //When the full_row_ptr and x vectors are resized, they will have
        //the correct size to receive the broadcasted data
        full_row_ptr.resize(n + 1);
        x.resize(n);
    }

    //The entire vector is required by each process for matrix multiplication
    //So the root process reads the vector and broadcasts it to all processes
    //We do not do this for the matrix, because only a section of the matrix is needed by each process
    MPI_Bcast(full_row_ptr.data(), n + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Build a metadata CSR to compute row partition locally on every rank
    CSRMatrix A_meta;
    A_meta.n = n;
    A_meta.nnz = nnz;
    A_meta.row_ptr = full_row_ptr;

    // Compute row offsets for all ranks locally (no point-to-point sends needed)
    vector<int> row_offset_of_each_process = get_row_offset_of_each_process(size, A_meta);
    vector<int> row_offset_for_this_process(2, 0);
    row_offset_for_this_process[0] = row_offset_of_each_process[rank];
    row_offset_for_this_process[1] = row_offset_of_each_process[rank + 1];


    int number_of_rows_i_have_to_process = row_offset_for_this_process[1] - row_offset_for_this_process[0];
    int global_nnz_offset = full_row_ptr[row_offset_for_this_process[0]];
    int number_of_nnz_i_have_to_process = full_row_ptr[row_offset_for_this_process[1]] - global_nnz_offset;

    //Read matrix data for the rows assigned to this process
    long long col_idx_file_start = (long long)sizeof(int) * 2 + (long long)sizeof(int) * (n + 1);
    long long values_file_start = col_idx_file_start + (long long)sizeof(int) * nnz;

    ifstream file(matrix_file, ios::binary);
    if (!file.is_open()) {
        cerr << "Rank " << rank << " failed to open matrix file." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //Read row_ptr for the rows assigned to this process
    A.col_idx.resize(number_of_nnz_i_have_to_process);

    //Read col_idx and values for the rows assigned to this process
    long long my_col_seek = col_idx_file_start + (long long)sizeof(int) * (long long)global_nnz_offset;
    
    //Set the file read position to the start of the col_idx data for this process
    file.seekg(my_col_seek);

    //Read the col_idx data for this process
    file.read(reinterpret_cast<char*>(A.col_idx.data()), number_of_nnz_i_have_to_process * sizeof(int));

    //Read values for the rows assigned to this process
    A.values.resize(number_of_nnz_i_have_to_process);

    //Set the file read position to the start of the values data for this process
    long long my_val_seek = values_file_start + (long long)sizeof(double) * (long long)global_nnz_offset;

    //Set the file read position to the start of the values data for this process
    file.seekg(my_val_seek);

    //Read the values data for this process
    file.read(reinterpret_cast<char*>(A.values.data()), number_of_nnz_i_have_to_process * sizeof(double));
    
    //Close the matrix file after reading
    file.close();

    //Set local matrix dimensions
    A.n = number_of_rows_i_have_to_process;
    A.nnz = number_of_nnz_i_have_to_process;
    A.row_ptr.resize(number_of_rows_i_have_to_process + 1);

    //Popular A.row_ptr for the local matrix
    for (int i = 0; i <= number_of_rows_i_have_to_process; ++i)
        A.row_ptr[i] = full_row_ptr[row_offset_for_this_process[0] + i] - global_nnz_offset;


    //Initialize the result vector for this process
    vector<double> local_result(number_of_rows_i_have_to_process, 0.0);

    //Initialise clock for computation timing
    double time = MPI_Wtime();
    //Perform local computation for the rows assigned to this process
    if (number_of_rows_i_have_to_process > 0)
    {
        for (int i = 0; i < number_of_rows_i_have_to_process; i++)
        {
            double sum = 0.0;
            int row_begin = A.row_ptr[i];
            int row_end = A.row_ptr[i + 1];

            //Compute the dot product for this row with the vector x
            for (int j = row_begin; j < row_end; j++)
                sum += A.values[j] * x[A.col_idx[j]];

            local_result[i] = sum;
        }
    }

    time = MPI_Wtime() - time;

    vector<int> recv_counts, recv_displs;
    int my_rows = number_of_rows_i_have_to_process;
    if (rank == 0) recv_counts.resize(size);
    MPI_Gather(&my_rows, 1, MPI_INT, (rank == 0 ? recv_counts.data() : nullptr), 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<double> y; 
    if (rank == 0) {
        recv_displs.resize(size, 0);
        for (int r = 1; r < size; ++r) recv_displs[r] = recv_displs[r-1] + recv_counts[r-1];
        y.resize(n, 0.0);
    }
    
    MPI_Gatherv(local_result.data(), my_rows, MPI_DOUBLE,
                (rank == 0 ? y.data() : nullptr),
                (rank == 0 ? recv_counts.data() : nullptr),
                (rank == 0 ? recv_displs.data() : nullptr), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    //Finalize MPI
    MPI_Finalize();

    if (rank == 0)
    {
        output_file << SEPARATOR << endl;
        cout << "Size: " << size << endl;
        cout << "Total execution time: " << MPI_Wtime() - total_time << " seconds." << endl;
    }

    output_file << size << COMMA << rank << COMMA << time << endl;

    return 0;
}

int serial()
{
    double total_time = MPI_Wtime();
    CSRMatrix A;
    vector<double> x, y;

    if (!load_matrix(matrix_file, A)) {
        cerr << "Failed to load matrix." << endl;
        return -1;
    }

    if (!load_vector(vector_file, x)) {
        cerr << "Failed to load vector." << endl;
        return -1;
    }

    y.resize(A.n, 0.0);
    
    double comp_time = MPI_Wtime();
    for (int i = 0; i < A.n; i++) {
        y[i] = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++)
            y[i] += A.values[j] * x[A.col_idx[j]];
    }

    // Print CSV format: no_of_processes,rank_of_process,computation_time_of_process
    // Store computation time into CSV
    output_file << 1 << COMMA << 0 << COMMA << (MPI_Wtime() - comp_time) << endl << SEPARATOR << endl;
    

    //Print total execution time to console
    cout << "Size: " << 1 << endl;
    cout << "Total execution time: " << MPI_Wtime() - total_time << " seconds." << endl;

    return 0;
}

int main(int argc, char** argv)
{
    //Use sequential mode by default
    int use_mpi = 0;

    if (argc >= 2)
        use_mpi = atoi(argv[1]) == 1 ? 1 : 0;

    if (use_mpi)
        mpi(argc, argv);
    else
        serial();


    return 0;
}