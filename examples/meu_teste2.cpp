#ifdef AMGX_DYNAMIC_LOADING
#include "amgx_capi.h"
#else
#include "amgx_c.h"
#endif

#include <string>
#include <stdio.h>
#include <amgx_c.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

template <typename T>
void printArray(T arr[], int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

void calcular(const char **argv)
{
    /*
        0 -> nome do programa
        1 -> -c
        2 -> arquivo config
    */

    // cout << "Config: " << argv[2] << endl;

    AMGX_initialize();

    // All of the objects are initialized using a default Resources.
    AMGX_matrix_handle matrix;
    AMGX_vector_handle rhs;
    AMGX_vector_handle soln;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_config_handle config;

    // arquivo de configuração passado como parametro: ../src/configs/FGMRES_AGGREGATION_JACOBI.json
    AMGX_config_create_from_file(&config, argv[2]);
    AMGX_resources_create_simple(&rsrc, config);

    AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);

    double *data = new double[100000000];

    /*
        1 d -> the first letter h or d specifies whether the matrix data (and subsequent linear solver algorithms) will run on the host or device
        2 D(or F) -> specifies the precision (double or float) of the Matrix data
        3 D The third D or F specifies the precision (double or float) of any Vector (including right-hand side or unknown vectors)
        4 I The last I specifies that 32-bit int types are used for all indices
    */
    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);

    // Next, data is uploaded from the application (or set, in the case of the solution vector which is initialized to all zeroes).
    //  If these are not specified than rhs=[1,...,1]^T and (initial guess) sol=[0,...,0]^T.
    // gerarMatrizS(matrix, rhs, soln, stepSize);
    // AMGX_read_system(matrix, rhs, soln, "../examples/matrix3.mtx");

     /*
    double values[] = {1, -2, -3, 1, 1, -4, -5, 1};
    int col_ind[] = {0, 1};
    int row_ptr[] = {0, 1, 2};
    double rhs_values[] = {1, 2, 3, 4};
    double soln_values[] = {0, 0, 0, 0};
    */

    int n = 0;
    int nnz = 0;
    double *values;
    int *col_indices;
    int *row_ptrs;
    double *rhs_values;
    double *soln_values;
    int block_size = 3;

    string problem = "problems/ve-mech";
    string path = problem + "/";
    std::string filename = path + "csr_matrix.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    std::string line;

    // Read ia: (number of values)
    std::getline(file, line);
    line = line.substr(3);
    n = std::stoi(line);
    n -= 1;
    cout << n << endl;

    // Allocate memory for the arrays
    row_ptrs = new int[n + 1];
    int aux_n = n;
    n = ceil(n / block_size);
    int *aux_row_ptrs = new int[n + 1];

    // Read ia: (integer values separated by spaces)
    std::getline(file, line);
    std::istringstream iss(line);
    for (int i = 0; i <= aux_n; i++)
    {
        iss >> row_ptrs[i];
        row_ptrs[i] -= 1;
    }

    int new_nnz = 0;
    aux_row_ptrs[0] = 0;

    for (int i = 0; i < n; i++)
    {
        new_nnz = 0;
        // if the next block_lines rows are not empty, count 1 to the nnz
        if (row_ptrs[(i + 1) * block_size] - row_ptrs[i * block_size] > 0)
        {
            new_nnz++;
        }
        aux_row_ptrs[i + 1] = aux_row_ptrs[i] + new_nnz;
    }

    row_ptrs = aux_row_ptrs;

    // Read ja: (integer values separated by spaces)
    std::getline(file, line);
    line = line.substr(3);
    nnz = std::stoi(line);
    cout << nnz << endl;
    col_indices = new int[nnz];
    values = new double[nnz];

    iss.clear();
    std::getline(file, line);
    iss.str(line);

    for (int i = 0; i < nnz; i++)
    {
        iss >> col_indices[i];
        col_indices[i] -= 1;
    }

    int aux_nnz = nnz;
    nnz = ceil(aux_nnz / (block_size * block_size));
    int *aux_col_indices = new int[nnz];

    for (int i = 0; i < nnz; i++)
    {
        aux_col_indices[i] = col_indices[i * block_size * block_size] / block_size;
    }

    col_indices = aux_col_indices;

    // Read values: (double values separated by spaces)
    iss.clear();
    std::getline(file, line);
    std::getline(file, line);
    iss.str(line);

    for (int i = 0; i < aux_nnz; i++)
    {
        iss >> values[i];
    }

    file.close();

    // Read rhs: (double values separated by newlines)
    rhs_values = new double[aux_n];
    soln_values = new double[aux_n];

    filename = path + "rhs.txt";
    std::ifstream rhs_file(filename);
    if (!rhs_file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    for (int i = 0; i < aux_n; i++)
    {
        rhs_file >> rhs_values[i];
    }

    rhs_file.close();

    // Soln is initialized to 0
    for (int i = 0; i < aux_n; i++)
    {
        soln_values[i] = 0.0;
    }

    cout << "n: " << n << endl;
    cout << "nnz: " << nnz << endl;


    AMGX_matrix_upload_all(matrix, n, nnz, block_size, block_size, row_ptrs, col_indices, values, NULL);
    AMGX_vector_upload(rhs, n, block_size, rhs_values);
    AMGX_vector_upload(soln, n, block_size, soln_values);

    AMGX_write_system(matrix, rhs, soln, "./output.system.mtx");

    AMGX_solver_setup(solver, matrix);

    AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln);

    AMGX_vector_download(soln, data);

    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(soln);
    AMGX_vector_destroy(rhs);
    AMGX_matrix_destroy(matrix);
    AMGX_resources_destroy(rsrc);

    AMGX_SAFE_CALL(AMGX_config_destroy(config));
    AMGX_SAFE_CALL(AMGX_finalize());

    delete[] data;
}

int main(int argc, const char **argv)
{
    // exemplo de chamada: examples/meu_teste2 -c ../src/configs/FGMRES_AGGREGATION.json -p cube-poisson
    // make && examples/meu_teste2 -c ../src/configs/GMRES_AMG_D2.json -s 0.002
    /*
        0 -> nome do programa
        1 -> -c
        2 -> arquivo config
    */

    // cout << "Config: " << argv[2] << endl;


    calcular(argv);
    return 0;
}

