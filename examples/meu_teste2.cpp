#ifdef AMGX_DYNAMIC_LOADING
#include "amgx_capi.h"
#else
#include "amgx_c.h"
#endif

#include <stdio.h>
#include <amgx_c.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

void gerarMatrizS(double stepSize)
{

    double a = -1, b = 1; // intervalo

    double alfa = 2.0;
    double beta = 2.0;

    int numeroPontos = ((b - a) / stepSize) + 1;
    cout << "numeroPontos: " << numeroPontos << endl;

    /*
        numeroPontos = 9 são 7 valores + 2 valores de contorno
        na resolução do sistema não uso o contorno então faço -2
    */

    int n = numeroPontos - 2;

    // double stepSize = (b - a) / (numeroPontos - 1); // valor de cada incremento, 9 pontos mas 8 intervalos
    //cout << "n: " << n << endl;
    //cout << "STEP SIZE: " << stepSize << endl;

    double h = stepSize * stepSize;
    //cout << "h: " << h << endl;

    int nnz = 3 * n - 2; // elementos não-nulos - n+(n-1)+(n-1) = 3n-2
    //cout << "NNZ: " << nnz << endl;
    int countI = 1, countF = 3;

    ofstream arquivoMtx;

    arquivoMtx.open("../examples/matrix3.mtx");
    arquivoMtx <<fixed<< scientific << setprecision(15);
    arquivoMtx << "%%MatrixMarket matrix coordinate real general symmetric" << endl;
    arquivoMtx << "%%AMGX 1 1 sorted rhs" << endl;
    arquivoMtx << n << " " << n << " " << nnz << endl;
    arquivoMtx << 1 << " " << 1 << " " << 2.0 << endl;
    arquivoMtx << 1 << " " << 2 << " " << -1.0 << endl;

    for (int i = 2; i <= n - 1; i++)
    {
        for (int j = countI; j <= countF; j++)
        {
            if (i == j)
            {
                arquivoMtx << i << " " << j << " " << 2.0 << endl;
            }
            else
            {
                arquivoMtx << i << " " << j << " " << -1.0 << endl;
            }
        }
        countI++;
        countF++;
    }

    arquivoMtx << n << " " << n - 1 << " " << -1.0 << endl;
    arquivoMtx << n << " " << n << " " << 2.0 << endl;

    double vd = -h * 2.0 + alfa;
    //cout << "RHS: " << vd << endl;

    arquivoMtx << vd << endl;
    for (int i = 0; i < n - 2; i++)
    {
        arquivoMtx << -2.0 * h << endl;
    }
    arquivoMtx << vd;

    arquivoMtx.close();
}

void calcular(const char **argv, double stepSize)
{
    //cout << "Arquivo matrix: " << argv[2] << endl;
    //cout << "Config: " << argv[4] << endl;
   
    //cout << "STEP SIZE: " << stepSize << endl;
    gerarMatrizS(stepSize);

    AMGX_initialize();
    
    // All of the objects are initialized using a default Resources.
    AMGX_matrix_handle matrix;
    AMGX_vector_handle rhs;
    AMGX_vector_handle soln;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_config_handle config;

    // arquivo de configuração passado como parametro: ../src/configs/FGMRES_AGGREGATION_JACOBI.json
    AMGX_config_create_from_file(&config, argv[4]);
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
    AMGX_read_system(matrix, rhs, soln, argv[2]);

    //AMGX_write_system(matrix, rhs, soln, "./output.system.mtx");
    AMGX_solver_setup(solver, matrix);

    AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln);

    AMGX_vector_download(soln, data);

    int sol_size, sol_bsize;
    AMGX_vector_get_size(soln, &sol_size, &sol_bsize);

    ofstream plotSol;

    plotSol.open("plotSol.csv");
    plotSol << 2 <<endl;
    for (int i = 0; i < sol_size; ++i)
    {
        //printf("%f \n",data[i]);
        plotSol << data[i]<<endl;
    }
    plotSol << 2 <<endl;
    plotSol.close();
    
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(soln);
    AMGX_vector_destroy(rhs);
    AMGX_matrix_destroy(matrix);
    AMGX_resources_destroy(rsrc);

    AMGX_SAFE_CALL(AMGX_config_destroy(config));
    AMGX_SAFE_CALL(AMGX_finalize());

    delete [] data;
}

int main(int argc, const char **argv)
{
    /*
     g++ gerarMatriz.cpp  -o gerarMtx && ./gerarMtx 9 && rm gerarMtx && make -j16 all && examples/meu_teste2 -m ../examples/matrix3.mtx -c ../src/configs/FGMRES_AGGREGATION.json
     g++ ../examples/gerarMatriz.cpp  -o gerarMtx && ./gerarMtx 1000 && rm gerarMtx && make -j16 all && examples/meu_teste2 -m ../examples/matrix3.mtx -c ../src/configs/FGMRES_AGGREGATION.json
     make -j16 all && examples/meu_teste2 -m ../examples/matrix3.mtx -c ../src/configs/FGMRES_AGGREGATION.json -n 9
     
     make -j16 all && examples/meu_teste2 -m ../examples/matrix3.mtx -c ../src/configs/FGMRES_AGGREGATION.json -s 0.025
     make -j16 all && examples/meu_teste2 -m ../examples/matrix3.mtx -c ../src/configs/FGMRES_AGGREGATION.json -n 9
     */

    //cout << "Arquivo matrix: " << argv[2] << endl;
    //cout << "Config: " << argv[4] << endl;
    int numeroPontos;
    double stepSize;
    istringstream ss(argv[6]);
    if (std::string(argv[5]) == "-n")
    {
        ss >> numeroPontos;
    }
    else if (std::string(argv[5]) == "-s")
    {
        ss >> stepSize;

    }

    calcular(argv, stepSize);
    stepSize = stepSize / 2.0;

    //system("python ../examples/plotSol.py");

    return 0;
}