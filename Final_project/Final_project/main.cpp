#include <iostream>
#include <ctime>
#include <math.h>
#include <time.h>
#include "Matrix.hpp"
#include "Matrix.cpp"
#include "CSRMatrix.hpp"
#include "CSRMatrix.cpp"

using namespace std;
template <class T>
void generate_dense(int dimension, T* result, T* rhs, bool dominant, bool symmetric)
{
    for (int i = 0; i < dimension; i++)
    {
        rhs[i] = rand() % 20 + 1;
    }
    if (dominant)
    {
        T* temp = new T[dimension];
        for (int i = 0; i < dimension; i++) { temp[i] = 0; }
        if (symmetric)
        {
            for (int i = 0; i < dimension; i++)
            {
                for (int j = i; j < dimension; j++)
                {
                    result[i * dimension + j] = rand() % 10 + 1;
                    //                    std::cout<<"tem_Val is "<<temp_val<<std::endl;
                }
                for (int j = 0; j < i; j++)
                {
                    result[i * dimension + j] = result[j * dimension + i];
                }
            };
        }
        else
        {
            for (int i = 0; i < dimension * dimension; i++)
            {
                result[i] = rand() % 10 + 1;
            };
        }
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++)
            {
                temp[i] += result[i * dimension + j];
            }
            result[i * dimension + i] = temp[i];
        }
        delete[] temp;
    }
    else
    {
        T temp_val;
        if (symmetric)
        {
            for (int i = 0; i < dimension; i++)
            {
                for (int j = i; j < dimension; j++)
                {
                    temp_val = rand() % 10 + 1;
                    result[i * dimension + j] = temp_val;
                }
                for (int j = 0; j < i; j++)
                {
                    result[i * dimension + j] = result[j * dimension + i];
                }
            };
        }
        else
        {
            for (int i = 0; i < dimension * dimension; i++)
            {
                temp_val = rand() % 10 + 1;
                result[i] = temp_val;
            };
        }
    }
}
template <class T>
void Solver4Dense(int rows, int cols, bool dominant, bool symmetric) {
    //define our dense matrix
    auto* dense_mat = new Matrix<T>(rows, cols, true);//auto here automatically detect that the class to be created is a Matrix class
    //define our rhs vector
    T* rhs_b = new T[rows];
    generate_dense<T>(rows, dense_mat->values, rhs_b, dominant, symmetric);
    //print the original matrix
    //cout << "The original marix is " << endl;
    //dense_mat->printMatrix();
    //sort matrix to make its diagonal values larger
    dense_mat->sort_mat(rhs_b);
    //print the matrix after sorting
    //cout << "The sorted marix is " << endl;
    //dense_mat->printMatrix();
    //define our result
    T* result = new T[rows];
    cout << "We only output the first ten x results of the solution because of limited space" << endl << endl;
    clock_t start;
    clock_t end;
    start = clock();
    dense_mat->Gauss_Seidel(rhs_b, result);
    end = clock();
    cout << "Gauss_Seidel(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << endl << endl;
    start = clock();
    dense_mat->Jacobi_solve(rhs_b, result);
    end = clock();
    std::cout << "Jacobi(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    dense_mat->BiCGSTAB(rhs_b, result);
    end = clock();
    std::cout << "BiCGSTAB(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    dense_mat->GMRES(rhs_b, result);
    end = clock();
    std::cout << "GMRES(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    dense_mat->LU_solve(rhs_b, result);
    end = clock();
    std::cout << "LU(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    dense_mat->Cholesky(rhs_b, result);
    end = clock();
    std::cout << "Cholesky(for dense) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    delete[] rhs_b;
    delete dense_mat;
    delete[] result;
}
template <class T>
void Solver4Sparse(T non_zero_values[], T rhs_b_data[], int col_pos[], int row_pos[], int nnzs, int rows, int cols, T* result) {
    //define our sparse matrix
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    // lets fill our matrix
    for (int i = 0; i < nnzs; i++) {
        sparse_mat->values[i] = non_zero_values[i];
        sparse_mat->col_index[i] = col_pos[i];
    }
    for (int i = 0; i < rows + 1; i++) sparse_mat->row_position[i] = row_pos[i];
    //define our rhs vector and fill it
    double* rhs_b = new double[rows];
    for (int i = 0; i < rows; i++)   rhs_b[i] = rhs_b_data[i];
    //print the original matrix
    //cout << "The original dense matrix is " << endl;
    //sparse_mat->printMatrix();
    //sort matrix to make its diagonal values larger
    //sparse_mat->sort_mat(rhs_b);
    //print the matrix after sorting
    //cout << "The sorted marix is " << endl;
    //sparse_mat->printMatrix();
    cout << "We only output the first ten x results of the solution because of limited space" << endl << endl;
    clock_t start;
    clock_t end;
    start = clock();
    sparse_mat->Gauss_seidel(rhs_b, result);
    end = clock();
    std::cout << "Gauss_Seidel(for sparse)takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    sparse_mat->Jacobi_solve(rhs_b, result);
    end = clock();
    std::cout << "Jacobi(for sparse) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    sparse_mat->BiCGSTAB(rhs_b, result);
    end = clock();
    std::cout << "BiCGSTAB(for sparse) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    sparse_mat->GMRES(rhs_b, result);
    end = clock();
    std::cout << "GMRES(for sparse) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    sparse_mat->Cholesky(rhs_b, result);
    end = clock();
    std::cout << "Cholesky(for sparse) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    start = clock();
    sparse_mat->LU_solve(rhs_b, result);
    end = clock();
    std::cout << "LU(for sparse) takes " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    delete[] rhs_b;
    delete sparse_mat;
}
int main()
{
    bool dominant, symmetric;
    int rows, cols;
    ///////////////////////////// First Dense Matrix Solver////////////////////////////////////
    //set type of matrix
    dominant = true;
    symmetric = true;
    rows = 100;
    cols = 100;
    //call solvers to solve the linear system
    Solver4Dense<double>(rows, cols, dominant, symmetric);
    ///////////////////////First Sparse Matrix Solver(10 * 10)/////////////////////////////////////////
    int nnzs = 164;
    rows = 100;
    cols = 100;
    double non_zero_values[164] = { 3001, 2680, 2991,   56,   57, 2920, 3018,   58, 2929,   57, 2947,
          56, 3150, 3056,   57, 3163,   56, 2925,   57, 2960,   56, 3070,
        2797, 2953,   56, 3145, 3021,   58, 2892, 2996,   56, 3159,   56,
        2943,   56, 2925, 2929,   56,   56, 2716, 3166,   56, 2872, 3093,
        2997, 2948,   56, 2802, 3164,   56, 3130, 2939,   56,   57, 2940,
        2962, 3015,   57, 3089,   56, 2893,   56, 2799, 3002, 2913, 2929,
        2823,   56,   57, 3121,   58,   57, 3158,   56, 2948, 3165,   56,
          57,   57, 3022, 2940,   58, 2876,   56, 2956,   57, 2939, 2880,
        2966, 3085, 2924, 2942, 3055,   57, 2854,   56, 2788, 2811, 2921,
          57,   57, 3004, 3055, 3042, 2917, 2817,   56, 2912,   56,   57,
          56, 3044, 2927, 2745, 3123, 2910,   57, 2886, 2794, 3058,   57,
        2861,   57, 3102, 2819,   58,   57, 3001, 2950,   56, 3163,   56,
        2933, 3119,   56, 3186,   56, 2906, 2906, 3140, 3087,   58,   56,
        3098,   58,   56,   56, 2845,   56, 3053,   57,   58, 2878, 2768,
          57, 3067, 2821,   57, 3048,   56,   56, 3164, 2826, 2930 };
    int col_pos[164] = { 0,  1,  2, 67, 77,  3,  4, 89,  5, 68,  6, 68,  7,  8, 92,  9,
        19, 10, 43, 11, 97, 12, 13, 14, 29, 15, 16, 90, 17, 18,  9, 19,
        38, 20, 84, 21, 22, 37, 90, 23, 24, 85, 25, 26, 27, 28, 14, 29,
        30, 89, 31, 32, 45, 76, 33, 34, 35, 73, 36, 22, 37, 19, 38, 39,
        40, 41, 42, 81, 10, 43, 92, 94, 44, 32, 45, 46, 50, 51, 96, 47,
        48, 79, 49, 46, 50, 46, 51, 52, 53, 54, 55, 56, 57, 62, 58, 90,
        59, 60, 61, 79, 57, 62, 63, 64, 65, 66,  2, 67, 91,  5,  6, 68,
        69, 70, 71, 72, 35, 73, 74, 75, 32, 76,  2, 77, 78, 48, 61, 79,
        80, 42, 81, 97, 82, 83, 20, 84, 24, 85, 86, 87, 88,  4, 30, 89,
        16, 22, 58, 90, 67, 91,  8, 43, 92, 93, 43, 94, 95, 46, 96, 11,
        81, 97, 98, 99 };
    int row_pos[101] = { 0,   1,   2,   5,   6,   8,  10,  12,  13,  15,  17,  19,  21,
        22,  23,  25,  26,  28,  29,  30,  33,  35,  36,  39,  40,  42,
        43,  44,  45,  46,  48,  50,  51,  54,  55,  56,  58,  59,  61,
        63,  64,  65,  66,  68,  72,  73,  75,  79,  80,  82,  83,  85,
        87,  88,  89,  90,  91,  92,  94,  96,  97,  98, 100, 102, 103,
       104, 105, 106, 109, 112, 113, 114, 115, 116, 118, 119, 120, 122,
       124, 125, 128, 129, 132, 133, 134, 136, 138, 139, 140, 141, 144,
       148, 150, 153, 154, 156, 157, 159, 162, 163, 164 };
    double rhs_b_data2[100] = { 13, 16, 15, 14, 29, 16, 15, 12, 27, 12, 27,  4, 15, 27, 29, 20,
         5,  3, 15, 10, 12,  6, 13, 24,  9,  3, 16,  9,  8, 29, 28, 15,
        14,  3,  4, 26, 28,  5, 29, 16, 12, 22, 16,  5, 14, 15, 28,  2,
        15, 12,  1, 20,  4, 26,  8, 11,  8, 11, 17, 14, 19,  1,  7,  3,
        26, 16, 29, 14,  2, 23, 21,  1, 19, 26, 24,  7, 19, 24,  5, 11,
         5, 19, 28,  1, 22, 22, 28, 22,  3, 13, 17, 14, 12,  8, 28, 23,
        27, 29,  9, 18 };
    double* result2 = new double[rows];
    Solver4Sparse(non_zero_values, rhs_b_data2, col_pos, row_pos, nnzs, rows, cols, result2);
    delete[] result2;
    
    
    std::cout<<std::endl<<"Let's test the sort_mat with CSRMatrix"<<std::endl<<std::endl;
    double mat_A_sort[100] =
    {140,  1,  2,  2,  5,  5,  7,  9,  7,  9,
        1, 0,  6,  1,  19,  7,  9,  9,  16,  3,
        7,  1, 150,  3,  2,  3,  2,  8,  8,  2,
        2,  3,  6, 190,  2,  3,  1,  3,  3,  3,
        1,  100,  3,  7, 0,  4,  2,  3,  1,  1,
        2,  2,  4,  1,  9, 170,  7,  3,  9,  8,
        5,  16,  1,  7,  3,  9, 120,  7,  9,  3,
        6,  0,  6,  5,  4,  4,  5, 160,  6,  8,
        7,  1,  5,  8,  2,  6,  3,  9, 170,  7,
        2,  2,  1,  7,  5,  1,  4,  3,  5, 140};
    rows = 10;
    cols = 10;
    auto *dense_mat = new Matrix<double>(rows, cols, true);
    for (int i = 0; i < rows * cols; i++)
    {
        dense_mat->values[i] = mat_A_sort[i];
    }
    dense_mat -> printMatrix();
    double rhs_b_sort[10] = {10., 20., 30., 40.,50,60,70,80,90,100};
    double *rhs_b_so = new double[10];
    for (int i = 0; i < rows; i++)
    {
        rhs_b_so[i] = rhs_b_sort[i];
    }
    dense_mat->sort_mat(rhs_b_so);
    std::cout<<std::endl<<"After pivoting"<<std::endl;
    dense_mat->printMatrix();
    
    std::cout<<std::endl<<"Let's test the sort_mat with CSRMatrix"<<std::endl<<std::endl;
    nnzs = 8;
    auto *m = new CSRMatrix<double>(4, 4, nnzs, true);
    int row_pos_sort[5] = { 0,2,3,5,8 };
    int col_pos_sort[8] = { 1,3,0,2,3,1,2,3 };
    double tempp_sort[8] ={14.,2.,10.,16.,3.,2.,3.,11.};
    for (int i = 0; i<8;i++)
    {
        m->col_index[i] =col_pos_sort[i];
        m->values[i] = tempp_sort[i];
    }
    for (int i = 0; i<8;i++)
    {
        m->row_position[i] = row_pos_sort[i];
    }
    double *rhs_b = new double[4];
    double rhs_b_data[4] = {10., 20., 30., 40.};
    for (int i = 0; i < 4; i++)
    {
        rhs_b[i] = rhs_b_data[i];
    }
    m->printMatrix();
    m->sort_mat(rhs_b);
    std::cout<<std::endl<<"After pivoting"<<std::endl;
    m->printMatrix();
    delete [] rhs_b;
    delete dense_mat;
    delete m;
    delete[] rhs_b_so;
   
    
    return 0;
    
}
