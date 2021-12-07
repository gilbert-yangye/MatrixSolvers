#pragma once
#include <iostream>
#include "Matrix.hpp"
#include <vector>

template <class T>
class CSRMatrix : public Matrix<T>
{
public:
    //constructor where we want to preallocate ourselves
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    //constructor where we already have allocated memory outside that is values_ptr
    CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
    
    ~CSRMatrix();
    //print the matrix
    virtual void printMatrix();
    //perform some matrix operations
    void matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);
    void matVecMult(T* input, T* output);
    //some solver
    void Gauss_seidel(T* mat_b, T* output);
    void Jacobi_solve(T* mat_b, T* output);
    void BiCGSTAB(T* mat_b, T* output);
    void GMRES(T* mat_b, T* output);
    void Cholesky(T* mat_b, T* output);
    void LU_solve(T* rhs, T* result);
    //pre-operate for matrix
    void sort_mat(double* rhs);
    // Variables
    int* row_position = nullptr;
    int* col_index = nullptr;
    //T* values = nullptr;
    //how many non-zero entries we have in the matrix
    int nnzs = -1;
protected:
    bool if_finish(T* mat_b, T* output);
    T get_A_value(int i, int j);
    void set_L_value(int i, int j, T value);
    T get_L_value(int i, int j);
    void Cholesky_backward_substitution(T* rhs, T* result);
    void Cholesky_forward_substitution(T* rhs, T* result);
    void update_A_value(int i, int j, T value);
    void LU_backward_substitution(T* rhs, T* result);
    void LU_forward_substitution(T* rhs, T* result);
    void find_unique_p(std::vector<bool> check_list, std::vector<int>& unique_list);
private:
    int it_max = 2000;
    std::vector<T> L_value;
    std::vector<T> L_col_index;
    std::vector<T> L_row_pos;
};
