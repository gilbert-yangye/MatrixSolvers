#pragma once
#include <iostream>
#include <vector>

template <class T>
class Matrix {
public:
    //constructor when we want to preallocate ourselves
    Matrix(int rows, int cols, bool preallocate);
    //constructor when we already have allocated memory outside
    Matrix(int rows, int cols, T* values_ptr);
    //destructor
    virtual ~Matrix();

    //print out the values in our matrix
    void printValues();
    virtual void printMatrix();

    //perform some operations with our matrix
    virtual void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);
    virtual void matVecMult(T* input, T* output);
    virtual void vecVecMult(T* vec1, T* vec2, T& result);
    
    // several solvers
    virtual void Gauss_Seidel(T* mat_b, T* output);
    virtual void Jacobi_solve(T* mat_b, T* output);
    virtual void BiCGSTAB(T* mat_b, T* output);
    virtual void GMRES(T* mat_b, T* output);
    virtual void Cholesky(T* mat_b, T* output);
    void LU_solve(T* rhs, T* result);

    // pre-opearate on matrix
    virtual void sort_mat(double* rhs);



    //explicitly using the c++11 nullptr here
    T* values = nullptr;
    int rows = -1;
    int cols = -1;

    //the subclass will know the protected contents
protected:
    bool preallocated = false;
    virtual bool if_finish(T* mat_b, T* output);
    virtual void cal_norm(T* vec, int k, T& normed);
    virtual void arnoldi(T** Q, T** H, int k);
    virtual void apply_givens_rotation(T** H, T* cs, T* sn, int k);
    virtual void givens_rotation(T& h_k, T& h_k_1, T& cs_k, T& sn_k);
    virtual void getCfactor(T** M, T** t, int p, int q, int n);
    virtual T det(T** M, int n);
    virtual void adj(T** M, T** adj, int n);
    virtual void inverse(T** M, T** inv, int n);
    void LU_decomposition(T* L, T* R);
    virtual void backward_substitution(T* rhs, T* result);
    virtual void forward_substitution(T* rhs, T* result);
    virtual void find_unique(std::vector<bool> check_list, std::vector<int>& unique_list);
private:
    int size_of_values = -1;
    int it_max = 2000;


};
