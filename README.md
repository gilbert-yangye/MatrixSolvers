# ACSE-5: Advanced Programming Assignment 2019/2020

## Introduction 
Implement different algorithms to solve the linear system **Ax=b**, where **A** is a positive definite matrix. We use **Jacobi, Gauss-Seidel, LU factorisation, Cholesky, GMRES, BiCGSTAB (Biconjugate gradient stabilized method)** to solve the linear system for either dense or sparse matrix and compare the effiency and scope of these methods.

## User instructions
In **main.cpp**, we use main() method to load matrix **A** and right vector **b**, and use methods "Solver4Dense" and "Solver4Sparse" to get the solution **x** to our linear system. And we have two classes **Matrix** and its subclass **CSRMatrix**.

The solver for the above 6 methods of the dense matrix and these methods are defined in **Matrix.hpp** and **Matrix.cpp**.
```
void Solver4Dense(Matrix<T>* dense_mat, T* rhs_b, T* result)
{
    dense_mat->sort_mat(rhs_b);
    dense_mat->Jacobi_solve(rhs_b, result);
    dense_mat->gauss_seidel(rhs_b, result);
    dense_mat->BiCGSTAB(rhs_b, result);
    dense_mat->GMRES(rhs_b, result);
    dense_mat->LU_solve(rhs_b, result);
    dense_mat->Cholesky(rhs_b, result);
}
```

The solver for  the above 6 methods of the sparse matrix and these methods are defined in **CSRMatrix.hpp** and **CSRMatrix.cpp**.
```
void Solver4Sparse(CSRMatrix<T>* sparse_mat, T* rhs_b, T* result)
{
    sparse_mat->sort_mat(rhs_b);
    sparse_mat->Jacobi_solve(rhs_b, result);
    sparse_mat->Gauss_seidel(rhs_b, result);
    sparse_mat->BiCGSTAB(rhs_b, result);
    sparse_mat->GMRES(rhs_b, result);
    sparse_mat->LU_solve(rhs_b, result);
    sparse_mat->Cholesky(rhs_b, result);
}
```
And in the main() method, we can easily call the above two function to get all 6 methods for dense matrix and sparse matrix respectively.
```
Solver4Dense(dense_mat, rhs_b, result);
Solver4Sparse(sparse_mat, rhs_b, result);
```
Other than the above public functions in **Matrix** and **CSRMatrix**, there are some **protected and private** functions which serve for these 6 methods.

In Matrix.h and Matrix.cpp:

The following functions are used to get the inverse of a matrix and serve for GMRES.
```
virtual void getCfactor(T** M, T** t, int p, int q, int n);
virtual T DET(T** M, int n);
virtual void ADJ(T** M, T** adj, int n);
virtual void INV(T** M, T** inv, int n);
virtual void cal_norm(T* vec,int k, T& normed);
virtual void arnoldi(T** Q, T** H, int k);
virtual void apply_givens_rotation(T** H, T* cs, T* sn, int k);
virtual void givens_rotation(T& h_k, T& h_k_1, T& cs_k, T& sn_k);
```
The following functions are used by LU methods and Cholesky methods.
```
void LU_decomposition(T* L, T* R);
virtual void backward_substitution(T* rhs, T* result);
virtual void forward_substitution(T* rhs, T* result);
```
In CSRMatrix.h and CSRMatrix.cpp:

The following are used to get the (i,j) value of a sparse matrix and store and update it in a new sparse matrix.
```
T get_A_value(int i, int j);
void update_A_value(int i, int j, T value);
void set_L_value(int i, int j, T value);
T get_L_value(int i, int j);
```
These are the functions needed by LU method and Cholesky method.
```
void LU_solve(T* rhs, T* result);
void LU_backward_substitution(T* rhs, T* result);
void LU_forward_substitution(T* rhs, T* result);
void Cholesky_backward_substitution(T* rhs, T* result);
void Cholesky_forward_substitution(T* rhs, T* result);
```
## Result
### Comparison and Discussion:

For these implemented functions, we tested their running time by passing in matrices with different properties (symmetric and diagonal-dominant) and size:

**<p align="center">Table 1: Run time (in ms) of Symmetric Diagonal-dominant Matrix</p >**

|       | Jacobi  |  Gauss-Seidel   | BiCGSTAB  | GMRES  | LU  |  Cholesky  | 
|  :----: | :----:    |  :----:     | :----:  |   :----:  | :----:  |:----:  |
| 10*10 | 0.21 |  0.08   | 0.03  | 16.091  | 0.006  | 0.005 |
| 100*100 | 19.29<br>(max iteration reached,<br>Invalid)  | 0.66   | 0.62    | 2.65  | 1.48  | 0.67  | 
| 500*500 | 464.12<br>(max iteration reached,<br>Invalid)  |   16.65   | 9.30   | 6.11  | 158.56  | 66.65  | 
| 1000*1000 | 1882.87<br>(max iteration reached,<br>Invalid)  |   63.91   | 38.33    | 19.37  | 1305.10  | 524.08  | 

The comparison between different solvers implicate the following points:
1. Though Gauss-Seidel and Jacobi are mostly similar, the convergence rate is much faster and the stability is better in Gauss-Seidel; 
2. Symmetricity does not influence the performance of Jacobi, Gauss-Seidel, and BiCGSTAB significantly, while all these methods do not converge and are possible to blow up when non-diagonal-dominant matrices are passed in;
3. LU decomposition with pivoting performs stably for all test conditions;
4. In contrast with all the other 5 algorithms we implemented, the GMRES algorithm is not sensitive to the size of matrices;

**<p align="center">Table 2: Run time (ms) of Non-symmetric Diagonal-dominant Matrix</p >**
    
|       | Jacobi  |  Gauss-Seidel   | BiCGSTAB  | GMRES  | LU  | 
|  :----: | :----:    |  :----:     | :----:  |   :----:  | :----:  |
| 10*10 | 0.34  |  0.21   | 0.08  | 39.61  | 0.01  | 
| 100*100 | 38.01<br>(max iteration reached)  |  1.62   | 0.93  | 5.29  | 2.86  | 
| 500*500 | 750.51<br>(max iteration reached)  |  25.52   | 16.10  | 10.87  | 347.52  | 
| 1000*1000 | 1954.75<br>(max iteration reached)  |  65.723   | 30.541  | 18.312  | 1299.1  | 
</br>

**<p align="center">Table 3: Run time (ms) of Symmetric Non-Diagonal-dominant Matrix</p >**

|       | Jacobi  |  Gauss-Seidel   | BiCGSTAB  | GMRES  | LU  |  Cholesky  | 
|  :----: | :----:    |  :----:     | :----:  |   :----:  | :----:  |:----:  |
| 10*10 | 0.214<br>(max iteration reached,<br>Invalid)  |  0.242<br>(max iteration reached,<br>Invalid)   | 0.082  | 14854.8  | 0.008  | 0.005<br>(Invalid) |
| 100*100 | 17.415<br>(max iteration reached,<br>Invalid)  | 12.925<br>(max iteration reached,<br>Invalid)   | 26.308<br>(max iteration reached,<br>Invalid)    | Time out  | 1.427  | 0.66<br>(Invalid)  | 
| 500*500 | 468.836<br>(max iteration reached,<br>Invalid)  |   325.054<br>(max iteration reached,<br>Invalid)   | 621.342<br>(max iteration reached,<br>Invalid)    | Time out  | 153.004  | 68.382<br>(Invalid)  | 
| 1000*1000 | 1285.99<br>(max iteration reached,<br>Invalid)  |   1806.2<br>(max iteration reached,<br>Invalid)   | 2598.51<br>(max iteration reached,<br>Invalid)    | Time out  | 1284.91  | 549.868<br>(Invalid)  | 
</br>

**<p align="center">Table 4: Run time (ms) of Non-Symmetric Non-Diagonal-dominant Matrix</p >**

|       | Jacobi  |  Gauss-Seidel   | BiCGSTAB  | GMRES  | LU  | 
|  :----: | :----:    |  :----:     | :----:  |   :----:  | :----:  |
| 10*10 | 0.212<br>(max iteration reached,<br>Invalid)  |  0.227<br>(max iteration reached,<br>Invalid)   | 0.05  | 14520.4  | 0.005  | 
| 100*100 | 14.578<br>(max iteration reached,<br>Invalid)  |  17.63<br>(max iteration reached,<br>Invalid)   | 26.687<br>(max iteration reached,<br>Invalid)  | Time out  | 1.766  | 
| 500*500 | 446.51<br>(max iteration reached,<br>Invalid)  |  321.944<br>(max iteration reached,<br>Invalid)   | 642.484<br>(max iteration reached,<br>Invalid)  | Time out  | 170.243  | 
| 1000*1000 | 1873.36<br>(max iteration reached,<br>Invalid)  |  1270.7<br>(max iteration reached,<br>Invalid)   | 2613.4<br>(max iteration reached,<br>Invalid)  | Time out  | 1301.75  | 
</br>

**<p align="center">Table 5: Comparison of run time (ms) between Sparse and Dense Matrix</p >**

|       | Jacobi  |  Gauss-Seidel   | BiCGSTAB  | GMRES  | LU  |  Cholesky  | 
|  :----: | :----:    |  :----:     | :----:  |   :----:  | :----:  |:----:  |
| 100*100<br>dense | 19.29<br>(max iteration reached,<br>Invalid)  |  0.66   | 0.62  | 2.65  | 1.48  | 0.67 |
| 100*100<br>sparse | 0.373  | 0.325   | 0.399    | 29.95  | 10.24  | 6.601  | 

Compared with dense matrix, time consumed in solving a sparse matrix with CSR format is less in the first three methods. In contrast, the latter three methods spend more time. This is caused by our complicate indexing method in CSR format, which caused over-complicity in the GMRES, LU, and Cholesky. This is consistent with the presumption that though CSR format is efficient in data storage, it may not be efficient in the context of data manipulation and treatment. It’s also presumed that when it comes to the larger size of matrices, the data manipulation will become one of the advantages of CSR format.

## Reference

- ACSE5-Assignment.pdf

