#include <iostream>
#include <algorithm>
#include "Matrix.hpp"

using namespace std;

template<class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) :rows(rows), cols(cols), size_of_values(rows* cols), preallocated(preallocate) {
    if (this->preallocated) {
        this->values = new T[size_of_values];
    }
}

template<class T>
Matrix<T>::Matrix(int rows, int cols, T* values_ptr) :rows(rows), cols(cols), size_of_values(rows* cols), values(values_ptr) {}

//destructor
template<class T>
Matrix<T>::~Matrix() {
    if (this->preallocated) {
        delete[] this->values;
    }
}

template<class T>
void Matrix<T>::printValues() {
    cout << "printing values" << endl;
    for (int i = 0; i < this->size_of_values; i++) {
        cout << this->values[i] << " ";
    }
    cout << endl << endl;
}

template<class T>
void Matrix<T>::printMatrix() {
    //cout << "printing matrix" << endl;
    for (int j = 0; j < this->rows; ++j) {
        cout << endl;
        for (int i = 0; i < this->cols; i++) {
            cout << this->values[i + j * this->cols] << " ";
        }
    }
    cout << endl << endl;
}


template<class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_right, Matrix<T>& output) {
    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != output.rows || this->cols != output.cols)
        {
            std::cerr << "Input dimensions for matrices don't match" << std::endl;
            return;
        }
    }
    // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->rows * mat_right.cols];
    }

    //set values to 0 before hand
    for (int i = 0; i < output.size_of_values; i++) {
        output.values[i] = 0;
    }

    //do matrix multiplication
    for (int i = 0; i < this->rows; i++) {
        for (int k = 0; k < this->cols; k++) {
            for (int j = 0; j < mat_right.cols; j++) {
                output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
            }
        }
    }
}

//do multiplication between matrix and vector
template<class T>
void Matrix<T>::matVecMult(T* input, T* output) {
    for (int i = 0; i < this->rows; i++) {
        output[i] = 0;
        for (int j = 0; j < this->cols; j++) {
            output[i] += this->values[i * this->cols + j] * input[j];
        }
    }
}

//do multiplication between vector and vector
template<class T>
void Matrix<T>::vecVecMult(T* vec1, T* vec2, T& result) {
    //do multiplication
    result = 0;
    for (int i = 0; i < this->rows; i++) result += vec1[i] * vec2[i];
}

//check if current A*x is close to b
template <class T>
bool Matrix<T>::if_finish(T* mat_b, T* output)
{
    T tol = 0.000001;
    T* cal_out = new T[this->rows];
    T res = 0;

    // calculate the resident
    for (int i = 0; i < this->rows; i++)
    {
        cal_out[i] = 0;
        for (int j = 0; j < this->cols; j++)
        {
            cal_out[i] += this->values[i * this->cols + j] * output[j];
        }
        res += abs(mat_b[i] - cal_out[i]);
    }

    // if the resident is small enough, stop the interation
    delete[] cal_out;
    if ((res / this->rows) < tol)
        return true;
    else
        return false;
}

template <class T>
void Matrix<T>::Jacobi_solve(T* mat_b, T* output) {

    T* temp = new T[this->rows];

    // set output to zero
    for (int i = 0; i < this->rows; i++) {
        output[i] = 0;
    }

    int k;

    //Iterate up to  user - specified max number of iterations, we will 'break'
    //this loop if we hit our convergence tolerance
    for (k = 0; k < this->it_max; k++) {
        for (int i = 0; i < this->rows; i++) {
            temp[i] = mat_b[i] / this->values[i * this->cols + i];
            for (int j = 0; j < this->cols; j++) {
                if (i == j) {
                    continue;
                }
                temp[i] = temp[i] - ((this->values[i * this->cols + j] / this->values[i * this->cols + i]) * output[j]);
            }
        }
        //after this iteration we update the values of x
        for (int i = 0; i < this->rows; i++) {
            output[i] = temp[i];
        }
        if (this->if_finish(mat_b, output)) {
            break;
        }
    }

    delete[] temp;

    //print the result
    if (k < this->it_max) {
        cout << "!!!!!!!! we got Jacobi solution for the dense matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for Jacobi solution for the dense matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, Jacobi solution cannot solve it!!" << endl << endl;
    }

}

template <class T>
void Matrix<T>::Gauss_Seidel(T* mat_b, T* output)
{
    // Set output to zero before hand
    for (int i = 0; i < this->rows; i++)  output[i] = 0;

    int i, j, k;

    //Iterate up to  user - specified max number of iterations, we will 'break'
    //this loop if we hit our convergence tolerance
    for (k = 0; k < this->it_max; k++)
    {
        for (i = 0; i < this->rows; i++)
        {
            T sum_ax = 0;
            for (j = 0; j < this->cols; j++)
                if (j != i)
                    sum_ax += this->values[i * this->cols + j] * output[j];
            //we will iterate and update at the same time
            output[i] = 1 / this->values[i * this->cols + i] * (mat_b[i] - sum_ax);
        }
        if (this->if_finish(mat_b, output)) {
            break;
        }
    }

    // print the result
    if (k < this->it_max) {
        cout << "!!!!!!!! we got Gauss-Seidel solution for the dense matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for Gauss-Seidel solution for the dense matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, Gauss-Seidel solution cannot solve it!!" << endl << endl;
    }
}

template <class T>
void Matrix<T>::BiCGSTAB(T* mat_b, T* output)
{
    //define related vector
    // the residual
    T* r_0_head = new T[this->rows];
    T* r_i_1 = new T[this->rows];
    T* r_i = new T[this->rows];
    // v=Ap
    T* v_i_1 = new T[this->rows];
    T* v_i = new T[this->rows];
    //the search directions
    T* p_i_1 = new T[this->rows];
    T* p_i = new T[this->rows];
    T* x_i_1 = new T[this->rows];
    T* x_i = new T[this->rows];
    T* h = new T[this->rows];
    T* s = new T[this->rows];
    T* t = new T[this->rows];
    // Initialize vectors to zero before hand
    for (int i = 0; i < this->rows; i++) {
        x_i_1[i] = 0;
        v_i_1[i] = 0;
        p_i_1[i] = 0;
    }
    //Ax_0=A*x_0
    T* Ax_0 = new T[this->rows];
    this->matVecMult(x_i_1, Ax_0);
    // initialize r_0_head=b-A*x0
    // initialize r_i-1=b-A*x0
    for (int i = 0; i < this->rows; i++) {
        r_0_head[i] = mat_b[i] - Ax_0[i];
        r_i_1[i] = mat_b[i] - Ax_0[i];
    }
    //define and initializa our parameters
    T rho_i_1 = 1., alpha = 1., omega_i_1 = 1.;
    T rho_i, beta, r_0_head_v_i, t_s, t_t, omega_i;
    // k is the interation number
    int k;
    for (k = 0; k < this->it_max; k++) {
        //rho_i=r_0_head * r_i_1
        this->vecVecMult(r_0_head, r_i_1, rho_i);
        beta = (rho_i / rho_i_1) * (alpha / omega_i_1);
        for (int i = 0; i < this->rows; i++) p_i[i] = r_i_1[i] + beta * (p_i_1[i] - omega_i_1 * v_i_1[i]);
        //v_i=A*p_i
        this->matVecMult(p_i, v_i);
        //alpha = rho_i / (r_0_head * v_i)
        this->vecVecMult(r_0_head, v_i, r_0_head_v_i);
        alpha = rho_i / r_0_head_v_i;
        for (int i = 0; i < this->rows; i++) h[i] = x_i_1[i] + alpha * p_i[i];
        // if finished, stop the iteration
        if (this->if_finish(mat_b, h)) {
            for (int i = 0; i < this->rows; i++) output[i] = h[i];
            break;
        }
        // s = r_i_1 - alpha * v
        for (int i = 0; i < this->rows; i++) s[i] = r_i_1[i] - alpha * v_i[i];
        // t = A*s
        this->matVecMult(s, t);
        // omega = (t * s)/(t * t)
        this->vecVecMult(t, s, t_s);
        this->vecVecMult(t, t, t_t);
        omega_i = t_s / t_t;
        //x = h + omega*s
        for (int i = 0; i < this->rows; i++) x_i[i] = h[i] + omega_i * s[i];
        if (this->if_finish(mat_b, x_i)) {
            for (int i = 0; i < this->rows; i++) output[i] = x_i[i];
            break;
        }
        //r_i = s - omega * t
        for (int i = 0; i < this->rows; i++) r_i[i] = s[i] - omega_i * t[i];
        // update values
        rho_i_1 = rho_i;
        omega_i_1 = omega_i;
        for (int i = 0; i < this->rows; i++) {
            r_i_1[i] = r_i[i];
            v_i_1[i] = v_i[i];
            p_i_1[i] = p_i[i];
            x_i_1[i] = x_i[i];
        }
    }
    // print results
    if (k < this->it_max) {
        cout << "!!!!!!!! we got BiCGSTAB solution for the dense matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for BiCGSTAB solution for the dense matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, BiCGSTAB solution cannot solve it!!" << endl << endl;
    }
    delete[] Ax_0;
    delete[] r_0_head;
    delete[] r_i_1;
    delete[] r_i;
    delete[] v_i_1;
    delete[] v_i;
    delete[] p_i_1;
    delete[] p_i;
    delete[] x_i;
    delete[] x_i_1;
    delete[] h;
    delete[] s;
    delete[] t;
}

//calculate the norm of a vector
template <class T>
void Matrix<T>::cal_norm(T* vec, int k, T& normed) {
    normed = 0;
    for (int i = 0; i < k; i++) {
        normed += pow(vec[i], 2);
    }
    normed = pow(normed, 0.5);
}

//used for GMRES to updata Q and H
template <class T>
void Matrix<T>::arnoldi(T** Q, T** H, int k) {

    this->matVecMult(Q[k], Q[k + 1]);

    for (int i = 0; i < k + 1; i++) {
        this->vecVecMult(Q[k + 1], Q[i], H[k][i]);

        for (int j = 0; j < this->rows; j++) Q[k + 1][j] -= H[k][i] * Q[i][j];
    }
    this->cal_norm(Q[k + 1], this->rows, H[k][k + 1]);
    for (int j = 0; j < this->rows; j++) Q[k + 1][j] = Q[k + 1][j] / H[k][k + 1];

}

//used for GMRES to determine the k-th value of sn and cs
template <class T>
void Matrix<T>::givens_rotation(T& h_k, T& h_k_1, T& cs_k, T& sn_k) {
    T t;
    if (h_k == 0) {
        cs_k = 0;
        sn_k = 1;
    }
    else {
        t = sqrt(h_k * h_k + h_k_1 * h_k_1);
        cs_k = abs(h_k) / t;
        sn_k = cs_k * h_k_1 / h_k;
    }
}

//used for GMRES to update H
template <class T>
void Matrix<T>::apply_givens_rotation(T** H, T* cs, T* sn, int k) {
    T temp;
    for (int i = 0; i < k; i++) {

        temp = cs[i] * H[k][i] + sn[i] * H[k][i + 1];
        H[k][i + 1] = -sn[i] * H[k][i] + cs[i] * H[k][i + 1];
        H[k][i] = temp;

    }
    //update the next sin cos values for rotation
    this->givens_rotation(H[k][k], H[k][k + 1], cs[k], sn[k]);

    //eliminate H[k][k + 1]
    H[k][k] = cs[k] * H[k][k] + sn[k] * H[k][k + 1];
    H[k][k + 1] = 0;
}

// used for get the inverse, get cofactors
// p and q are the given row and column
template <class T>
void Matrix<T>::getCfactor(T** M, T** t, int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        //copy only those elements which are not in given row r and column c
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                t[i][j++] = M[row][col];
                //if row is filled go to next row and reset column index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

//to find determinant
template <class T>
T Matrix<T>::det(T** M, int n) {
    T D = 0;
    if (n == 1) { return M[0][0]; }
    //store cofactors
    T** t = new T * [n];
    for (int i = 0; i < n; i++) {
        t[i] = new T[n];
    }
    //store sign multiplier
    T s = 1;
    for (int f = 0; f < n; f++) {
        //for getting cofactor of M[0][f]
        this->getCfactor(M, t, 0, f, n);
        D += s * M[0][f] * this->det(t, n - 1);
        s = -s;
    }
    for (int i = 0; i < n; i++) {
        delete[] t[i];
    }
    delete[] t;
    return D;
}

//to find adjoint matrix
template <class T>
void Matrix<T>::adj(T** M, T** adj, int n) {
    if (n == 1) {
        adj[0][0] = 1; return;
    }
    int s = 1;
    ////store cofactors
    T** t = new T * [n];
    for (int i = 0; i < n; i++) {
        t[i] = new T[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            //to get cofactors
            this->getCfactor(M, t, i, j, n);
            //add sign multiplier
            s = ((i + j) % 2 == 0) ? 1 : -1;
            //to get the transpose of the cofactor matrix
            adj[j][i] = (s) * (this->det(t, n - 1));
        }
    }
    for (int i = 0; i < n; i++) {
        delete[] t[i];
    }
    delete[] t;
}

//get the inverse of a matrix
template <class T>
void Matrix<T>::inverse(T** M, T** inv, int n) {
    T det = this->det(M, n);
    if (det == 0) {
        cout << "can't find its inverse";
    }
    T** adj = new T * [n];
    for (int i = 0; i < n; i++) {
        adj[i] = new T[n];
    }
    this->adj(M, adj, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) inv[i][j] = adj[i][j] / float(det);
    for (int i = 0; i < n; i++) {
        delete[] adj[i];
    }
    delete[] adj;
}


template <class T>
void Matrix<T>::GMRES(T* mat_b, T* output) {
    int k;
    //store sin values for rotation
    T* sn = new T[this->it_max];
    //store cos values for rotation
    T* cs = new T[this->it_max];
    T* e1 = new T[this->rows];
    //initial residual vector for initial output
    T* r = new T[this->rows];
    //residual norm vector through iterations
    T* beta = new T[this->it_max + 1];
    T error;
    T tol = 0.000001;

    // initialize
    for (int i = 0; i < this->rows; i++) {
        output[i] = 0;
        e1[i] = 0;
    }
    e1[0] = 1;
    //initialize sn and cs
    for (int i = 0; i < this->it_max; i++) {
        sn[i] = 0;
        cs[i] = 0;
    }

    T* A_x = new T[this->rows];
    this->matVecMult(output, A_x);
    for (int i = 0; i < this->rows; i++) r[i] = mat_b[i] - A_x[i];
    T b_norm;
    this->cal_norm(mat_b, this->rows, b_norm);
    T r_norm;
    this->cal_norm(r, this->rows, r_norm);

    //Q is a it_max*rows matrix
    T** Q;
    Q = (T**) new T * [this->it_max];
    for (int i = 0; i < this->it_max; i++) {
        Q[i] = new T[this->rows];
    }
    // initialize Q to 0
    for (int i = 0; i < this->it_max; i++)
        for (int j = 0; j < this->rows; j++)
            Q[i][j] = 0;

    //H is a it_max*it_max matrix
    T** H;
    H = (T**) new T * [this->it_max];
    for (int i = 0; i < this->it_max; i++) {
        H[i] = new T[this->it_max];
    }

    // initialize H to 0
    for (int i = 0; i < this->it_max; i++)
        for (int j = 0; j < this->it_max; j++)
            H[i][j] = 0;

    //initialize to Q and beta
    for (int i = 0; i < this->rows; i++) {
        Q[0][i] = r[i] / r_norm;
        beta[i] = r_norm * e1[i];
    }

    for (k = 0; k < this->it_max; k++) {

        //generate Q[k+1] and H[k+1]
        this->arnoldi(Q, H, k);

        //updata H[k+1]
        this->apply_givens_rotation(H, cs, sn, k);

        //update the residual vector
        beta[k + 1] = -sn[k] * beta[k];
        beta[k] = cs[k] * beta[k];
        error = abs(beta[k + 1]) / b_norm;

        if (error <= tol) break;
    }

    delete[] r;
    delete[] A_x;
    delete[] sn;
    delete[] cs;
    delete[] e1;

    // calculate the result
    T** HH = new T * [k + 1];
    for (int i = 0; i < k + 1; i++) {
        HH[i] = new T[k + 1];
    }

    for (int i = 0; i < k + 1; i++) {
        for (int j = 0; j < k + 1; j++) {
            HH[j][i] = H[i][j];
        }
    }

    //store the inverse of HH
    T** inv = new T * [k + 1];
    for (int i = 0; i < k + 1; i++) {
        inv[i] = new T[k + 1];
    }


    //store y y=H^(-1)*beta
    T* y = new T[k + 1];
    this->inverse(HH, inv, k + 1);
    for (int i = 0; i < k + 1; i++) {
        y[i] = 0;
        for (int j = 0; j < k + 1; j++) {
            y[i] += inv[i][j] * beta[j];
        }
    }

    //get the result
    T* Q_y = new T[this->rows];
    for (int i = 0; i < this->rows; i++) {
        Q_y[i] = 0;
        for (int j = 0; j < k + 1; j++) {
            Q_y[i] += Q[j][i] * y[j];
        }
        output[i] = output[i] + Q_y[i];
    }
    if (k < this->it_max) {
        cout << "!!!!!!!! we got GMRES solution for the dense matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for GMRES solution for the dense matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, GMRES solution cannot solve it!!" << endl << endl;
    }

    for (int i = 0; i < this->it_max; i++) delete[] Q[i];
    delete[] Q;
    delete[] beta;
    for (int i = 0; i < this->it_max; i++) delete[] H[i];
    delete[] H;
    for (int i = 0; i < k + 1; i++) delete[] HH[i];
    delete[] HH;
    for (int i = 0; i < k + 1; i++) delete[] inv[i];
    delete[] inv;
    delete[] y;
    delete[] Q_y;
}

template <class T>
void Matrix<T>::LU_solve(T* rhs, T* result)
{
    //LU_decomposition(T & lmat, T & umat)
    auto* lmat = new T[this->size_of_values];
    auto* umat = new T[this->size_of_values];
    // L*y=b, for_result means y
    auto* for_result = new T[this->rows];
    this->LU_decomposition(lmat, umat);
    //define matrix to store lower triangular matrix and upper triangular matrix
    auto* lmatrix = new Matrix<T>(this->rows, this->cols, lmat);
    auto* umatrix = new Matrix<T>(this->rows, this->cols, umat);
    // L*y=b, U*x=y, solve these two system in the following step
    lmatrix->forward_substitution(rhs, for_result);
    umatrix->backward_substitution(for_result, result);
    // print result
    cout << "!!!!!!!! we got LU solution for the dense matrix(only the first ten x values):" << endl;
    for (int i = 0; i < 10; i++)
    {
        cout << result[i] << " ";
    }
    cout << endl << endl;
    delete[] lmat;
    delete[] umat;
    delete[] for_result;
    delete lmatrix;
    delete umatrix;
}

// this function is used to decomposite A=LU
template <class T>
void Matrix<T>::LU_decomposition(T* lmat, T* umat)
{
    // judge if row number equals to col number
    if (this->cols != this->rows) std::cout << "rows != cols, cannot be decomposed with LU method";
    // initialize lower triangular matrix and upper triangular matrix
    for (int i = 0; i < this->size_of_values; i++)
    {
        lmat[i] = 0;
        umat[i] = this->values[i];
    }
    for (int k = 0; k < this->cols - 1; k++)
    {
        for (int i = k + 1; i < this->cols; i++)
        {
            // store U[i,k] / U[k,k] as temp_s, which will be used to update U
            double temp_s = (umat[i * this->cols + k] / umat[k * this->cols + k]);
            // update U[i,j]=U[i,j]- temp_s*U[k,j]
            for (int j = k; j < this->cols; j++)
            {
                umat[i * this->cols + j] = umat[i * this->cols + j] - temp_s * umat[k * this->cols + j];
            }
            // set L[i,k] to be temp_s
            lmat[i * this->cols + k] = temp_s;
        }
    }
    for (int i = 0; i < this->rows; i++)
    {
        //L[i,i] add 1
        lmat[i * this->cols + i] += 1;
    }
}

template <class T>
void Matrix<T>::forward_substitution(T* rhs, T* result)
{
    // initialize result as 0
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }
    // L* y = b
    // the y calculated at this step will be used in the next step
    double s;
    for (int k = 0; k < this->rows; k++) {
        s = 0.0;
        for (int j = 0; j < k; j++)
        {
            s += this->values[k * this->cols + j] * result[j];
        }
        result[k] = (rhs[k] - s) / this->values[k * this->cols + k];
    }
}
template <class T>
void Matrix<T>::backward_substitution(T* rhs, T* result)
{
    // initialize result as 0
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }
    // U*x=y
    // the x calculated at this step will be used in the next step
    for (int k = this->rows - 1; k > -1; k--) {
        double s(0.);
        for (int j = k + 1; j < this->rows; j++)
        {
            s += this->values[k * this->cols + j] * result[j];
        }
        result[k] = (rhs[k] - s) / this->values[k * this->cols + k];
    }
}

template <class T>
void Matrix<T>::Cholesky(T* rhs, T* result) {

    //determine whether this matrix is symmetric or not
    int x = 0;
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            if (this->values[i * this->cols + j] != this->values[j * this->cols + i]) {
                x = 1;
            }
        }
    }
    if (x == 1) {
        cout << "!!!!!!!! Input matrix is not symmetric, we can't use Cholesky method!!" << endl << endl;
    }

    else {
        //initialize a lower triangular matrix
        auto* L = new T[this->size_of_values];
        //initialize the transpose of above matrix
        auto* L_T = new T[this->size_of_values];
        T L1, L2;

        for (int i = 0; i < this->size_of_values; i++) L[i] = 0;

        //decompose original matrix and generate L and L_T
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < i + 1; j++) {
                if (i == j) {
                    L1 = 0;
                    //set L1=sum(L[j,k]^2) when k<j
                    for (int k = 0; k < j; k++) L1 += L[j * this->rows + k] * L[j * this->rows + k];
                    //set L[j,j]=(A[j,j]-L1)^2
                    L[j * this->rows + j] = pow((this->values[j * this->rows + j] - L1), 0.5);
                    //set L_T[j,j]=L[j,j]
                    L_T[j * this->rows + j] = L[j * this->rows + j];
                }
                else {
                    L2 = 0;
                    //set L2=sum(L[i,k]*L[j,k]) when k<j
                    for (int k = 0; k < j; k++) L2 += L[i * this->rows + k] * L[j * this->rows + k];
                    //set L[i,j]=1/L[j,j]*(A[i,j]-L2)
                    L[i * this->rows + j] = 1 / L[j * this->rows + j] * (this->values[i * this->rows + j] - L2);
                    //set L_T[j,i]=L[i,j]
                    L_T[j * this->rows + i] = L[i * this->rows + j];
                }
            }
        }


        auto* for_result = new T[this->rows];

        auto* lmatrix = new Matrix<T>(this->rows, this->cols, L);
        auto* umatrix = new Matrix<T>(this->rows, this->cols, L_T);
        //apply forward substitution to lower triangular matrix
        lmatrix->forward_substitution(rhs, for_result);
        //apply backward substitution to transpose of L
        umatrix->backward_substitution(for_result, result);

        cout << "!!!!!!!! we got Cholesky solution for the dense matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << result[i] << " ";
        }
        cout << endl << endl;

        delete[] for_result;
        delete lmatrix;
        delete umatrix;
        delete[] L;
        delete[] L_T;
    }

}

template <class T>
void Matrix<T>::find_unique(std::vector<bool> check_list, std::vector<int>& unique_list)
{
    //loop over all the columns in °∞check_list°± to check how many non-zeros are inside;
    //if only one element in that column, mark down the corresponding row index in °∞unique_list°±
    int row_index = 0;
    int count = 0;
    for (int i = 0; i < this->cols; i++)
    {
        if (check_list[i] == false) continue;
        for (int j = 0; j < this->rows; j++)
        {
            if (this->values[i + j * this->cols] != 0)
            {
                row_index = j;
                count++;
            }
        }
        if (count >= 2)
        {
            unique_list[i] = -1;
        }
        else
        {
            unique_list[i] = row_index;
        }
        count = 0;
    }
}
template <class T>
void Matrix<T>::sort_mat(double* rhs) {
    // sort the matrix and corresponding rhs in order to eliminate the non-zeros on diagonal;
    // and also, do partial pivoting to make the pivot element more significant;
    auto* temp_mat = new Matrix<double>(this->rows, this->cols, true);
    auto* temp_rhs = new double[this->rows];
    //use °∞check_list°± to track columns to be sorted
    std::vector <bool> check_list(this->cols, true);

    //when check_list is not empty (all false), continue the while loop
    while (!(std::none_of(check_list.begin(), check_list.end(), [](bool v) { return v; })))
    {
        //use °∞unique_list°± to find columns that have only one non-zero
        std::vector <int> unique_list(this->cols, -1);
        //update unique_list with unique function
        this->find_unique(check_list, unique_list);
        //if column j has a unique entry on row i (equals to °∞unique_list[j]°∞)
        //then in temp_mat, set row j equals to (row i in original matrix)
        //so that in temp_mat, the entry on [i,j] is the unique one;
        //set this column j as false in the while loop to be excluded
        if (!(std::all_of(unique_list.begin(), unique_list.end(), [](int v) { return v == -1; })))
        {

            for (int j = 0; j < this->cols; j++)
            {
                if (unique_list[j] != -1)
                {
                    for (int col = 0; col < this->cols; col++)
                    {
                        temp_mat->values[j * this->cols + col] = this->values[unique_list[j] * this->cols + col];
                        this->values[unique_list[j] * this->cols + col] = 0;
                    }
                    temp_rhs[j] = rhs[unique_list[j]];
                    check_list[j] = false;
                }
            }
        }

        //if no unique element found, fill the 1st available column with the row
        //where value that will be the diagonal element is maximum in that row
        //and remove the column from check_list;
        else {
            for (int j = 0; j < this->cols; j++)
            {
                if (check_list[j])
                {
                    int index_row(-1);
                    double max_value(-std::numeric_limits<int>::max());
                    double diff;
                    //now we fill temp_mat row j with value in row x; let°Øs find out x
                    //row x should have the most significance element at position j
                    // significance means Axj/sum(Aij) is maximum
                    for (int row = 0; row < this->rows; row++)
                    {
                        if (this->values[row * this->cols + j] != 0)
                        {
                            diff = 0;
                            for (int aa = 0; aa < this->cols; aa++)
                            {
                                diff += abs(this->values[row * this->cols + aa]);
                            }

                            diff = abs(this->values[row * this->rows + j]) / (diff - abs(this->values[row * this->rows + j]));
                            if (diff >= max_value)
                            {
                                index_row = row;
                                max_value = diff;
                            }
                        }

                    }
                    // now index_row takes the index of row x
                    // fill and exclude
                    if (index_row != -1) {
                        for (int kk = 0; kk < this->cols; kk++)
                        {
                            temp_mat->values[j * this->cols + kk] = this->values[index_row * this->cols + kk];
                            this->values[index_row * this->cols + kk] = 0;
                        }
                        temp_rhs[j] = rhs[index_row];
                        check_list[j] = false;
                    }

                    if (index_row == -1)
                    {
                        std::cerr << std::endl << "bad matrix input" << std::endl;
                    }
                    break;
                }
            }
        }

    }
    //fill the values to original matrix
    for (int i = 0; i < this->size_of_values; i++)
    {
        this->values[i] = temp_mat->values[i];
    }

    for (int i = 0; i < this->rows; i++)
    {
        rhs[i] = temp_rhs[i];
    }

    delete temp_mat;
    delete[] temp_rhs;

}
