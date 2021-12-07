#include <iostream>
#include "CSRMatrix.hpp"
#include "Matrix.hpp"
#include <vector>
#include <algorithm>

using namespace std;
template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate) :Matrix<T>(rows, cols, false), nnzs(nnzs) {
    this->preallocated = preallocate;
    if (this->preallocated) {
        //values are the non-zero values
        this->values = new T[this->nnzs];
        //the next minus the previous equals the number of non-zero entries we have in this row
        this->row_position = new int[this->rows + 1];
        //the non-zero entries column index
        this->col_index = new int[this->nnzs];
    }
}

template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index) :Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index) {}

template<class T>
CSRMatrix<T>::~CSRMatrix() {
    if (this->preallocated) {
        //delete[] values;
        //we didn't need to do this step,because the super class has done this
        delete[] row_position;
        delete[] col_index;
    }
}



template<class T>
void CSRMatrix<T>::printMatrix() {
    cout << "non-zero values" << endl;
    //print the non-zero entries
    for (int j = 0; j < this->nnzs; j++) {
        cout << this->values[j] << " ";
    }
    cout << endl;
    cout << "row_position" << endl;
    for (int j = 0; j < this->rows + 1; j++) {
        cout << this->row_position[j] << " ";
    }
    cout << endl;
    cout << "col_index:" << endl;
    for (int j = 0; j < this->nnzs; j++) {
        cout << this->col_index[j] << " ";
    }
    cout << endl << endl;
}

//do multiplication between matrix and vector
template<class T>
void CSRMatrix<T>::matVecMult(T* input, T* output) {
    if (input == nullptr || output == nullptr)
    {
        std::cerr << "Input or output haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

//    int val_counter = 0;
    for (int i = 0; i < this->rows; i++) {
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++) {
            output[i] += this->values[val_index] * input[this->col_index[val_index]];
        }
    }
}

//do multiplication between matrix and matrix
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output)
{

    // Check our dimensions match
    if (this->cols != mat_right.rows)
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }

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
        std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;

    }

    // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
}

//check if current A*x is close to b
template <class T>
bool CSRMatrix<T>::if_finish(T* mat_b, T* output) {
    T tol = 0.000001;
    T* cal_out = new T[this->rows];
    T res = 0;

    // calculate the resident
    for (int i = 0; i < this->rows; i++) {
        cal_out[i] = 0;
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++) {
            cal_out[i] += this->values[val_index] * output[this->col_index[val_index]];
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
void CSRMatrix<T>::Gauss_seidel(T* mat_b, T* output) {

    // Set output to zero before hand
    for (int i = 0; i < this->rows; i++)  output[i] = 0;

    T temp;
    int k;

    //start iterations
    for (k = 0; k < this->it_max; k++) {
        for (int i = 0; i < this->rows; i++) {
            temp = mat_b[i];
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++) {
                if (i == this->col_index[val_index]) {
                    continue;
                }
                temp = temp - this->values[val_index] * output[this->col_index[val_index]];
            }
            //we will update and iterate at the same time
            output[i] = temp / this->get_A_value(i, i);
        }

        if (this->if_finish(mat_b, output)) {
            break;
        }
    }

    // print the result
    if (k < this->it_max) {
        cout << "!!!!!!!! we got Gauss-Seidel solution for the sparse matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for Gauss-Seidel solution for the sparse matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, Gauss-Seidel solution cannot solve it!!" << endl << endl;
    }

}

template <class T>
void CSRMatrix<T>::Jacobi_solve(T* mat_b, T* output) {

    T* temp = new T[this->rows];

    // set output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

    int k;

    //start iterations
    for (k = 0; k < this->it_max; k++) {
        for (int i = 0; i < this->rows; i++) {
            temp[i] = mat_b[i];
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++) {
                if (i == this->col_index[val_index]) {
                    continue;
                }
                temp[i] = temp[i] - this->values[val_index] * output[this->col_index[val_index]];
            }

            temp[i] = temp[i] / this->get_A_value(i, i);
        }
        //update the x values after this iteration finishes
        for (int i = 0; i < this->rows; i++) {
            output[i] = temp[i];
        }
        // judge if finished
        if (this->if_finish(mat_b, output)) {
            break;
        }
    }

    delete[] temp;

    //print the result
    if (k < this->it_max) {
        cout << "!!!!!!!! we got Jacobi solution for the sparse matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for Jacobi solution for the sparse matrix: " << k << endl << endl;
    }
    else {
        cout << "!!!!!!!! Bad input matrix, Jacobi solution cannot solve it!!" << endl << endl;
    }

}

template <class T>
void CSRMatrix<T>::BiCGSTAB(T* mat_b, T* output) {
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
        cout << "!!!!!!!! we got BiCGSTAB solution for the sparse matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << output[i] << " ";
        }
        cout << endl << "the iterations for BiCGSTAB solution for the sparse matrix: " << k << endl << endl;
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

template <class T>
void CSRMatrix<T>::GMRES(T* mat_b, T* output) {
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

    cout << "!!!!!!!! we got GMRES solution for the sparse matrix(only the first ten x values):" << endl;
    for (int i = 0; i < 10; i++)
    {
        cout << output[i] << " ";
    }
    cout << endl << "the iterations for GMRES solution for the sparse matrix: " << k << endl << endl;

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


// get A value following the sparse storing rule
template <class T>
T CSRMatrix<T>::get_A_value(int i, int j) {
    for (int k = this->row_position[i]; k < this->row_position[i + 1]; k++) {
        if (this->col_index[k] == j) return this->values[k];
    }
    return 0;
}
// add L value following the sparse storing rule
template <class T>
void CSRMatrix<T>::set_L_value(int i, int j, T value) {
    if (value != 0) {
        this->L_value.push_back(value);
        this->L_col_index.push_back(j);
        for (int k = i + 1; k <= this->rows; k++) this->L_row_pos[k] += 1;
    }
}
// get L value following the sparse storing rule
template <class T>
T CSRMatrix<T>::get_L_value(int i, int j) {
    for (int k = this->L_row_pos[i]; k < this->L_row_pos[i + 1]; k++) {
        if (this->L_col_index[k] == j) return this->L_value[k];
    }
    return 0;
}
template <class T>
void CSRMatrix<T>::Cholesky_forward_substitution(T* rhs, T* result)
{
    // initialize result as 0
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }
    // L * y = b
    // the y calculated at this step will be used in the next step
    double s;
    for (int k = 0; k < this->rows; k++) {
        s = 0.0;
        for (int j = 0; j < k; j++)
        {
            s += this->get_L_value(k, j) * result[j];
        }
        result[k] = (rhs[k] - s) / get_L_value(k, k);
    }
}
template <class T>
void CSRMatrix<T>::Cholesky_backward_substitution(T* rhs, T* result)
{
    // initialize result as 0
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }
    // L.T*x=y
    // the x calculated at this step will be used in the next step
    for (int k = this->rows - 1; k > -1; k--) {
        double s(0.);
        for (int j = k + 1; j < this->rows; j++)
        {
            s += get_L_value(j, k) * result[j];
        }
        result[k] = (rhs[k] - s) / get_L_value(k, k);
    }
}

template <class T>
void CSRMatrix<T>::Cholesky(T* rhs, T* result) {

    //determing whether this matrix is symmetric or nor
    int x = 0;
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            if (this->get_A_value(i, j) != this->get_A_value(j, i)) {
                x = 1;
            }
        }
    }
    if (x == 1) {
        cout << "input matrix is not symmetric, we can't use Cholesky method!!" << endl;
    }

    else {
        T L1, L2, A;

        //clear L
        this->L_value.clear();
        this->L_col_index.clear();
        this->L_row_pos.clear();
        for (int k = 0; k <= this->rows; k++) this->L_row_pos.push_back(0);

        //decompose original matrix
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < i + 1; j++) {
                if (i == j) {
                    L1 = 0;
                    //set L1=sum(L[j,k]^2) when k<j
                    for (int k = 0; k < j; k++) L1 += pow(this->get_L_value(j, k), 2);
                    //set L[j,j]=(A[j,j]-L1)^2
                    this->set_L_value(j, j, pow((this->get_A_value(j, j) - L1), 0.5));
                }
                else {
                    L2 = 0;
                    //set L2=sum(L[i,k]*L[j,k]) when k<j
                    for (int k = 0; k < j; k++) L2 += this->get_L_value(i, k) * this->get_L_value(j, k);
                    //set L[i,j]=1/L[j,j]*(A[i,j]-L2)
                    this->set_L_value(i, j, 1 / this->get_L_value(j, j) * (this->get_A_value(i, j) - L2));
                }
            }
        }

        auto* for_result = new T[this->rows];
        //apply forward substitution to lower triangular matrix
        this->Cholesky_forward_substitution(rhs, for_result);
        //apply backward substitution to upper triangular matrix
        this->Cholesky_backward_substitution(for_result, result);

        cout << "!!!!!!!! we got Cholesky solution for the sparse matrix(only the first ten x values):" << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << result[i] << " ";
        }
        cout << endl << endl;

        delete[] for_result;
    }

}

template <class T>
void CSRMatrix<T>::update_A_value(int i, int j, T value) {
    for (int k = this->row_position[i]; k < this->row_position[i + 1]; k++) {
        if (this->col_index[k] == j)  this->values[k] = value;
    }

}


template <class T>
void CSRMatrix<T>::LU_forward_substitution(T* rhs, T* result)
{
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }

    double s;
    for (int k = 0; k < this->rows; k++) {
        s = 0.0;
        for (int j = 0; j < k; j++)
        {

            s += this->get_L_value(k, j) * result[j];
        }
        result[k] = (rhs[k] - s) / get_L_value(k, k);
    }
}

template <class T>
void CSRMatrix<T>::LU_backward_substitution(T* rhs, T* result)
{
    for (int i = 0; i < this->rows; i++) {
        result[i] = 0;
    }
    for (int k = this->rows - 1; k > -1; k--) {
        double s(0.);
        for (int j = k + 1; j < this->rows; j++)
        {
            s += get_A_value(k, j) * result[j];
        }
        result[k] = (rhs[k] - s) / get_A_value(k, k);
    }
}

template <class T>
void CSRMatrix<T>::LU_solve(T* rhs, T* result)
{

    this->L_value.clear();
    this->L_col_index.clear();
    this->L_row_pos.clear();

    for (int k = 0; k <= this->rows; k++) {
        this->L_row_pos.push_back(0);
    }

    for (int k = 0; k < this->rows; k++)
    {

        for (int i = 0; i < k + 1; i++)
        {
            if (i == k) this->set_L_value(k, i, 1);
            else {

                this->set_L_value(k, i, get_A_value(k, i) / get_A_value(i, i));
            }
        }
    }

    for (int k = 0; k < this->rows - 1; k++)
    {
        for (int i = k + 1; i < this->cols; i++)
        {
            double temp_s = get_A_value(i, k) / get_A_value(k, k);
            for (int j = k; j < this->cols; j++)
            {
                this->update_A_value(i, j, (get_A_value(i, j) - temp_s * get_A_value(k, j)));
            }
        }

    }
    auto* for_result = new T[this->rows];
    this->LU_forward_substitution(rhs, for_result);
    this->LU_backward_substitution(for_result, result);

    cout << "!!!!!!!! we got LU solution for the sparse matrix(only the first ten x values):" << endl;
    for (int i = 0; i < 10; i++)
    {
        cout << result[i] << " ";
    }
    cout << endl << endl;

    delete[] for_result;

}

template <class T>
void CSRMatrix<T>::find_unique_p(std::vector<bool> check_list, std::vector<int>& unique_list)
{
    //loop over all the columns in °∞check_list°± to check how many non-zeros are inside;
    //if only one element in that column, mark down the corresponding row index in °∞unique_list°±
    auto* count = new int[this->cols];
    for (int i = 0; i < this->cols; i++) //set initial count array values 0
    {
        count[i] = 0;
    }
    for (int i = 0; i < this->nnzs; i++)//count the number of non-zero value of each cols
    {
        if (col_index[i] != -1)//make sure the col_index[i] is valid
        {
            if (check_list[col_index[i]])
            {
                count[this->col_index[i]] += 1;
            }
        }
    }
    for (int n = 0; n < this->cols; n++)//transfer to unique_list
    {
        if (count[n] == 1)
        {
            for (int j = 0; j < this->nnzs; j++)
            {
                if (this->col_index[j] == n)
                {
                    unique_list[n] = j;
                    break;
                }
            }
        }
    }
    delete[] count;
}
template <class T>
void CSRMatrix<T>::sort_mat(double* rhs) {
    // sort the matrix and corresponding rhs in order to eliminate the non-zeros on diagonal;
    // and also, do partial pivoting to make the pivot element more significant;
    auto* temp_mat = new CSRMatrix<T>(this->rows, this->cols, this->nnzs, true);
    auto* temp_rhs = new double[this->rows];
    for (int i = 0; i < nnzs; i++)
    {
        temp_mat->values[i] = 0;
        temp_mat->col_index[i] = 0;
        temp_mat->row_position[i] = 0;
    }

    int* c = new int[this->rows]; // number of values in a row
    for (int i = 0; i < this->rows; i++)
    {
        c[i] = 0;
    }

    double** a = new double* [this->rows]; // values in a row
    int** b = new int* [this->rows]; // column index of a row

    std::vector <bool> check_list(this->cols, true);//use °∞check_list°± to track columns to be sorted

    //when check_list is not empty (all false), continue the while loop
    while (!(std::none_of(check_list.begin(), check_list.end(), [](bool v) { return v; })))
    {
        //use °∞unique_list°± to find columns that have only one non-zero
        std::vector <int> unique_list(this->cols, -1);

        this->find_unique_p(check_list, unique_list);//update unique_list
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
                    // attention!!!
                    for (int i = 0; i < this->rows + 1; i++)
                    {

                        if (this->row_position[i] < unique_list[j] + 1 && unique_list[j] < this->row_position[i + 1])
                        {

                            a[j] = new double[this->row_position[i + 1] - this->row_position[i]];
                            b[j] = new int[this->row_position[i + 1] - this->row_position[i]];
                            c[j] = this->row_position[i + 1] - this->row_position[i];
                            //              std::cout<<°∞numbers per row is °∞<<c[j]<<°± at row °∞<<j<<std::endl;
                            for (int idx = this->row_position[i]; idx < this->row_position[i + 1]; idx++)
                            {
                                // add the col_index, exactly the same as original col_index, but on different location
                                // add the value value[i*cols+col] to value°Ø[j*cols+col] (this expression is in dense format)
                                // all the row_position value (row[j] and after) add 1
                                for (int itr_row_pos = j + 1; itr_row_pos < this->cols + 1; itr_row_pos++)
                                {
                                    temp_mat->row_position[itr_row_pos] += 1;
                                }
                                a[j][idx - row_position[i]] = this->values[idx];
                                b[j][idx - row_position[i]] = this->col_index[idx];

                                this->col_index[idx] = -1;

                            }
                            temp_rhs[j] = rhs[i];
                            //only run once; break
                            break;
                        }

                    }
                    // column [j] is finished
                    check_list[j] = false;
                }
            }
        }

        //if no unique element found, fill the 1st available column with the row
        //where value that will be the diagonal element is maximum in that row
        //and remove the column from check_list;
        else
        {
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

                        if (this->get_A_value(row, j) != 0 && this->get_A_value(row, j) != -1)
                        {
                            diff = 0;
                            for (int aa = 0; aa < this->cols; aa++)
                            {
                                diff += abs(this->get_A_value(row, aa));
                            }

                            diff = abs(this->get_A_value(row, j)) / (diff - abs(this->get_A_value(row, j)));
                            if (diff >= max_value)
                            {
                                index_row = row;
                                max_value = diff;
                            }
                        }
                    }
                    //now we fill temp_mat row j with value in row °∞index_row°±
                    a[j] = new double[this->row_position[index_row + 1] - this->row_position[index_row]];
                    b[j] = new int[this->row_position[index_row + 1] - this->row_position[index_row]];
                    c[j] = this->row_position[index_row + 1] - this->row_position[index_row];

                    for (int idx = this->row_position[index_row]; idx < this->row_position[index_row + 1]; idx++)
                    {
                        // all the row_position value (row[j] and after) add 1
                        for (int itr_row_pos = j; itr_row_pos < this->cols; itr_row_pos++)
                        {
                            temp_mat->row_position[itr_row_pos + 1] += 1;
                        }
                        // add the col_index, exactly the same as original col_index, but on different location
                        a[j][idx - this->row_position[index_row]] = this->values[idx];
                        b[j][idx - this->row_position[index_row]] = this->col_index[idx];

                        this->col_index[idx] = -1;



                    }
                    temp_rhs[j] = rhs[index_row];
                    //only run once;
                    check_list[j] = false;
                    break;
                }
            }
        }
    }
    //now set all the values in origin matrix

    //row_position is straight forward
    //col and value needs to be read and then stored
    std::vector<double> temp_valuess;
    std::vector<int> temp_colss;
    for (int i = 0; i < this->rows; i++)
    {
        int temp_c = c[i];

        for (int j = 0; j < temp_c; j++)
        {
            temp_valuess.push_back(a[i][j]);
            temp_colss.push_back(b[i][j]);
        }
        this->row_position[i] = temp_mat->row_position[i];
    }

    this->row_position[this->rows] = temp_mat->row_position[this->rows];

    //now store values and col_index
    for (int i = 0; i < this->nnzs; i++)
    {
        this->values[i] = temp_valuess[i];
        this->col_index[i] = temp_colss[i];
    }

    //store rhs
    for (int i = 0; i < this->rows; i++)
    {
        rhs[i] = temp_rhs[i];
    }


    delete temp_mat;
    delete[] temp_rhs;
    delete[] c;
    for (int i = 0; i < this->rows; i++) {
        delete[] a[i];
        delete[] b[i];
    }
    delete[] a;
    delete[] b;
}

