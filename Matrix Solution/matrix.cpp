#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

const int n = 4;
const double val[4][4] = {{1,-2,1,-2},{1,-2,1,0},{1,-2,0,0},{1,0,0,0}};
const double RHS[4][4] = {{1,-.5,1,-.5},{-.5,-2,-.5,0},{1,-.5,0,0},{-.5,0,0,0}};
int LU_decomposition(const int n, std::vector<std::vector<double>> & a, std::vector<int> & swap_indices, int & swap_count);
int LU_solve(const int n, const std::vector<std::vector<double>> & LU, const std::vector<int> & swap_indices, const int swap_count, const std::vector<double> & rhs, std::vector<double> & x);
int Tridiagonal_solve(const int n, const vector<double> & a, const vector<double> & b, const vector<double> & c, const vector<double> & rhs, vector<double> & x);
void print(vector<double>& a, int x);
void fill(vector<double>& a, vector<double>& b, vector<double>& c);

int main(){
    vector<vector<double>> LU(n);
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n ; j++){
            LU[i].push_back(val[i][j]);
        }
    }
    
    vector<int> swap_indices;
    int swap_count = 0;
    LU_decomposition(n, LU, swap_indices, swap_count);
    
    vector<double> rhs1, rhs2, rhs3, rhs4;
    for(int i = 0; i< n; i++){
        rhs1.push_back(RHS[i][0]);
        rhs2.push_back(RHS[i][1]);
        rhs3.push_back(RHS[i][2]);
        rhs4.push_back(RHS[i][3]);
    }
    vector<double> x(4);
    LU_solve(n,LU,swap_indices,swap_count,rhs1,x);
    print(x, n);
    LU_solve(n,LU,swap_indices,swap_count,rhs2,x);
    print(x,n);
    LU_solve(n,LU,swap_indices,swap_count,rhs3,x);
    print(x,n);
    LU_solve(n,LU,swap_indices,swap_count,rhs4,x);
    print(x,n);

    vector<double> a,b,c, y(5);
    x.clear();
    x.resize(5);
    for(int i = 0; i<5; i++){
        y[i] = i+1;
    }
    fill(a, b, c);
    Tridiagonal_solve(5, a, b, c, y, x);
    print(x, 5);
}

void fill(vector<double>& a, vector<double>& b, vector<double>& c){
    b.push_back(0);
    a.push_back(3);
    c.push_back(1);
    for(int i = 0;i<3;i++){
        b.push_back(1);
        a.push_back(3);
        c.push_back(1);
    }
    b.push_back(1);
    a.push_back(3);
    c.push_back(0);
}

void print(vector<double>& a, int x){
    for(int i =0; i<x; i++){
        cout << a[i] << "  ";
    }
    cout << endl;
}

int LU_decomposition(const int n, std::vector<std::vector<double>> & a, std::vector<int> & swap_indices, int & swap_count){
    const double tol = 1.0e-14; // tolerance for pivot
    swap_indices.clear();
    swap_count = 0;
    if (n < 1) return 1; // fail
    if (a.size() < n) return 1; // fail
    int i=0;
    int j=0;
    swap_indices.reserve(n);
    for (i = 0; i < n; ++i) {
        swap_indices.push_back(i); // initialize to (0,...n,-1)
    }
    for (i = 0; i < n; ++i) {
        // step 1: calculate scaled coeffs, compare for max pivot
        int pivot_row = i;
        double max_pivot = 0;
        for (j=i; j < n; ++j) {
            int k = i;
            double a_hat = std::abs(a[j][k]);
            if (a_hat > 0.0) {
            for (k=i+1; k < n; ++k) {
                double tmp = std::abs(a[j][k]);
                if (a_hat < tmp) a_hat = tmp;
                }
                double a_scaled = std::abs(a[j][i])/a_hat; // note: a_ji not a_ij
                if (max_pivot < a_scaled) {
                    max_pivot = a_scaled;
                    pivot_row = j;
                }
            }
        }
        // step 2: if max pivot <= tol return fail (inconsistent or not linearly independent)
        if (std::abs(a[pivot_row][i]) <= tol) {
            swap_indices.clear();
            swap_count = 0;
            return 1; // fail
        }
        // step 3: found the max pivot, swap rows if row != i
        if (pivot_row != i) {
            // swap the whole row, including lower triangular entries
            // keep track of number of swaps
            // update the swap index
            ++swap_count;
            int si = swap_indices[i];
            swap_indices[i] = swap_indices[pivot_row];
            swap_indices[pivot_row] = si;
            for (j=0; j < n; ++j) {
                double tmp = a[i][j];
                a[i][j] = a[pivot_row][j];
                a[pivot_row][j] = tmp;
            }
        }
        double inv_pivot = 1.0/a[i][i]; // this is nonzero
        // step 4: eliminate column i from rows i+1 <= j < n, overwrite make LU matrix
        for (j=i+1; j < n; ++j) {
            double multiplier = a[j][i]*inv_pivot;
            a[j][i] = multiplier; // store multiplier in lower triangular matrix
            for (int k=i+1; k < n; ++k) {
                a[j][k] -= multiplier*a[i][k]; // calculate upper triangular matrix
            }
        }
    }
    return 0;
}

int LU_solve(const int n, const std::vector<std::vector<double>> & LU, const std::vector<int> & swap_indices, const int swap_count, const std::vector<double> & rhs, std::vector<double> & x){
    x.clear();
    if (n < 1) return 1; // fail
    std::vector<double> y(n, 0.0); // temporary storage
    // forward substitution Ly = rhs
    int i=0;
    for (i = 0; i < n; ++i) {
        int si = swap_indices[i];
        double sum = rhs[si];
        for (int j = 0; j < i; ++j) {
            sum -= LU[i][j] * y[j];
        }
        y[i] = sum;
    }
    // backward substitution Ux = y
    for (i = n-1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i+1; j < n; ++j) {
            sum -= LU[i][j] * x[j];
        }
        x[i] = sum / LU[i][i];
    }
    return 0;
}

int Tridiagonal_solve(const int n, const vector<double> & a, const vector<double> & b, const vector<double> & c, const vector<double> & rhs, vector<double> & x){
    const double tol = 1.0e-14;
    x.clear();

    if ((n < 1) || (a.size() < n) || (b.size() < n) || (c.size() < n) || (rhs.size() < n)){
        return 1;
    }

    double alpha[n]; // temporary storage
    x.resize(n, 0.0);

    // initial equation i = 0
    int i = 0;
    double gamma = a[i];
    if (abs(gamma) <= tol){
        return 1;
    }


    x[i] = rhs[i]/gamma;
    alpha[i] = c[i]/gamma;

    // forward pass: elimination
    for (i = 1; i < n-1; ++i) {
        gamma = a[i] - b[i]*alpha[i-1];
        if (abs(gamma) <= tol){
            return 1;
        }

        x[i] = (rhs[i] - b[i]*x[i-1])/gamma;
        alpha[i] = c[i]/gamma;
    }

    // solve final equation i = n-1
    i = n-1;
    gamma = a[i] - b[i]*alpha[i-1];

    if (abs(gamma) <= tol){
        return 1;
    }

    x[i] = (rhs[i] - b[i]*x[i-1])/gamma;

    // backward substitution
    for (i = n-2; i >= 0; --i) {
        x[i] -= alpha[i]*x[i+1];
    }
    
    return 0;
}