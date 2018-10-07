#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
using namespace std;

const double alpha = 0.2302;
const double beta = 0.7199;

int Tridiagonal_solve(const int n, const vector<double> & a, const vector<double> & b, const vector<double> & c, const vector<double> & rhs, vector<double> & x);
void fillLHS(const int n, vector<double> & a, vector<double> & b, vector<double> & c);
void fillRHS(const int n, vector<double> & rhs);
void toFile(const int n, const vector<double> & y, ofstream& output );

int main(){
    ofstream output;
    output.open("output.txt");
    const int n = 19000;
    vector<double> a, b, c ;
    fillLHS(n, a, b, c);

    vector<double> rhs (n+1, 0.0);
    fillRHS(n, rhs);

    vector<double> y (n+1, 0.0);
    cout << Tridiagonal_solve(n+1, a, b, c, rhs, y);
    toFile(n, y, output);
    output.close();
}

void toFile(const int n, const vector<double> & y, ofstream& output ){
    const double x_0 = 1;
    const double x_n = 20;
    const double h = (x_n - x_0) / n;
    double x_i;
    for(int i = 0; i<=n; i++){
        x_i = x_0 + i*h;
        output << x_i << ", " << y[i] << endl;
    }
}

void fillRHS(const int n, vector<double> & rhs){
    const double x_0 = 1;
    const double x_n = 20;
    const double h = (x_n - x_0) / n;

    rhs[0] = alpha;
    rhs[n] = h*beta;
}

void fillLHS(const int n, vector<double> & a, vector<double> & b, vector<double> & c){
    const double x_0 = 1;
    const double x_n = 20;
    const double h = (x_n - x_0) / n;

    a.push_back(1.0);
    b.push_back(0.0);
    c.push_back(0.0);

    double x_i = x_0;
    for(int i = 1; i<n; i++){
        x_i += h;
        double a_i = (h*h)*(x_i*x_i) - (h*h) - (2*x_i*x_i);
        double b_i = (x_i*x_i) - (h/2.0)*x_i;
        double c_i = (x_i*x_i) + (h/2.0)*x_i;
        a.push_back(a_i);
        b.push_back(b_i);
        c.push_back(c_i);
    }

    b.push_back(-1.0);
    a.push_back(1.0);
    c.push_back(0.0);

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