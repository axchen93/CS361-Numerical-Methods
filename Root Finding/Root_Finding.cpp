#include<iostream>
#include<cmath>
using namespace std;

double func(double x);
void func(double x, double &f, double &fprime);
int root_bisection(double target, double tol_f, double tol_x, int max_iter, double x_low, double x_high, double & x, int & num_iter);
int root_NR(double target, double tol_f, double tol_x, int max_iter, double x0,
double & x, int & num_iter);


int main(){
    double x = 0;
    int num_iter=0;

    root_bisection(3.0, 1.0e-6 , 1.0e-6, 100, 0.0, 5.0, x, num_iter);
    cout << "Bisection: " << x << "\nNumber of iteration: " << num_iter << endl;
    cout <<  "f(x) = " <<func(x) << endl;

    root_NR(3.0, 1.0e-6, 1.0e-6, 100, 1, x, num_iter);
    cout << "Newton: " << x << "\nNumber of iteration: " << num_iter << endl;
    cout <<  "f(x) = " <<func(x) << endl;
}


//Root finding by Bisection
int root_bisection(double target, double tol_f, double tol_x, int max_iter, double x_low, double x_high, double & x, int & num_iter){
    double y_low, y_high, diff_y_low, diff_y_high; 
    x = 0, num_iter = 0;

    y_low = func(x_low);
    diff_y_low = y_low - target;
    if(abs(diff_y_low) <= tol_f){
        x = x_low;
        return 0;
    }

    y_high = func(x_high);
    diff_y_high = y_high - target;
    if(abs(diff_y_high) <= tol_f){
        x = x_high;
        return 0;
    }

    if(diff_y_high * diff_y_low > 0.0){
        x = 0;
        return 1;
    }

    double y, diff_y;
    for(num_iter = 1; num_iter < max_iter; num_iter++){
        x = (x_low + x_high) / 2.0;
        y = func(x);
        diff_y = y - target;

        if(abs(diff_y) <= tol_f) return 0;
        
        if(diff_y * diff_y_low > 0.0) x_low = x;
        else x_high = x;

        if(x_high - x_low <= tol_x) return 0;
    }
    
    x = 0;
    num_iter = max_iter;
    return 1;
}

//Root finding by Newton-Raphson
int root_NR(double target, double tol_f, double tol_x, int max_iter, double x0,
double & x, int & num_iter){
    const double tol_fprime = 1.0e-12;
    double f = 0;
    double fprime = 0;
    x = x0;

    double diff_f, delta_x;
    for (num_iter = 1; num_iter < max_iter; ++num_iter) {
        func(x, f, fprime);

        diff_f = f - target;
        if(abs(diff_f) <= tol_f) return 0;
        if(abs(fprime) <= tol_fprime) {
            x = 0;
            return 1;
        }

        delta_x = diff_f/fprime;
        if(abs(delta_x) <= tol_x) return 0;

        x -= delta_x;
    }

    x = 0;
    num_iter = max_iter;
    return 1;
}

//#define func_quadratic
#define func_cubic
//#define func_cum_norm
#ifdef func_quadratic
double func(double x){
return x*x;
};
void func(double x, double &f, double &fprime){
f=x*x;
fprime=2.0*x;
};
#endif // func_quadratic
#ifdef func_cubic
double func(double x){
return x*x*x;
};
void func(double x, double &f, double &fprime){
f=x*x*x;
fprime=3.0*x*x;
};
#endif // func_cubic
#ifdef func_cum_norm
double cum_norm(double x)
{
const double root = sqrt(0.5);
return 0.5*(1.0 + erf(x*root));
}
double func(double x)
{
return cum_norm(x);
}
void func(double x, double &f, double &fprime)
{
const double pi = 4.0*atan2(1.0,1.0);
f = cum_norm(x);
fprime = exp(-0.5*x*x)/sqrt(2.0*pi);
}
#endif // func_cum_norm

