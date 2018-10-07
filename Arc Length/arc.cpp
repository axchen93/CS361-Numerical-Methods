#include<iostream>
#include<cmath>
using namespace std;

typedef double(*Function)(const double x);
const double alpha = 0.2302;
const double beta = 0.7199;

double func(const double x);
void funcPrime(const double x, const double y, double &fprime);
double arcLengthFunc(const double x);
void trapezoid(Function f, const double a, const double b, const int n, double& length);
void trapezoid2(Function f, const double a, const double b, const int n, double& length);
void bisection(double x_low, double x_high, const double target, double & x);

int main(){
    double x_high = 10;
    double x_low = 20;
    double x_target;
    bisection(x_low, x_high, 1.0+beta, x_target);
    cout << "a = " << x_target << endl << "b = " << func(x_target) << endl;
    cout << "Target arc lenght  = " << 1.0 + beta << endl;
    double arc_f;
    trapezoid(arcLengthFunc, x_target, func(x_target), 1000, arc_f);
    cout << "Arc Length from a  to b is " << arc_f << endl;

    double area = 0;
    double previousArea = 0;
    double tol = 1.0e-5;
    double a = x_target;
    double b = func(x_target);
    for(int i = 1; i<100; i++){
        int n = pow(2,i);
        trapezoid(func, a, b, n, area);
        cout << "n = " << n << " " << area <<endl;
        if(abs(area-previousArea) <= tol) break;
        previousArea = area;
    }

    double rectangleArea = (b-a) * b;
    cout << "Area of the curve from a to b is " << area << endl;
    cout << "Area of the rectangle arouny the integral is " << rectangleArea << endl;
    cout << "Area of the region enclosed by rectangle and curve is " << rectangleArea - area << endl; 

    previousArea = 0;
    for(int i = 1; i<100; i++){
        int n = pow(2,i);
        trapezoid2(func, a, b, n, area);
        cout << "n = " << n << " " << area <<endl;
        if(abs(area-previousArea) <= tol) break;
        previousArea = area;
    }

}

double func(const double x){
    double a = 1 - pow(x, -alpha);
    return pow(a, 1.0/-alpha);
}
void funcPrime(const double x, const double y, double &fprime){
    double a = -alpha * pow(y, alpha+1);
    double b = alpha * pow(x, alpha+1);
    fprime = a/b;
}
double arcLengthFunc(const double x){
    double y, fprime;
    y = func(x);
    funcPrime(x, y, fprime);

    fprime*=fprime;
    return sqrt(1 + fprime);
}
void trapezoid(Function f, const double a, const double b, const int n, double& length){
    double h = (b-a) / n;
    double sum = 0.5 * (f(a) + f(b));
    double x_i = a;
    for(int i = 1; i<n; i++){
        x_i += h;
        sum += f(x_i); 
    }
    length = sum*h;
}
void trapezoid2(Function f, const double a, const double b, const int n, double& length){
    double h = (b-a) / n;
    double sum = 0.5 * ( (b-f(a)) + (b-f(b)));
    double x_i = a;
    for(int i = 1; i<n; i++){
        x_i += h;
        sum += (b - f(x_i)); 
    }
    length = sum*h;
}
void bisection(double x_low, double x_high, const double target, double & x){
    double tol = 1.0e-5;
    int maxIteration = 100;
    double f_low, f_high, y_high, y_low, diff_low, diff_high;

    y_high = func(x_high);
    trapezoid(arcLengthFunc, x_high, y_high, 1000, f_high);
    diff_high = f_high - target;
    if(abs(diff_high) <= tol){
        x = x_high;
        return;
    }

    y_low = func(x_low);
    trapezoid(arcLengthFunc, x_low, y_low, 1000, f_low);
    diff_low = f_low - target;
    if(abs(diff_low) <= tol){
        x = x_low;
        return;
    }
    
    double y_mid, f_mid, diff;
    for(int i = 1; i<maxIteration; i++){
        x = (x_high + x_low) / 2.0;
        y_mid = func(x);
        trapezoid(arcLengthFunc, x, y_mid, 1000, f_mid);

        diff = f_mid - target;

        if(abs(diff) <= tol) return;

        if(diff*diff_high > 0.0 ) x_high = x;
        else x_low = x;

        if( abs(x_high - x_low) <= tol) return;
    }

    x=0;
    return;

}