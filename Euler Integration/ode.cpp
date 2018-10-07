#include<iostream>
#include <iomanip>
#include<cmath>
using namespace std;

namespace _15_1{
    void f(double x, double y, double & g);
    void exact (double x, double & exact);
    void Euler_explicit(double h, double n, double& y_in, double & y_out);
    void Euler_explicit_reverse(double h, double n, double& y_in, double & y_out);
    void print(double explict, double exact, int n);
    void run();
}

namespace _15_2{
    void Euler_implicit(double h, double n, double& y_in, double & y_out);
    void Euler_implicit_reverse(double h, double n, double& y_in, double & y_out);
    void run();
}

namespace _15_3{
    void Euler_explicit(double h, double n, double& y_in, double & y_out);
    void Euler_implicit(double h, double n, double& y_in, double & y_out);
    void exact(double x, double& ex);
    void run();
}

namespace _15_4{
    void Euler_explicit(double h, double n, double& y_in, double & y_out);
    void Euler_implicit(double h, double n, double& y_in, double & y_out);
    void exact(double x, double& ex);
    void run();
}

int main(){
    _15_1::run();
    _15_2::run();
    _15_3::run();
    _15_4::run();
}

namespace _15_1{
    void f(double x, double y, double & g){
        const double constant = 2;
        g = x - (constant * y);
    }

    void exact(double x, double & exact){
        const double constant = 2;
        double a = (1 + (1 / (constant*constant))) * exp( (-1)*constant*x );
        double b = ( x / constant);
        double c = (1 / (constant*constant));
        exact = a + b - c;
    }

    void Euler_explicit(double h, double n, double& y_in, double & y_out){
        double g = 0;
        for(int i = 0 ; i<n; i++){
            double x = i * h;
            f(x, y_in, g);
            y_out = y_in + h*g;
            y_in = y_out;
        }
    }
    void Euler_explicit_reverse(double h, double n, double& y_in, double & y_out){
        double g = 0;
        double x = 1.0;
        for(int i = n ; i<2*n; i++){
            x+=h;
            f(x, y_in, g);
            y_out = y_in + h*g;
            y_in = y_out;
        }
    }
    void print(double explict, double exact, int n){
        cout << "n = " << n << " " <<fixed << setprecision(5) << exact << " ";
        cout << fixed << setprecision(5) << explict << " ";
        double diff = exact - explict;
        cout << fixed << setprecision(5) << diff << " ";
        cout << fixed << setprecision(2) << n*diff << endl;
    }

    void run(){
        double ex;
        double y_out, y_in;
        int n = 10;
        for(int i = 0; i< 3; i++){
            double h = 1.0 / n;
            exact(1, ex);
            y_out = 0.0;
            y_in  = 1.0;
            cout << "Eluer Explicit Foward: n = " << n << endl;
            Euler_explicit(h, n, y_in, y_out);
            print(y_out, ex, n);

            exact(0, ex);
            h = -h;
            cout << "Eluer Explicit Reverse: n = " << 2*n << endl;
            Euler_explicit_reverse(h, n, y_in, y_out);
            print(y_out, ex, 2*n);

            n*=10;
            cout << endl;
        }
    }
}

namespace _15_2{
    void Euler_implicit(double h, double n, double& y_in, double & y_out){
        const double c = 2;
        for(int i = 0 ; i<n; i++){
            double x = (i+1) * h;
            y_out = ( y_in + (h*x) ) / ( 1 + (h*c));
            y_in = y_out;
        }
    }
    void Euler_implicit_reverse(double h, double n, double& y_in, double & y_out){
        const double c = 2;
        double x = 1.0;
        for(int i = 0 ; i<n; i++){
            x += h;
            y_out = ( y_in + (h*x) ) / ( 1 + (h*c));
            y_in = y_out;
        }
    }
    void run(){
        cout << endl;
        double ex;
        double y_out, y_in;
        int n = 10;
        for(int i = 0; i< 3; i++){
            double h = 1.0/n;
            _15_1::exact(1, ex);
            y_out = 0.0;
            y_in = 1.0;
            cout << "Eluer Implicit Foward: n = " << n << endl;
            Euler_implicit(h, n, y_in, y_out);
            _15_1::print(y_out, ex, n);

            h = -h;
            _15_1::exact(0, ex);
            cout << "Eluer Implicit Reverse: n = " << 2*n << endl;
            Euler_implicit_reverse(h, n, y_in, y_out);
            _15_1::print(y_out, ex, 2*n);
            n*=10;
            cout << endl;
        }
    }
}

namespace _15_3{
    void exact(double x, double& ex){
        ex = x*x;
    }

    void Euler_explicit(double h, double n, double& y_in, double & y_out){
        double x = 0.0;
        for(int i = 0; i < n; i++){
            x+=h;
            y_out = y_in + (2*h*x);
            y_in = y_out;
        }
    }

    void Euler_implicit(double h, double n, double& y_in, double & y_out){
        double x = h;
        for(int i = 0; i < n; i++){
            x+=h;
            y_out = y_in + (2*h*x);
            y_in = y_out;
        }
    }

    void run(){
        double ex;
        exact(1, ex);
        int n = 10;
        double y_in, y_out;
        cout << endl << "dy/dx = x^2" << endl;
        for(int i = 0; i<3; i++){
            y_in = y_out = 0;
            double h = 1.0/n;
            Euler_explicit(h, n, y_in, y_out);
            double explicitResult = y_out;
            y_in = y_out = 0;
            Euler_implicit(h, n, y_in, y_out);
            cout << "n = " << n << " "<<  fixed <<  setprecision(5) << explicitResult << " " <<  fixed <<  setprecision(5) << y_out << " " <<  fixed <<  setprecision(2) << n*(ex-explicitResult)<< " " <<  fixed <<  setprecision(2) << n*(ex-y_out) << endl;
            n*=10;  
        }
    }
}

namespace _15_4{
    void exact(double x, double& ex){
        ex = x*x;
    }

    void Euler_explicit(double h, double n, double& y_in, double & y_out){
        for(int i = 0; i < n; i++){
            double a = 2 * h * (sqrt(y_in));
            y_out = y_in + a;
            y_in = y_out;
        }
    }

    void Euler_implicit(double h, double n, double& y_in, double & y_out){
        for(int i = 0; i < n; i++){
            double a = sqrt( (h*h) + y_in );
            a = h + a;
            y_out = a*a;
            y_in = y_out;
        }
    }
    void run(){
        double ex;
        exact(1, ex);
        int n = 10;
        double y_in, y_out;
        cout << endl << "dy/dx = 2*sqrt(y)" << endl;
        cout << "Explicit\n";
        for(int i = 0; i<3; i++){
            y_in = y_out = 0;
            double h = 1.0/n;
            Euler_explicit(h, n, y_in, y_out);
            cout << "n = " << n << "  " << y_out << endl;
            y_in = y_out = 0;
            n*=10;  
        }
        cout << "Implicit\n";
        n = 10;
        for(int i = 0; i<3; i++){
            y_in = y_out = 0;
            double h = 1.0/n;
            Euler_implicit(h, n, y_in, y_out);
            cout << "n = " << n << "  " << fixed << setprecision(5) << y_out << " " << fixed <<  setprecision(2) << n * (ex - y_out)<< endl;
            y_in = y_out = 0;
            n*=10;  
        }
    }
}
