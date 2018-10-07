#include<iostream>
#include<cmath>
using namespace std;

typedef void (*FUNC)(double x, double& f, double& fprime);
double func1(double x);
double func2(double x);
void funcIntersect(double x, double& f, double& fprime);
void func1Root(double x, double& f, double& fprime);
void func2Root(double x, double& f, double& fprime);
void intersectNewtonRaphson(FUNC func, double x0, double& answer);
void midPointIntegration(const double root1,const double root2,  const double intersection, const int n, double& area);

int main(){
    double intersection;
    intersectNewtonRaphson(funcIntersect, 0, intersection);

    double root1, root2;
    intersectNewtonRaphson(func1Root, 0, root1);
    intersectNewtonRaphson(func2Root, 0, root2);
    
    double area = 0;
    double previousArea = 0;
    double tol = 1.0e-5;
    for(int i = 1; i<20; i++){
        int n = pow(2,i);
        midPointIntegration(root1, root2, intersection, n, area);
        cout << "n = " << n << " " << area <<endl;
        if(abs(area-previousArea) <= tol) break;
        previousArea = area;
    }
    
}

double func1(double x){
    const double p = 0.2302;
    double a = ( 1 - (x-p)*(x-p) ) / 2.0;
    return sqrt(a);
}

double func2(double x){
    const double q = 0.7199;
    double a = ( 1 - (x+q)*(x+q) ) / 4.0;
    return sqrt(a);
}

void funcIntersect(double x, double& f, double& fprime){
    const double p = 0.2302, q = 0.7199;
    f = 2 - 4*(x-p)*(x-p) + 2*(x+q)*(x+q);
    fprime = -8*(x-p) + 4*(x+q);
}

void func1Root(double x, double& f, double& fprime){
    const double p = 0.2302;
    f = (x - p)*(x-p) -1.0;
    fprime = 2*(x-p);
}

void func2Root(double x, double& f, double& fprime){
    const double q = 0.7199;
    f = (x + q)*(x+q) -1.0;
    fprime = 2*(x+q);
}


void intersectNewtonRaphson(FUNC func, double x0, double& answer){
    const double tol = 1.0e-5;
    const int max_iteration = 100;
    double f = 0, fprime = 0;

    double diff;
    double x1 = x0;
    double x2 = 0;
    for(int i = 0; i<max_iteration; i++){
        func(x1, f, fprime);
        x2 = x1 - (f/fprime);
        diff = x2-x1;
        if(abs(diff) <= tol){
            answer = x2;
            return;
        }
        x1 = x2;
    }
    answer = 0;
    return;
}

void midPointIntegration(const double root1, const double root2,  const double intersection, const int n, double& area){
    double h1 = (intersection - root1)/n;
    double h2 = (root2 - intersection)/n;

    double sum1, sum2;
    sum1 = 0.5 * (func1(intersection));
    sum2 = 0.5 * (func2(intersection));

    for(int i = 1; i<n; i++){
        double x_1 = root1 + i * h1; 
        sum1+= func1(x_1);

        double x_2 = intersection + i * h2; 
        sum2+= func2(x_2);
    }

    double integral1 = sum1*h1;
    double integral2 = sum2*h2;

    area = (integral1 + integral2) *2.0;
}
