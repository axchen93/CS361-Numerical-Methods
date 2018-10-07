#include<iostream>
#include<cmath>
using namespace std;

double sum_prod_ab1(int n, const double a[], const double b[]);
double sum_prod_ab2(int n, const double a[], const double b[]);
double sum_prod1(int n, const double a[], const double b[], const double c[]);
double sum_prod11(int n, const double a[], const double b[], const double c[]);
double sum_prod2(int n, const double a[], const double b[], const double c[]);
double sum_prod22(int n, const double a[], const double b[], const double c[]);


int main(){
    double a[] = {1.1, 2.15, 3.14519, -1.369, 0.002, 9.632, 158.25584, 15.66854, 254.22541, -996.36};
    double b[] = {2.366, 2.3355, -66.65, 100.1236, 1.3362, 6.36, 1.002, -.366, 6.332, 2.36};
    double c[] = {1, 2, 3, 4 , 5, 6, 7, 8 , 9, 10};
    int n =10;
    cout << sum_prod_ab1(n, a, b) << endl;
    cout << sum_prod_ab2(n, a, b) << endl;
    cout << sum_prod1(n, c, c, c) << endl;
    cout << sum_prod2(n, a, c, c) << endl;
}

double sum_prod2(int n, const double a[], const double b[], const double c[]){
    double totalSum = 0.0, currentC = 0.0, partSum = 0.0;
    //Complexity O(N)
    for(int i =0; i<n; i++){
        if(i%2 == 0){
            currentC += c[i];
            partSum += (b[i]*currentC);
            totalSum += (a[i]*partSum);
        }
        else{
            currentC -= c[i];
            partSum += (b[i]*currentC);
            totalSum -= (a[i]*partSum);
        }
    }
    return totalSum;
} 

double sum_prod1(int n, const double a[], const double b[], const double c[]){
    //Complexity O(N)
    double sumA = 0.0, sumB = 0.0, sumC = 0.0;
    for(int i = 0; i<n; i++){
        sumB += b[i];
        if(i%2 == 0){
            sumA += a[i];
            sumC += b[i];
        }else{
            sumA -= a[i];
            sumC -= c[i];
        }
    }
    return sumA * sumB * sumC;
}

double sum_prod_ab1(int n, const double a[], const double b[]){
    //Complexity O(N)
    double sumA = 0.0, sumB = 0.0;
    for(int i = 0; i<n; i++){
        sumA += a[i];
        sumB += b[i];
    }
    return sumA * sumB;
}

double sum_prod_ab2(int n, const double a[], const double b[]){
    //Complexity O(N)
    double totalSum = 0.0, currentB = 0.0;
    for(int i = 0; i<n; i++){
        currentB += b[i];
        totalSum += (a[i]*currentB);
    }
    return totalSum;
}

