#include<iostream>
#include <string>
using namespace std;

int Dyck_path(int n, int k);
void Dyck_path(int i, int j, const int n, const int k, int &numPaths, string path);
bool isOutbound(int i, int j, int n, int k);

int main(){
    int n = 7;
    for(int i = 1; i<=4; i++){
        for(int j =3; j<=n; j++){
            cout << i << " " << j << " : " << Dyck_path(j, i) << endl;
        }
        n--;
    }
}

int Dyck_path(int n, int k){
    int numPath = 0;
    Dyck_path(0,0,n,k,numPath, "");
    return numPath;
}

void Dyck_path(int i, int j, const int n, const int k, int &numPaths, string path){
    if( i==n && j==(k*n)){
        numPaths++;
        path = path + "(" + to_string(i) + ", " + to_string(j) +") ";
        cout << path << endl;
        return;
    }
    path = path + "(" + to_string(i) + ", " + to_string(j) +") ";

    if(!isOutbound(i, j+1, n, k)){
        Dyck_path(i, j+1, n, k, numPaths, path);
    }
    if(!isOutbound(i+1, j, n, k)){
        Dyck_path(i+1, j, n, k, numPaths, path);
    }
    return;
}

bool isOutbound(int i, int j, int n, int k){
    if(i>n) return true;
    if(j>k*n) return true;
    if(j>k*i) return true;
    return false;
}