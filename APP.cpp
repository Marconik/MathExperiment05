#include <bits/stdc++.h>
#include "LA.h"

int main()
{
    std::freopen("result.txt","w",stdout);
    int n=10;
    std::vector<double> b(n,1);
    LinearEquation M(n);
    M.Hilbert();
    M.b=b;
    printf("H:\n");
    M.displayA();
    printf("b:\n");
    M.displayVec(b);

    std::vector<double> x(n),y(n);
    printf("Gauss-Seidel:\n");
    M.GSIter(x);
    M.displayVec(x);
    printf("Gradient Descent:\n");
    M.GD(y);
    M.displayVec(y);
    return 0;
}