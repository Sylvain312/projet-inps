/** 
* \file basis.cpp
* This file contains the declaration of the basis class
*/

#include "../headers/basis.h"
#include <math.h>
#include "../headers/Poly.h"
#include <iostream>
using namespace std;

Basis::Basis()
{

}

double Basis::v(int N, double Q, int i) 
{
    return (N + 2) * pow(Q, 2./3.) + (1./2.) - i * Q;
}


Basis::Basis(double br, double bz, int N, double Q)
{
    int i = 0;
    float iMax = ((N + 2) * pow(Q, 2./3.) - 1./2.) / Q;
    while (i < iMax - 1)
    {
        i++;
    }
    mMax = i;

    nMax = arma::ivec(mMax);
    for (int m = 0; m < mMax; m++)
    {
        nMax[m] = (1./2.) * (mMax - m - 1) + 1;
    }

    n_zMax = arma::imat(mMax, nMax[0]);
    //nMax[0] is the maximum value that can reach n.
    for (int m = 0; m < mMax; m++)
    {
        for (int n = 0; n < nMax[0]; n++)
        {
            float value_of_v = v(N, Q, m + 2 * n + 1);
            if (value_of_v > 0)
            {
                n_zMax(m, n) = value_of_v;
            }
        }
    }
    this->bz = bz;
    this->br = br;
}

arma::vec Basis::zPart(arma::vec z, int nz)
{
    double pi = arma::datum::pi;

    double calc = pow(pi, -0.25);

    Poly poly = Poly();
    poly.calcHermite(nz+1,z/bz);

    for (int i=1; i<=nz; i++){
        calc *= pow(2 * i, -0.5);
    }

    return (1.0 / sqrt(bz)) * calc * arma::exp((-1.0 / (2 * bz * bz))* z % z) % poly.hermite(nz);
}

arma::vec Basis::rPart(arma::vec r, int m, int n)
{
    arma::vec res;
    double pi = arma::datum::pi;

    double calc = 1.0 / (br*sqrt(pi));
    double fact = 1.0;

    for (int i=n+1; i <= n + abs(m); i++){
        fact *= sqrt(1.0/i);
    }

    res = (-1.0 / (2* br * br)) * r % r;
    res = arma::exp(res);

    Poly poly = Poly();
    poly.calcLaguerre(m+1,n+1,r%r/(br*br));

    return calc * fact * res % arma::pow((1.0/br)*r,m) % poly.laguerre(m,n);
}

arma::mat Basis::psi(int m, int n, int n_z, arma::vec zVals, arma::vec rVals)
{
    return Basis::rPart(rVals, m, n) * Basis::zPart(zVals, n_z).t();
}