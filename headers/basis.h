#ifndef BASIS_H
#define BASIS_H

#include <armadillo>

class Basis
{
    public:
/**
 * @brief values for each quantum number
*/
        int mMax;
        arma::ivec nMax;
        arma::imat n_zMax;
        double bz;
        double br;

        // Constructor
        Basis();
        Basis(double br, double bz, int N, double Q);

        // Methods

/**
 * @brief Function which is used to calculate mMax and n_zMax values of the quantum number
 * 
 * @param N
 * @param Q
 * @param i
 * 
 * @return a double
*/
        double v(int N, double Q, int i);

/**
 * @brief The function that returns the r part of the psi function.
 * 
 * @param r a vector of radius values
 * @param m a quantum number
 * @param n a quantum number
 * 
 * @return a vector containing the r part of the psi function.
*/
        arma::vec rPart(arma::vec r, int m, int n);

/**
 * @brief The function that returns the z part of the psi function.
 * 
 * @param z a vector of z-axis values
 * @param nz a quantum number
 * 
 * @return a vector containing the z part of the psi function.
*/
        arma::vec zPart(arma::vec z, int nz);

/**
 * @brief Function which calculates the psi-function on all couple(z, r) given in the parameters zVals and rVals
 * 
 * @param m a quantum number
 * @param n a quantum number
 * @param n_z a quantum number
 * @param zVals the set of values on the z-axis
 * @param rVals the set of values on the radius
 * 
 * @return a matrix containing the psi function for every (r,z) couple.
*/
        arma::mat psi(int m, int n, int n_z, arma::vec zVals, arma::vec rVals);
};

#endif