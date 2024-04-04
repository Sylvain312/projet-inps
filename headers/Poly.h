#ifndef POLY_H
#define POLY_H

#include <armadillo>

/** 
 * \brief Class of matrices that contain the solutions.
 * \param row number of row
 * \param col number of column
*/
class Poly 
{
    public:
        //attributs
        //arma::vec z; //Vector of the values (points to evaluates)

        arma::mat psi_matrix; // the matrix of the psi solutions evaluated on a mesh
        arma::mat H; //Hermite part of the solution
        arma::cube L; //Laguerre part of the solution
        arma::vec E; //exponentionnal part of the solution
        arma::colvec F; //part of the solution dependent on n except the hermit polynomial
        arma::mat scal; // the matrix of the scalar product of computed solutions


        //Constructor
        Poly();


        //Methods
/**
 * @brief Function that evaluates the Hermit matrix -> We want to evaluate the matrix at points sqrt(m*w/hbarre)*z
 *
 * @assigns the value of the Hermit matrix to H (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the hermit matix : H .
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
        void calcHermite(int, arma::vec);

/**
 * @brief Function which return the values of H for a specific value of n
 *
 * @assigns the value of the Hermit matrix to H (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the hermit matrix : H .
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
		arma::vec hermite(int n);

/**
 * @brief Function that evaluates the Hermit matrix -> We want to evaluate the matrix at points sqrt(m*w/hbarre)*z
 *
 * @assigns the value of the Hermit matrix to H (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the hermit matix : H .
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
        void calcLaguerre(int, int, arma::vec);

/**
 * @brief Function which return the values of H for a specific value of n
 *
 * @assigns the value of the Hermit matrix to H (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the hermit matrix : H .
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
		arma::vec laguerre(int, int);
};

#endif