#include <math.h>
#include "../headers/Poly.h"

/** 
 * \brief Class of matrices that contain the solutions.
*/
Poly::Poly()
{
    arma::mat a;
    arma::vec b;
    arma::colvec c;
	arma::cube d;
    psi_matrix = a;
	H = a;
	L = d;
    E = b;
    F = c;
    scal = a;
}

/**
 * @brief Function that evaluates the Hermit matrix -> We want to evaluate the matrix at points sqrt(m*w/hbarre)*z
 *
 * @assigns the value of the Hermit matrix to H (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the hermit matrix : H .
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
void Poly::calcHermite(int n, arma::vec z)
{
    H.arma::mat::ones(z.n_elem, n+1);

    if (n==0) return;

    //H_1(z)
    H.col(1) = 2 * H.col(1) % z;

    //H_{n+1}(z)
    for (int k = 1; k < n - 1; k++) 
    {
        H.col(k+1) = 2 * z % H.col(k) - 2 * k * H.col(k-1);
    }
}

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
arma::vec Poly::hermite(int n)
{
    return H.col(n);
}

/**
 * @brief Function that evaluates the Laguerre cube
 *
 * @assigns the value of the Laguerre cube to L (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the Laguerre matrix : L.
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
void Poly::calcLaguerre(int m, int n, arma::vec z) 
{
    
    // We set L_m_0(z) to 1
    L.arma::cube::ones(z.n_rows, n + 1, m + 1);

	for (int i = 0; i < m; i++)
	{
		//We set L_i_1(z) to 1 + i - z
        L.slice(i).col(1) = 1 + i - z;

		//L_m_{n+1}(z)
        for (int k = 1; k < n - 1; k++)
        {
            L.slice(i).col(k+1) = (2 + (i - 1 - z) / (double)(k + 1)) % L.slice(i).col(k) - (1 + (i - 1) / (double)(k + 1)) * L.slice(i).col(k - 1);
        }
	}
}

/**
 * @brief Function which return the values of L for a specific value of n and m
 *
 * @assigns the value of the Laguerre cube to L (attribute of the class)
 * @attention requires that the constructor has been called to set the size of the Laguerre cube : L.
 * 
 * @param z the vector with the points where we evaluate the function
 * 
 * @return nothing
 */
arma::vec Poly::laguerre(int m, int n)
{
    return L.slice(m).col(n);
}