/**
 * @file main.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //for the use of input and output stream

#include "../headers/Poly.h"
#include "../headers/basis.h"


//to modify files
#include <fstream> 

using namespace std;

arma::mat rho;
//     br = 1.935801664793151, bz = 2.829683956491218, N = 14, Q = 1.3
Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);

/**
 * @brief Naive algorithme given in the project presentation
 * 
 * @param r r- axis values
 * @param z z- axis values
 * 
 * @return the nuclear local density matrix
*/
arma::mat naive_algo(arma::vec r, arma::vec z)
{
    int a = 0, b = 0;
    arma::mat result = arma::zeros(r.n_rows, z.n_rows);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat psiA = basis.psi(m, n, n_z, z, r);
                            arma::mat psiB = basis.psi(mp, np, n_zp, z, r);

                            int i;
                            i = 0;
                            for (int m_loop = 0; m_loop < basis.mMax; m_loop++)
                            {
                                for (int n_loop = 0; n_loop < basis.nMax(m_loop); n_loop++)
                                {
                                    for (int n_z_loop = 0; n_z_loop < basis.n_zMax(m_loop, n_loop); n_z_loop++)
                                    {
                                        if (m_loop == m && n_loop == n && n_z_loop == n_z)
                                        {
                                            a = i;
                                        }
                                        i++;
                                    }
                                }
                            }

                            i = 0;
                            for (int m_loop = 0; m_loop < basis.mMax; m_loop++)
                            {
                                for (int n_loop = 0; n_loop < basis.nMax(m_loop); n_loop++)
                                {
                                    for (int n_z_loop = 0; n_z_loop < basis.n_zMax(m_loop, n_loop); n_z_loop++)
                                    {
                                        if (m_loop == mp && n_loop == np && n_z_loop == n_zp)
                                        {
                                            b = i;
                                        }
                                        i++;
                                    }
                                }
                            }

                            result += psiA % psiB * rho(a, b);
                        }
                    }
                }
            }
        }
    }
    return result;
}

arma::icube init_rho_indices()
{ 
    arma::icube indices = arma::icube(basis.mMax, basis.nMax(0), basis.n_zMax(0, 0));
    int i = 0;
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                indices(m, n, n_z) = i;
                i++;
            }
        }
    }

    return indices;
}

arma::icube rho_indice = init_rho_indices();

// first acceleration -> stored all of the indices
arma::mat accel1_algo(arma::vec r, arma::vec z)
{
    int a, b;

    arma::mat result = arma::zeros(r.n_rows, z.n_rows);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                

                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat psiA = basis.psi(m, n, n_z, z, r);
                            arma::mat psiB = basis.psi(mp, np, n_zp, z, r);
                            a = rho_indice(m, n, n_z);
                            b = rho_indice(mp, np, n_zp);

                            result += psiA % psiB * rho(a, b);
                        }
                    }
                }
            }
        }
    }

    return result;
}

// Second acceleration -> moved anything a related as high as possible
arma::mat accel2_algo(arma::vec r, arma::vec z)
{
    int a, b;

    arma::mat result = arma::zeros(r.n_rows, z.n_rows);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat psiA = basis.psi(m, n, n_z, z, r);
                a = rho_indice(m, n, n_z);

                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat psiB = basis.psi(mp, np, n_zp, z, r);
                            b = rho_indice(mp, np, n_zp);

                            result += psiA % psiB * rho(a, b);
                        }
                    }
                }
            }
        }
    }

    return result;
}


// Thrid acceleration -> mp is the same as m so we can remove a loop
arma::mat accel3_algo(arma::vec r, arma::vec z)
{

    int a, b;

    arma::mat result = arma::zeros(r.n_rows, z.n_rows);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat psiA = basis.psi(m, n, n_z, z, r);
                a = rho_indice(m, n, n_z);

                for (int np = 0; np < basis.nMax(m); np++)
                {
                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++)
                    {
                        arma::mat psiB = basis.psi(m, np, n_zp, z, r);
                        b = rho_indice(m, np, n_zp);
                        result += psiA % psiB * rho(a, b);
                    }
                }
            }
        }
    }

    return result;
}

// Fourth acceleration -> removed redudant call to basis.nMax
arma::mat accel4_algo(arma::vec r, arma::vec z)
{
    int a, b;
    arma::mat result = arma::zeros(r.n_rows, z.n_rows);
    for (int m = 0; m < basis.mMax; m++)
    {
        int nmax = basis.nMax(m);
        for (int n = 0; n < nmax; n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat psiA = basis.psi(m, n, n_z, z, r);
                a = rho_indice(m, n, n_z);

                for (int np = 0; np < nmax; np++)
                {
                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++)
                    {
                        arma::mat psiB = basis.psi(m, np, n_zp, z, r);
                        b = rho_indice(m, np, n_zp);
                        result += psiA % psiB * rho(a, b);
                    }
                }
            }
        }
    }

    return result;
}


std::string cubeToDf3(const arma::cube &m)
{
  std::stringstream ss(std::stringstream::out | std::stringstream::binary);
  int nx = m.n_rows;
  int nz = m.n_cols;
  int ny = m.n_slices;
  ss.put(nx >> 8);
  ss.put(nx & 0xff);
  ss.put(ny >> 8);
  ss.put(ny & 0xff);
  ss.put(nz >> 8);
  ss.put(nz & 0xff);
  double theMin = 0.0;
  double theMax = m.max();
  for (uint k = 0; k < m.n_cols; k++)
  {
    for (uint j = 0; j < m.n_slices; j++)
    {
      for (uint i = 0; i < m.n_rows; i++)
      {
        uint v = 255 * (fabs(m(i, k, j)) - theMin) / (theMax - theMin);
        ss.put(v);
      }
    }
  }
  return ss.str();
}


int main(int argc, char *argv[])
{
    rho.load("rho.arma", arma::arma_ascii); 

    if (argc != 2)
    {
        printf("Error : Usage : \n - ./bin/exec -a for acceleration computations \n - ./bin/exec -r for density computation\n");
        return 1;
    }

    if (!strcmp(argv[1],"-a"))
    {
        struct timespec tstart, tend;
        double naiveTime, accelTime1, accelTime2, accelTime3, accelTime4;
        double accel1, accel2, accel3, accel4;
        
        arma::vec r_accel = {3.1, 2.3, 1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
        arma::vec z_accel = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};

        // Acceleration Calculations
        
        //Naive Version
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat naive_result = naive_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        naiveTime = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        printf("Time for the naive implementation %fs \n", naiveTime);

        //Accelerated 1 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel1_result = accel1_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime1 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel1 = naiveTime / accelTime1;
        printf("Factor of acceleration for the first optimisation %f \n", accel1);

        //Accelerated 2 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel2_result = accel2_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime2 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel2 = naiveTime / accelTime2;
        printf("Factor of acceleration for the second optimisation %f \n", accel2);

        //Accelerated 3 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel3_result = accel3_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime3 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel3 = naiveTime / accelTime3;
        printf("Factor of acceleration for the third optimisation %f \n", accel3);

        //Accelerated 4 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel4_result = accel4_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime4 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel4 = naiveTime / accelTime4;
        printf("Factor of acceleration for the fourth optimisation %f \n", accel4);
        return 0;
    }
    
    
    if (!strcmp(argv[1],"-r"))
    {
        // 2D plot and csv generation
        arma::vec r_calc = arma::vec(251).zeros();
        arma::vec z_calc = arma::vec(501).zeros();
        arma::vec x_calc = arma::vec(501);

        double h = 0.1;
        for (int i = 1; i < 251; i++)
        {
            r_calc(i) = r_calc(i - 1) + h;
            z_calc(250 + i) = z_calc(250 + i - 1) + h;
            z_calc(250 - i) = z_calc(250 - i + 1) - h;
        }

        arma::mat result = accel3_algo(r_calc, z_calc);

        arma::mat result_2D = arma::join_cols(arma::reverse(result.submat(1, 0, 250, 500)), result);


        double min_2D = 0.0;
        double max_2D = result_2D.max();

        result_2D.transform(([min_2D, max_2D](double val)
                            { return 255 * (val - min_2D) / (max_2D - min_2D); }));


        result_2D.save("result.csv", arma::csv_ascii);

        // 3D plot and df3 generation
        arma::vec x_3D_plot = arma::vec(32);
        arma::vec y_3D_plot = arma::vec(32);
        arma::vec z_3D_plot = arma::vec(64);

        x_3D_plot(0) = -1 * 10;
        y_3D_plot(0) = -1 * 10;
        z_3D_plot(0) = -2 * 10;
        z_3D_plot(63) = 2 * 10;
        double h_xy = (2. / (32. - 1.)) * 10;
        double h_z = (4. / (64. - 1.)) * 10;
        for (int i = 1; i < 32; i++)
        {
            x_3D_plot(i) = x_3D_plot(i - 1) + h_xy;
            y_3D_plot(i) = y_3D_plot(i - 1) + h_xy;
            z_3D_plot(2 * i - 1) = z_3D_plot(2 * i - 2) + h_z;
            z_3D_plot(2 * i) = z_3D_plot(2 * i - 1) + h_z;
        }

        arma::mat y_squared_3D = y_3D_plot;
        y_squared_3D.transform([](double val)
                            { return (val * val); });

        arma::cube m = arma::cube(32, 64, 32);
        arma::vec r_3D_plot = arma::vec(32);

        for (int i = 0; i < 32; i++)
        {
            r_3D_plot = x_3D_plot(i) * x_3D_plot(i) + y_squared_3D;

            r_3D_plot.transform(([](double val)
                                { return (pow(val, 0.5)); }));

            m.slice(i) = accel3_algo(r_3D_plot, z_3D_plot);
        }


        //To plot the graph: povray +A0.0001 -W800 -H600 +P +Q11 visu.pov 
        std::string input = cubeToDf3(m);
        std::ofstream out("sphere.df3");
        out << input;
        out.close();
        printf("Computations Done\n");
        system("povray visu.pov");
        return 0;
    }

    if (!strcmp(argv[1],"-ra") || !strcmp(argv[1],"-ar"))
    {
        struct timespec tstart, tend;
        double naiveTime, accelTime1, accelTime2, accelTime3;
        double accel1, accel2, accel3;
        
        arma::vec r_accel = {3.1, 2.3, 1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
        arma::vec z_accel = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};

        // Acceleration Calculations
        
        //Naive Version
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat naive_result = naive_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        naiveTime = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        printf("Time for the naive implementation %fs \n", naiveTime);

        //Accelerated 1 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel1_result = accel1_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime1 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel1 = naiveTime / accelTime1;
        printf("Factor of acceleration for the first optimisation %f \n", accel1);

        //Accelerated 2 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel2_result = accel2_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime2 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel2 = naiveTime / accelTime2;
        printf("Factor of acceleration for the second optimisation %f \n", accel2);

        //Accelerated 3 Time
        clock_gettime(CLOCK_REALTIME, &tstart);
        arma::mat accel3_result = accel3_algo(r_accel, z_accel);
        clock_gettime(CLOCK_REALTIME, &tend);
        accelTime3 = (1000. * (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1000000.) / 1000.;
        accel3 = naiveTime / accelTime3;
        printf("Factor of acceleration for the third optimisation %f \n", accel3);
        
        // 2D plot and csv generation
        arma::vec r_calc = arma::vec(251).zeros();
        arma::vec z_calc = arma::vec(501).zeros();
        arma::vec x_calc = arma::vec(501);

        double h = 0.1;
        for (int i = 1; i < 251; i++)
        {
            r_calc(i) = r_calc(i - 1) + h;
            z_calc(250 + i) = z_calc(250 + i - 1) + h;
            z_calc(250 - i) = z_calc(250 - i + 1) - h;
        }

        arma::mat result = accel4_algo(r_calc, z_calc);
        arma::mat result_2D = arma::join_cols(arma::reverse(result.submat(1, 0, 250, 500)), result);
        double min_2D = 0.0;
        double max_2D = result_2D.max();
        result_2D.transform(([min_2D, max_2D](double val)
                            { return 255 * (val - min_2D) / (max_2D - min_2D); }));
        result_2D.save("result.csv", arma::csv_ascii);
        system("python3 plot.py");

        // 3D plot and df3 generation
        arma::vec x_3D_plot = arma::vec(32);
        arma::vec y_3D_plot = arma::vec(32);
        arma::vec z_3D_plot = arma::vec(64);

        x_3D_plot(0) = -1 * 10;
        y_3D_plot(0) = -1 * 10;
        z_3D_plot(0) = -2 * 10;
        z_3D_plot(63) = 2 * 10;
        double h_xy = (2. / (32. - 1.)) * 10;
        double h_z = (4. / (64. - 1.)) * 10;
        for (int i = 1; i < 32; i++)
        {
            x_3D_plot(i) = x_3D_plot(i - 1) + h_xy;
            y_3D_plot(i) = y_3D_plot(i - 1) + h_xy;
            z_3D_plot(2 * i - 1) = z_3D_plot(2 * i - 2) + h_z;
            z_3D_plot(2 * i) = z_3D_plot(2 * i - 1) + h_z;
        }

        arma::mat y_squared_3D = y_3D_plot;
        y_squared_3D.transform([](double val)
                            { return (val * val); });

        arma::cube m = arma::cube(32, 64, 32);
        arma::vec r_3D_plot = arma::vec(32);

        for (int i = 0; i < 32; i++)
        {
            r_3D_plot = x_3D_plot(i) * x_3D_plot(i) + y_squared_3D;

            r_3D_plot.transform(([](double val)
                                { return (pow(val, 0.5)); }));

            m.slice(i) = accel4_algo(r_3D_plot, z_3D_plot);
        }
        
        std::string input = cubeToDf3(m);
        std::ofstream out("sphere.df3");
        out << input;
        out.close();
        printf("Computations Done\n");
        system("povray visu.pov");
        return 0;
    }

    printf("Error : Invalid Option : \n - ./bin/exec -a for acceleration computations \n - ./bin/exec -r for density computation\n");
    return 1;
}