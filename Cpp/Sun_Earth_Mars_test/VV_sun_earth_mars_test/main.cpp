#include "include/armadillo"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
using namespace std;
using namespace arma;

// generating mass with gaussian distribution arounf 10 with standard deviation 1
void gaussian_mass_generator(vec (&mass), int number_of_particles)
{
  srand(time(NULL));
  for (int i = 0; i < number_of_particles; i++)
  {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
    do{
      v1 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      v2 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    mass(i) = v2*fac;
    mass(i) += 10;
  }
} // end function for gaussian deviates


void Derivative(mat r, mat v, vec m, mat (&drdt), mat (&dvdt), double G, int number_of_particles){
    double acc_x = 0, acc_y = 0, acc_z = 0;
    double distance_squared = 0;
    double force_between_particles = 0;
    for (int i=0; i<number_of_particles; i++)
    {
        force_between_particles = 0;
        acc_x = 0;
        acc_y = 0;
        acc_z = 0;
        for (int j=0; j<number_of_particles; j++)
        {
            if (j!=i)
            {
                distance_squared = 0;
                for (int k=0; k<3; k++)
                {
                    distance_squared += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
                }
                force_between_particles = m(j)*pow(distance_squared,-1.5);
                acc_x += G*force_between_particles*(r(j,0)-r(i,0));
                acc_y += G*force_between_particles*(r(j,1)-r(i,1));
                acc_z += G*force_between_particles*(r(j,2)-r(i,2));
            }
        }
        dvdt(i,0) = acc_x;
        dvdt(i,1) = acc_y;
        dvdt(i,2) = acc_z;
    }

    for (int i = 0; i<number_of_particles; i++)
    {
        drdt(i,0) = v(i,0);
        drdt(i,1) = v(i,1);
        drdt(i,2) = v(i,2);
    }
}

void making_ks(mat &pos_temp, mat &vel_temp, mat r, mat v, mat drdt, mat dvdt, double time_step, mat &k, double number, int number_of_particles){
    for (int i = 0; i< number_of_particles; i++)
    {
        for (int j=0; j<3; j++)
        {
            k(i,j) = time_step*drdt(i,j);
            k(i,j+3) = time_step*dvdt(i,j);
            pos_temp(i,j) = r(i,j) + k(i,j)/number;
            vel_temp(i,j) = v(i,j) + k(i,j+3)/number;
        }
    }
}


int main()
{
    ofstream myfile ("VV_sun_earth_mars_test.txt");
            if (myfile.is_open())
            {
    int number_of_particles = 3;
    int n = 3;
    double G = 2.96e-4;     //Grav const in the units of AU^3 / ( days^3 * mass_sun )
    double t0 = 0.0;
    double t_final = 20*365;     //time in the unit of days
    int number_of_time_step = 20*365;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    mat r(number_of_particles,3);
    mat v(number_of_particles,3);
    mat r_initial(number_of_particles,3);
    mat v_initial(number_of_particles,3);
    //position of Sun
    r(0,0) = 1.0;
    r(0,1) = 1.0;
    r(0,2) = 1.0;
    //position of Earth
    r(1,0) = 2.0;
    r(1,1) = 1.0;
    r(1,2) = 1.0;
    //position of Mars
    r(2,0) = -0.5;
    r(2,1) = 1.0;
    r(2,2) = 1.0;
    //velocity of Sun
    v(0,0) = 0.0;
    v(0,1) = 0.0;
    v(0,2) = 0.0;
    //velocity of Earth
    v(1,0) = 0.0;
    v(1,1) = 0.017;
    v(1,2) = 0.0;
    //velocity of Mars
    v(2,0) = 0.0;
    v(2,1) = 0.014;
    v(2,2) = 0.0;
    vec m(number_of_particles);
    m(0) = 1.0;
    m(1) = 3.0e-6;
    m(2) = 3.2e-7;
    myfile << "Sun_x" << setw(20) << "Sun_y" << setw(20) << "Sun_z" << setw(20) << "Earth_x" << setw(20) << "Earth_y" << setw(20) << "Earth_z" << setw(20) << "Mars_x" << setw(20) << "Mars_y" << setw(20) << "Mars_z" << endl;
    myfile << r(0,0) << setw(20) << r(0,1) << setw(20) << r(0,2) << setw(20) << r(1,0) << setw(20) << r(1,1) << setw(20) << r(1,2) << setw(20) << r(2,0) << setw(20) << r(2,1) << setw(20) << r(2,2) << endl;

    mat v_partly(number_of_particles,3);
    v_partly.zeros();
    mat drdt(number_of_particles,3);
    drdt.zeros();
    mat dvdt(number_of_particles,3);
    dvdt.zeros();
    while(time<=t_final){
         Derivative(r,v,m,drdt,dvdt,G,number_of_particles);

        for(int i=0; i<number_of_particles ; i++){
            for(int j=0; j<3 ; j++){
               r(i,j) = r(i,j)+dt*drdt(i,j)+0.5*dt*dt*dvdt(i,j);
               v_partly(i,j) = drdt(i,j)+0.5*dt*dvdt(i,j);
               dvdt(i,j) = v_partly(i,j);
            }
        }
        Derivative(r,v,m,drdt,dvdt,G,number_of_particles);
        for(int i=0; i<number_of_particles ; i++){
            for(int j=0; j<3 ; j++){
        v(i,j) = v_partly(i,j)+0.5*dt*dvdt(i,j);
            }
        }
        myfile << r(0,0) << setw(20) << r(0,1) << setw(20) << r(0,2) << setw(20) << r(1,0) << setw(20) << r(1,1) << setw(20) << r(1,2) << setw(20) << r(2,0) << setw(20) << r(2,1) << setw(20) << r(2,2) << endl;

    time += dt;
    }
    }
}

