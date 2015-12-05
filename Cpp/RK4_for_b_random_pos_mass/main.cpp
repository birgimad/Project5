#include <iostream>
#include "include/armadillo"
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

int main()
{
    int number_of_particles = 100000;
    vec m(number_of_particles);
    m.zeros();
    gaussian_mass_generator(m,number_of_particles);
    ofstream myfile ("random_mass_generator_test.txt");
        if (myfile.is_open())
        {
        myfile << m << endl;
        }
    cout << "hello" << endl;
    return 0;
}

