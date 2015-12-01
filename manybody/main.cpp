#include <iostream>
#include <armadillo>
#include <star_cluster.h>
#include <stars.h>
#include <cmath>

using namespace arma;
using namespace std;

int main()
{
    star_cluster mysystem;
    stars Sun(1,0,0,0,0,0,0);

    mysystem.add(Sun);

}

