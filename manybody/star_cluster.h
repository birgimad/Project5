#ifndef STAR_CLUSTER_H
#define STAR_CLUSTER_H
#include <armadillo>
#include <vector>
#include "stars.h"

using std::vector;
using namespace arma;


class star_cluster
{
public:
    int number_stars = 0;
    star_cluster();
    vector<stars> all_stars;
    void add(stars n);
};

#endif // STAR_CLUSTER_H
