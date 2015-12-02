#ifndef STARS_H
#define STARS_H
class stars{
public: // Setting up the public that is accessible to all of the classes
    double position[3];
    double velocity[3];
    double mass;
    stars(double mas, double x, double y, double z, double vx, double vy, double vz);
    stars();
};

#endif // STARS_H


