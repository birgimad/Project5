#include <iostream>
#include "include/armadillo"
using namespace std;
using namespace arma;

void Derivative(double (&r)[3], double (&v)[3], double (&drdt)[3], double (&dvdt)[3], double G, double mass){
    //int n = 3;       // use evt: sizeof(drdt)/sizeof(drdt[0])
    drdt[0] = v[0];
    drdt[1] = v[1];
    drdt[2] = v[2];

    double distance_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double newtonian_force = -G*mass/pow(distance_squared,1.5);
    dvdt[0] = newtonian_force*r[0];
    dvdt[1] = newtonian_force*r[1];
    dvdt[2] = newtonian_force*r[2];
}
/*
void updating_dummies(double dt, double drdt, double dvdt, double (&r_dummy)[3], double (&v_dummy)[3], int number)
{
    double kr[n];
    double kv[n];
    for (int i = 0; i<n; i++){
        kr[i] = dt*drdt[i];
        kv[i] = dt*dvdt[i];
        r_dummy[i] = r[i] +kr[i]/number;
        v_dummy[i] = v[i] +kv[i]/number;
    }
} */


int main()
{
    int n=3;

    double G = 2.96e-4;
    double mass_earth = 3e-6;
    double mass = 1;
    double t0 = 0.0;
    double t_final = 365.0;
    int number_of_time_step = 100000;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;

    double r[3];
    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;
    double v[3];
    v[0] = 0;
    v[1] = 0.017;
    v[2] = 0;
    double drdt[3];
    drdt[0] = 0;
    drdt[1] = 0;
    drdt[2] = 0;
    double dvdt[3];
    drdt[0] = 0;
    drdt[1] = 0;
    drdt[2] = 0;
    cout << "initial position:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r[i] << endl;
    }
    cout << "initial velocity:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v[i] << endl;
    }
    double energy_initial = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    // energy calculated in the units of J/mass_earth (orbiting object)
    cout << "Energy_initial = " << energy_initial << endl;


    while(time<=t_final){
    Derivative(r,v,drdt,dvdt,G,mass);

    double k1r[n];
    double k1v[n];
    double r_dummy[3];
    double v_dummy[3];
    for (int i = 0; i<n; i++){
        k1r[i] = dt*drdt[i];
        k1v[i] = dt*dvdt[i];
        r_dummy[i] = r[i] +k1r[i]/2.0;
        v_dummy[i] = v[i] +k1v[i]/2.0;
    }
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    double k2r[n];
    double k2v[n];
    for (int i = 0; i<n; i++){
        k2r[i] = dt*drdt[i];
        k2v[i] = dt*dvdt[i];
        r_dummy[i] = r[i] +k2r[i]/2.0;
        v_dummy[i] = v[i] +k2v[i]/2.0;
    }
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    double k3r[n];
    double k3v[n];
    for (int i = 0; i<n; i++){
        k3r[i] = dt*drdt[i];
        k3v[i] = dt*dvdt[i];
        r_dummy[i] = r[i] +k3r[i];
        v_dummy[i] = v[i] +k3v[i];
    }
    Derivative(r_dummy,v_dummy,drdt,dvdt,G,mass);
    double k4r[n];
    double k4v[n];
    for (int i = 0; i<n; i++){
        k4r[i] = dt*drdt[i];
        k4v[i] = dt*dvdt[i];

    }

    for (int i=0; i<n;i++){
        r[i] = r[i] +(1.0/6.0)*(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i]);
        v[i] = v[i] +(1.0/6.0)*(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);
    }
    time += dt;
    }
    cout << "final position:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r[i] << endl;
    }
    cout << "final velocity:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v[i] << endl;
    }

    double energy_final = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_final = " << energy_final<< endl;
    return 0;
}



