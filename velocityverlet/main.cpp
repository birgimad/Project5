#include <iostream>
#include <armadillo>
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


int main()
{
    int n=3;

    double G = 2.96e-4;
    //double mass_of_earth = 3e-6;
    double t0 = 0.0;
    double mass =  1.0;
    double t_final = 182.0;
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
    double v_partly[3];
    v_partly[0] = 0;
    v_partly[1] = 0;
    v_partly[2] = 0;
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
     cout << "Energy_initial = " << energy_initial << endl;

    while(time<=t_final){
    Derivative(r,v,drdt,dvdt,G,mass);

    for(int i=0; i<n ; i++){
    r[i] = r[i]+dt*drdt[i]+0.5*dt*dt*dvdt[i];
    v_partly[i] = drdt[i]+0.5*dt*dvdt[i];
    dvdt[i] = v_partly[i];
    //va[i] = drdt[i];
    }
    Derivative(r,v,drdt,dvdt,G,mass);

    for(int i=0; i<n ; i++){
    v[i] = v_partly[i]+0.5*dt*dvdt[i];
    //v1[i] = drdt[i];
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


        cout << drdt[i] << endl;
    }

    double energy_final = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]*0.5+G*mass*pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],-0.5);
    cout << "Energy_final = " << energy_final<< endl;
    return 0;
}

