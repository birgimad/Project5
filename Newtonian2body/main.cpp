#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;



//void Rungekutta4_3dim(f, g, h, double t0, double dt, double r, double thet, double phi);
//void Rungekutta4(double t0, double dt, double &r, double &drdt, double &v, double &dvdt);
void Derivative(double &r, double &v, double &drdt, double &dvdt, double G, double mass);
int main()
{
    double G = 1;
    double mass = 1.0;
    double t0 = 0.0;
    double t_final = 1.0;
    int number_of_time_step = 10;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    int n = 3;
    double r[n];
    r[0] = 2.0;
    r[1] = 3.0;
    r[2] = 0.0;
    /*for (int i = 0; i<n;i++){
        cout<<r[i]<<endl;
    }*/
    double v[n];
    v[0] = 0.0;
    v[1] = 1.0;
    v[2] = 0.0;

    double drdt[n];

    double dvdt[n];
    while(time<=t_final){
    Derivative(r,v,drdt,dvdt, G, mass);
    double k1r[n];
    double k1v[n];
    double r_dummy[n];
    double v_dummy[n];
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
        r[i] = r[i] +(1.0/6.0)*(k1r[i]+k2r[i]+k3r[i]+k4r[i]);
        v[i] = v[i] +(1.0/6.0)*(k1v[i]+k2v[i]+k3v[i]+k4v[i]);
    }
    time += dt;
    }



   
}

void Derivative(double &r, double &v, double &drdt, double &dvdt, double G, double mass){
    drdt[0] = v[0];
    drdt[1] = v[1];
    drdt[2] = v[2];

    double distance_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double newtonian_force = G*mass/pow(distance_squared,1.5);
    dvdt[0] = newtonian_force*r[0];
    dvdt[1] = newtonian_force*r[1];
    dvdt[2] = newtonian_force*r[2];
}
/*
void Rungekutta4( double t0, double dt, double &r,double &drdt, double &v, double &dvdt){
    double n = r.size();
    for (int i = 0; i<n ;i++){
        double k1 = dt*drdt[i];
        double k2 = dt*f(t[i]+(dt/2.0), r[i] + (k1/2.0));
        double k3 = dt*f(t[i]+(dt/2.0), r[i] + (k2/2.0));
        double k4 = dt*f(t[i]+(dt), r[i] + (k3/2.0));
        r[i] = r[i] +(1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
}
}
*/
/*
void Rungekutta4_3dim( double t0, double dt, double &r, double &thet, double &phi)
{
    double t = t0;
    const int n = r.size();


    for (int i = 0; i<n; i++){
        double k1 = dt*f(t[i],r[i], thet[i], phi[i]);
        double kvx = dt*fv(vx[i], vy[i], vz[i]);
        double l1 = dt*g(t[i],r[i], thet[i], phi[i]);
        double o1 = dt*h(t[i],r[i], thet[i], phi[i]);
        double k2 = dt*f(t[i]+(dt/2.0), r[i] + (k1/2.0), thet[i] + (l1/2.0), phi[i] + (o1/2.0));
        double l2 = dt*g(t[i]+(dt/2.0), r[i] + (k1/2.0), thet[i] + (l1/2.0), phi[i] + (o1/2.0));
        double o2 = dt*h(t[i]+(dt/2.0), r[i] + (k1/2.0), thet[i] + (l1/2.0), phi[i] + (o1/2.0));
        double k3 = dt*f(t[i]+(dt/2.0), r[i] + (k2/2.0), thet[i] + (l2/2.0), phi[i] + (o2/2.0));
        double l3 = dt*g(t[i]+(dt/2.0), r[i] + (k2/2.0), thet[i] + (l2/2.0), phi[i] + (o2/2.0));
        double o2 = dt*h(t[i]+(dt/2.0), r[i] + (k2/2.0), thet[i] + (l2/2.0), phi[i] + (o2/2.0));
        double k4 = dt*f(t[i]+(dt), r[i] + (k3/2.0), thet[i] + (l3/2.0), phi[i] + (o3/2.0));
        double l4 = dt*g(t[i]+(dt), r[i] + (k3/2.0), thet[i] + (l3/2.0), phi[i] + (o3/2.0));
        double o4 = dt*h(t[i]+(dt), r[i] + (k3/2.0), thet[i] + (l3/2.0), phi[i] + (o3/2.0));
        r[i+1] = r[i] +(1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        thet[i+1] = thet[i] +(1.0/6.0)*(l1 + 2*l2 + 2*l3 + l4);
        phi[i+1] = phi[i] +(1.0/6.0)*(o1 + 2*o2 + 2*o3 + o4);
        t += dt;
    }



}
*/

