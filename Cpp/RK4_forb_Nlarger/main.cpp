#include <iostream>
#include "include/armadillo"
using namespace std;
using namespace arma;

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
    int n=3;

    double G = 2.96e-4;     //Grav const in the units of AU^3 / ( days^3 * mass_sun )
    double t0 = 0.0;
    double t_final = 20*365.0;     //time in the unit of days
    int number_of_time_step = 100000;
    double time = t0;
    double dt = (t_final - t0)/number_of_time_step;
    int number_of_particles = 3;
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
    r_initial = r;
    v_initial = v;

    mat pos_temp(number_of_particles,3);
    mat vel_temp(number_of_particles,3);
    mat k1(number_of_particles,7);
    mat k2(number_of_particles,7);
    mat k3(number_of_particles,7);
    mat k4(number_of_particles,7);
    pos_temp.zeros();
    vel_temp.zeros();
    k1.zeros();
    k2.zeros();
    k3.zeros();
    k4.zeros();

    cout << "initial position sun:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(0,i) << endl;
    }
    cout << "initial velocity sun:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(0,i) << endl;
    }

    cout << "initial position earth:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(1,i) << endl;
    }
    cout << "initial velocity earth:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(1,i) << endl;
    }
    cout << "initial position mars:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(2,i) << endl;
    }
    cout << "initial velocity mars:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(2,i) << endl;
    }
    mat distance_between_particles(number_of_particles,number_of_particles);
    distance_between_particles.zeros();

    for (int i = 0; i<number_of_particles; i++)
    {
        for (int j = 0; j<number_of_particles; j++)
        {
            for (int k = 0; k<3; k++)
            {
                distance_between_particles(i,j) += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            distance_between_particles(i,j) = pow(distance_between_particles(i,j),0.5);
        }

    }
    cout << "distance between earth and sun: " << distance_between_particles(0,1) << " AU" << endl;
    cout << "distance between mars and sun: " << distance_between_particles(0,2) << " AU" << endl;

    vec kin_en(number_of_particles);
    kin_en.zeros();
    vec pot_en(number_of_particles);
    pot_en.zeros();
    vec tot_en(number_of_particles);
    tot_en.zeros();
    for (int i=0; i<number_of_particles; i++)
    {
        for (int k=0; k<3; k++)
        {
            kin_en(i) += v(i,k)*v(i,k);
        }
        kin_en(i) = 0.5*m(i)*kin_en(i);
        for (int j=0; j<number_of_particles; j++)
        {
            if (j != i)
            {
                pot_en(i) += pow(distance_between_particles(i,j),-1.0)*m(j);
            }
        }
        pot_en(i) = pot_en(i)*G*m(i);
        tot_en(i) = kin_en(i)+pot_en(i);
    }
    cout << "Energy Sun: " << tot_en(0) << endl;
    cout << "Energy Earth: " << tot_en(1) << endl;

    mat drdt(number_of_particles,3);
    drdt.zeros();
    mat dvdt(number_of_particles,3);
    dvdt.zeros();

    while(time<=t_final){
    Derivative(r,v,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k1,2.0,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k2,2.0,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    making_ks(pos_temp,vel_temp,r,v,drdt,dvdt,dt,k3,1,number_of_particles);
    Derivative(pos_temp,vel_temp,m,drdt,dvdt,G,number_of_particles);
    for (int j=0; j<number_of_particles; j++)
    {
        for (int i = 0; i<3; i++){
            k4(j,i) = dt*drdt(j,i);
            k4(j,i+3) = dt*dvdt(j,i);
        }
        for (int i=0; i<3; i++)
        {
            r(j,i) += (1.0/6.0)*(k1(j,i)+2*k2(j,i)+2*k3(j,i)+k4(j,i));
            v(j,i) += (1.0/6.0)*(k1(j,i+3)+2*k2(j,i+3)+2*k3(j,i+3)+k4(j,i+3));
        }
    }
    time += dt;
    }

    cout << "final position sun:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(0,i) << endl;
    }
    cout << "final velocity sun:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(0,i) << endl;
    }

    cout << "final position earth:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(1,i) << endl;
    }
    cout << "final velocity earth:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(1,i) << endl;
    }
    cout << "final position mars:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << r(2,i) << endl;
    }
    cout << "final velocity mars:" << endl;
    for (int i = 0; i<n; i++)
    {
        cout << v(2,i) << endl;
    }
    mat distance_between_particles_final(number_of_particles,number_of_particles);
    distance_between_particles_final.zeros();

    for (int i = 0; i<number_of_particles; i++)
    {
        for (int j = 0; j<number_of_particles; j++)
        {
            for (int k = 0; k<3; k++)
            {
                distance_between_particles_final(i,j) += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            distance_between_particles_final(i,j) = pow(distance_between_particles_final(i,j),0.5);
        }

    }
    cout << "distance between earth and sun: " <<  distance_between_particles_final(1,0) << " AU" << endl;
    cout << "distance between mars and sun: " <<  distance_between_particles_final(2,0) << " AU" << endl;

    vec kin_en_final(number_of_particles);
    kin_en_final.zeros();
    vec pot_en_final(number_of_particles);
    pot_en_final.zeros();
    vec tot_en_final(number_of_particles);
    tot_en_final.zeros();
    for (int i=0; i<number_of_particles; i++)
    {
        for (int k=0; k<3; k++)
        {
            kin_en_final(i) += v(i,k)*v(i,k);
        }
        kin_en_final(i) = 0.5*m(i)*kin_en(i);
        for (int j=0; j<number_of_particles; j++)
        {
            if (j != i)
            {
                pot_en_final(i) += pow(distance_between_particles_final(i,j),-1.0)*m(j);
            }
        }
        pot_en_final(i) = pot_en_final(i)*G*m(i);
        tot_en_final(i) = kin_en_final(i)+pot_en_final(i);
    }
    cout << "Energy Sun: " << tot_en_final(0) << endl;
    cout << "Energy Earth: " << tot_en_final(1) << endl;
/*
    ofstream myfile ("RungeKutta4_2body3D_newscript.txt");
        if (myfile.is_open())
        {
            myfile << "Runge-Kutta Method, 2 body, 3D, new script" << endl;
            myfile << "Time: " << t_final << " days" << endl;
            myfile << "Number of time steps: " << number_of_time_step << endl;
            myfile << "Time step: " << dt << " days" << endl;
            myfile << " " << endl;
            myfile << "mass Sun" << m(0) << endl;
            myfile << "mass Earth" << m(1) << endl;
            myfile << " " << endl;
            myfile << "initial position Sun:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r_initial(0,i) << endl;
            }
            myfile << "initial velocity Sun:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v_initial(0,i) << endl;
            }
            myfile << "initial position Earth:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r_initial(1,i) << endl;
            }
            myfile << "initial velocity Earth:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v_initial(1,i) << endl;
            }
            myfile << " " << endl;
            myfile << "Initial distance earth-sun: " << distance_between_particles(1,0) << " AU" << endl;
            myfile << "Initial energy Sun:  " <<  tot_en(0) << endl;
            myfile << "Initial energy Earth:  " <<  tot_en(1) << endl;
            myfile << " " << endl;
            myfile << "final position Sun:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r(0,i) << endl;
            }
            myfile << "final velocity Sun:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v(0,i) << endl;
            }
            myfile << "final position Earth:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << r(1,i) << endl;
            }
            myfile << "final velocity Earth:" << endl;
            for (int i = 0; i<n; i++)
            {
                myfile << v(1,i) << endl;
            }
            myfile << " " << endl;
            myfile << "Final distance earth-sun: " << distance_between_particles_final(1,0) << " AU" << endl;
            myfile << "Final energy Sun:  " <<  tot_en_final(0) << endl;
            myfile << "Final energy Earth:  " <<  tot_en_final(1) << endl;
            myfile << " " << endl;
        }
*/

    return 0;
}

