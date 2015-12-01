#include "star_cluster.h"


star_cluster::star_cluster()
{
}

 void star_cluster::add(stars n){
    number_stars++;
    all_stars.push_back(n);
}
void star_cluster::print_position(ofsteam &output, ofstream &output2, vector<stars> vec){
    print_position(output, output2, vec, 3)
}
void star_cluster::print_position(ofsteam &output, ofstream &output2, vector<stars> vec, int n){
    if (n>3 || n<=0) n = 3;
    for (int i = 0; i<vec.size(); i++){
        stars &this = vec[i];
        std::cout << std::scientific;
        for(int j = 0; j<n; j++){
            std::cout << this.position[j] << "   ";
            output << std::scientific << this.position[j] << "   ";
            output2 << std::scientific << this.velocity[j] << "   ";
        }
        std::cout << "     ";
        output << "       ";
        output2 << "        ";
        std::cout << std::endl;
        output << endl;
        output2 << endl;



    }
}


void star_cluster::insert_data(vector<stars> vec, mat &ma){
    for(int i=0; i<vec.size(); i++){
            stars &this = vec[i];
            ma(i,6)=this.mass;

            for(int k=0; k<3;k++){
                ma(i,k)=this.position[k];
                ma(i,k+3)=this.velocity[k];
            }
        }
}

