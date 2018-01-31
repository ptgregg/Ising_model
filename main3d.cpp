//
//  main.cpp
//  Ising3d
//
//  Created by Parisa on 03/01/2018.
//  Copyright © 2018 Parisa. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <random>
#include<fstream>


using namespace std;

const double J=1;

class Ising{
    
public:
    int* list;
    int D,Steps;
    double M,E,T;
    double averageM(void);
    double averageE(void);
    void totE(void);
    double averageMsq(void);
    Ising(int d=5, double t=4,int steps=100000);
    ~Ising()  {  if(list) delete [] list; }
    
private:
    void step (void);
    int & spin(int x, int y,int z);
    bool flip(int x,int y,int z);
};

Ising::Ising(int dim, double t,int steps){
    D=dim;
    T=t;
    M=0;
    E=0;
    Steps=steps;
    int temp = D*D*D;
    list= new int[temp];
    for(int i=0; i<temp; i++)
    {
        if (drand48()>0.5)
        {  list[i] = 1;
            M++;
        }
        else
        {  list[i] = -1;
            M--;
        }
    }
    totE();
}



void Ising::totE(void){
      int n = 0;    // nearest neighbour
        
        for(int i=0; i<D; i++)
            for(int j=0; j<D; j++)
                for(int k=0; k<D; k++)
            {  if(spin(i,j,k)==spin(i+1, j,k))
                n++;
            else
                n--;
                
                if(spin(i, j,k)==spin(i, j+1,k))
                    n++;
                else
                    n--;
                
                if(spin(i, j,k)==spin(i, j,k+1))
                    n++;
                else
                    n--;
            }
        
        E = -J*n;
}

int & Ising::spin(int x, int y, int z){
    if(x<0) x+=D; //periodic boundary conditions
    if(x>D) x-=D;
    if(y<0) y+=D;
    if(y>D) y-=D;
    if(z<0) z+=D;
    if(z>D) z-=D;
    int pos=(x*D+y)*D+z; //maps x,y,z grid position to a list index
    return list[pos];//returns the spin of x,y grid postion in list
    
}



bool Ising::flip(int x,int y,int z){
    const double kb=1;
    double Eflip=2*J*spin(x,y,z)*(spin(x+1,y,z)+spin(x-1,y,z)+spin(x,y+1,z)+spin(x,y-1,z)+spin(x,y,z+1)+spin(x,y,z-1));
    double deltaM= -2*spin(x,y,z);

    if (Eflip<13 && Eflip>-13){
        if (Eflip<0){
            spin(x, y,z) *= -1;
            M += deltaM;
            E += Eflip;
            return 1;
        }else{ if(drand48()<exp(-Eflip/(kb*T)))
            {spin(x, y,z) *= -1;
                M += deltaM;
                E += Eflip;
                return 1;
            }else{
                M+=0;
                E+=0;
                return 0;
            }
        }
        
    }else{
        return 0;
    }
    
    
}



void Ising::step(void){
    for(int i=0;i<D;i++)
        for(int j=0;j<D;j++)
            for(int k=0;k<D;k++){
                flip(i,j,k);
        }
}


double Ising::averageM(void){
    double totM=0;
    for(int i=0;i<Steps;i++){
        totM += M/(D*D*D);
        step();
    }
    return totM/Steps;
}

double Ising::averageMsq(void){
    double totMsq=0;
    for(int i=0;i<Steps;i++){
        totMsq += pow(M/(D*D*D),2);
        step();
    }
    return totMsq/Steps;
}

double Ising::averageE(void){
    double totE=0;
    for(int i=0;i<Steps;i++){
        totE += E/(D*D*D);
        step();
    }
    return totE/Steps;
}








int main(int argc, const char * argv[]) {
    ofstream myfile;
    myfile.open("Mcritdatazoom2.txt");
    int lattice=30;
    int iter=10000;
    myfile << "lattice" << lattice;
    myfile << "steps" << iter << '\n';
    for(float i=1.5;i<6;i+=0.1){
        Ising I(lattice,i,iter);
        myfile << i << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }
    
    myfile.close();
    
    
    cout << "done"<<'\n';
    return 0;
}

