//
//  main.cpp
//  Ising
//
//  Created by Parisa on 12/12/2017.
//  Copyright © 2017 Parisa. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <random>
#include<fstream>


using namespace std;

const double J=1;
const int u=1;
const int d=-1;

class Ising{

public:
    int* list;
    int D,Steps,Dim;
    double M,E,T,Mu,H;
    double averageM(void);
    void totE(void);
    double averageE(void);
    double averageMsq(void);
    Ising(int size=5, double t=4,int steps=1000,int dim=2,double mu=0,double h=0);
    ~Ising()  {  if(list) delete [] list; }
    
private:
    
    void step (void);
    bool flip(int x,int y);
    int & spin(int x, int y);

};

Ising::Ising(int size, double t,int steps,int dim,double mu,double h){
    D=size;
    T=t;
    M=0;
    E=0;
    Mu=mu;
    H=h;
    Steps=steps;
    Dim=dim;
    int temp = pow(D,Dim);
    list= new int[temp];
    for(int i=0; i<temp; i++)
    {
        if (drand48()>0.5) //creates initial random configuration of spins
        {  list[i] = u;
            M++;
        }
        else
        {  list[i] = d;
            M--;
        }
    }
    totE();

}


void Ising::totE(void)
{  int n = 0;    // nearest neighbour
    
    for(int i=0; i<D; i++)
        for(int j=0; j<D; j++)
        {  if(spin(i,j)==spin(i+1, j))
            n++;
        else
            n--;
            if(spin(i, j)==spin(i, j+1))
                n++;
            else
                n--;
        }
    
    E = -J*n-Mu*H*M;
}


int & Ising::spin(int x, int y){
    if(x<0) x+=D; //periodic boundary conditions
    if(x>D) x-=D;
    if(y<0) y+=D;
    if(y>D) y-=D;
    int pos=x*D+y;
    return list[pos]; //returns the value of the list indexed position x,y
}



bool Ising::flip(int x,int y){
    const double kb=1;
    double Eflip=2*J*spin(x,y)*(spin(x,y+1)+spin(x+1,y)+spin(x,y-1)+spin(x-1,y))+2*Mu*H*spin(x,y);
    double deltaM= -2*spin(x,y);
    
    if (Eflip<9 && Eflip>-9){   //ensures spins do not go out of the range of values they can take
        if (Eflip<0){
            spin(x, y) *= -1; //filps spin if energy change is negative
            M += deltaM; //updates the magnetisation and energy of the system
            E += Eflip;
            return 1;
            }else{ if(drand48()<exp(-Eflip/(kb*T)))
                {spin(x, y) *= -1;// flips spin(x,y) if the boltzmann factor is greater than a random number between 0-1
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
        for(int j=0;j<D;j++){
            flip(i,j); //carry out one monte carlo time step so each spin in the lattice has 1 chance to flip
        }
}


double Ising::averageM(void){
    double totM=0;
    for(int i=0;i<Steps;i++){
        totM += M/(D*D);
        step();
    }
    return totM/Steps;
}

double Ising::averageMsq(void){
    double totMsq=0;
    for(int i=0;i<Steps;i++){
        totMsq += pow(M/(D*D),2);
        step();
    }
    return totMsq/Steps;
}


double Ising::averageE(void){
    double totE=E;
    for(int i=0;i<Steps;i++){
        totE += E/(D*D);
        step();
    }
    return totE/Steps;
}




int main(int argc, const char * argv[]) {
    //float sum = 0;
    ofstream myfile;
    int lattice=200;
    int iter=1000000;
    myfile.open("varMdata100.txt");
    myfile << "lattice" << lattice;
    myfile << "steps" << iter << '\n';
    for(float k=1.5;k<4;k+=0.05){
        Ising I(lattice,k,iter);
        myfile << k << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }
    for(float k=1.9;k<2.25;k+=0.01){
        Ising I(lattice,k,iter);
        myfile << k << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }
    for(float k=2.3;k<4;k+=0.1){
        Ising I(lattice,k,iter);
        myfile << k << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }
    myfile.close();
    
    int temp =0.1;
    myfile.open("hyster2.txt");
    myfile << "lattice" << lattice;
    myfile << "Temp" << temp;
    myfile << "steps" << iter << '\n';
    
    for(float k=-10;k<10;k+=0.05){
        Ising I(lattice,temp,iter,2,1,k);
        myfile << k << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }
    myfile.close();
    myfile.open("hyster_minus2.txt");
    myfile << "lattice" << lattice;
    myfile << "steps" << iter << '\n';
    for(float k=10;k>-10;k-=0.05){
        Ising I(lattice,temp,iter,2,1,k);
        myfile << k << "  ";
        myfile << I.averageM() << "  ";
        myfile << I.averageMsq() << '\n';
    }


    
    
    myfile.close();
    

    cout << "done"<<'\n';
    return 0;
}

