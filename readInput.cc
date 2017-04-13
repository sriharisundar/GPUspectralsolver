#include "readInput.h"
#include "matrixOperations.h"
#include "globalVariables.h"
#include "printFunctions.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <iterator>
#include <stdlib.h> //for atoi()
#include <iomanip>

void readtexture(std::string filename){

    int i,j,k,n,m,p,count;
    int gid,pid;
    double a[3][3];
    double phi1,Phi,phi2,dummy,det=0;
    fourthOrderTensor cvoxel3333;
    double cvoxel66[6][6];
    double aux6[6],aux33[3][3];
    double saux[6][6],taux[6][6];
    double prodDim=n1*n2*n3;
    double volumeVoxel=1.0/prodDim;

    std::string line;

    std::fstream textureIn;
    textureIn.open(filename.c_str(), std::ios::in);

    for(count=0;count<prodDim;count++){
        textureIn>>phi1>>Phi>>phi2>>i>>j>>k>>gid>>pid;
        euler[(k-1)*n2*(n1)+(j-1)*(n1)+i-1+0]=phi1;
        euler[(k-1)*n2*(n1)+(j-1)*(n1)+i-1+1]=Phi;
        euler[(k-1)*n2*(n1)+(j-1)*(n1)+i-1+2]=phi2;
        grainID[(k-1)*n2*(n1)+(j-1)*(n1)+i-1]=gid;
        phaseID[(k-1)*n2*(n1)+(j-1)*(n1)+i-1]=pid;

        if(pid==2){
            for(m=0;m<6;m++)
                for(n=0;n<6;n++)
                    cloc[((k-1)*n2*n1+(j-1)*n1+i-1)*36+m*6+n]=0;                    
        }
        else{
            transformationMatrix(a,&euler[(k-1)*n2*(n1)+(j-1)*(n1)+i-1],2);
            transformFourthOrderTensor(cmat3333.tensor,cvoxel3333.tensor,a,1);
            
            change_basis(aux6,aux33,cvoxel66,cvoxel3333.tensor,4);

            for(m=0;m<6;m++)
                for(n=0;n<6;n++){
                    cloc[((k-1)*n2*n1+(j-1)*n1+i-1)*36+m*6+n]=cvoxel66[m][n];
                    C0_66[m][n]=C0_66[m][n]+cvoxel66[m][n]*volumeVoxel;
            }
        }
    }

    change_basis(aux6,aux33,C0_66,C0_3333.tensor,3);

    for(k=0;k<n3;k++)
        for(j=0;j<n2;j++)
            for(i=0;i<n1;i++){
                
                for(n=0;n<6;n++)
                    for(m=0;m<6;m++)
                        saux[n][m]=cloc[(k*n2*(n1)+j*(n1)+i)*36+6*n+m];
        
                findInverse((double *)saux,det,6);

                for(n=0;n<6;n++)
                    for(m=0;m<6;m++){
                        dummy=0.0;
                        for(p=0;p<6;p++)
                            dummy+=C0_66[n][p]*saux[p][m];
                        taux[n][m]=((n+1)/(m+1))*((m+1)/(n+1))+dummy;
                    }

                findInverse((double *)taux,det,6);

                for(n=0;n<6;n++)
                    for(m=0;m<6;m++){
                        dummy=0.0;
                        for(p=0;p<6;p++)
                            dummy+=saux[n][p]*taux[p][m];
                        fsloc[(k*n2*(n1)+j*(n1)+i)*36+6*n+m]=dummy;
                    }
    
    }
}

void readprops(std::string filename){
    char line[100];
    int i,j,k,l;
    double youngs,poisson,mu,lambda;
    double cmat66[6][6];

    std::fstream propsIn;

    propsIn.open(filename.c_str(), std::ios::in);

    propsIn.getline(line,99);
    if(atoi(line)==0){
        std::cout<<"Not isotropic"<<std::endl;

        for(i=0;i<6;i++)
            for(j=0;j<6;j++)
                propsIn>>cmat66[i][j];

        voigt(cmat66,cmat3333.tensor,1);
    }

    else{
        std::cout<<"Isotropic"<<std::endl;
        propsIn>>youngs>>poisson;

        mu=youngs/(2.0*(1+poisson));
        lambda=2.0*mu*poisson/(1.0-2.0*poisson);

        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                for(k=0;k<3;k++)
                    for(l=0;l<3;l++)
                        cmat3333.tensor[i][j][k][l]=lambda*identityR2[i][j]*identityR2[k][l]+2.0*mu*identityR4[i][j][k][l];

    }
    
}

void readinput(char filename[100]){
    
    std::string textureFile,propsFile;
    std::string line;
    double aux66[6][6],aux3333[3][3][3][3];
    int i;


    std::fstream maininputIn;
    maininputIn.open(filename, std::ios::in);

    // Read dimensions
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        n1=std::stoi(tokens[0]);
        n2=std::stoi(tokens[1]);
        n3=std::stoi(tokens[2]);
    }

    initglobal();

    //Read in propsfile name and execute reading that file
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
        
        propsFile=tokens[0];
        readprops(propsFile);
    }
    
    
    //Read in texturefile name and execute reading that file
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
        
        textureFile=tokens[0];
        readtexture(textureFile);

    }

    //Read in ouputfile name and store it in outfile variable
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
        
        std::strcpy(outputFile, tokens[0].c_str());
    }

    // Read in RVE dimensions
    {    
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        RVEdim[0]=std::stod(tokens[0]);
        RVEdim[1]=std::stod(tokens[1]);
        RVEdim[2]=std::stod(tokens[2]);
    }
    
    // Read in boundary conditions
    {
        for(i=0;i<3;i++){
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        velgrad33[i][0]=std::stod(tokens[0]);
        velgrad33[i][1]=std::stod(tokens[1]);
        velgrad33[i][2]=std::stod(tokens[2]);
        }

        symmetric((double *) velgrad33,(double *) straingradrate33);
        symmetric((double *) velgrad33,(double *) rotationrate33);

        for(i=0;i<3;i++)
            IDstraingradrate[i]=straingradrate33[i][i];

        IDstraingradrate[3]=0;
        if(velgrad33[1][2]==1 && velgrad33[3][2]==1)IDstraingradrate[3]=1;
        IDstraingradrate[4]=0;
        if(velgrad33[0][2]==1 && velgrad33[2][0]==1)IDstraingradrate[4]=1;
        IDstraingradrate[5]=0;
        if(velgrad33[0][1]==1 && velgrad33[1][0]==1)IDstraingradrate[5]=1;

        change_basis(strainbar,straingradrate33,aux66,aux3333,2);
    }

    // Read in control parameters, no idea what this does
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        ictrl=std::stoi(tokens[0]);
    
        if(ictrl<=3){
            ictrl1=ictrl-1;
            ictrl2=ictrl-1;
        }
        else if(ictrl==4){
            ictrl1=1;
            ictrl2=2;
        }
        else if(ictrl==5){
            ictrl1=0;
            ictrl2=2;
        }
        else{
            ictrl1=1;
            ictrl2=2;
        }

        strainref=straingradrate33[ictrl1][ictrl2];


    }

    // Error criteria for convergence
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        error=std::stod(tokens[0]);
    }

    // Total number of steps
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        nsteps=std::stoi(tokens[0]);
    }

    // Maximum iterations per step
    {
        std::getline(maininputIn,line);
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
        
        itermax=std::stoi(tokens[0]);
    }
}
