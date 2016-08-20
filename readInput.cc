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
void readtexture(std::string filename){

	int i,j,k,n,m,p;
	int gid,pid;
	double a[3][3];
	double phi1,Phi,phi2,dummy;
	double cvoxel3333[3][3][3][3];
	double cvoxel66[6][6];
	double aux6[6],aux33[3][3];
	double saux[6][6],taux[6][6];

	std::string line;

	std::fstream textureIn;
	textureIn.open(filename.c_str(), std::ios::in);

	while(!textureIn.eof()){
		textureIn>>phi1>>Phi>>phi2>>i>>j>>k>>gid>>pid;
		euler[k-1][j-1][i-1][0]=phi1;
		euler[k-1][j-1][i-1][1]=Phi;
		euler[k-1][j-1][i-1][2]=phi2;
		grainID[k-1][j-1][i-1]=gid;
		phaseID[k-1][j-1][i-1]=pid;

		transformationMatrix(a,euler[k-1][j-1][i-1],2);
		transformFourthOrderTensor(cmat3333,cvoxel3333,a,1);
	
		change_basis(aux6,aux33,cvoxel66,cvoxel3333,4);

		for(n=0;n<6;n++)
			for(m=0;m<6;m++){
				cloc[k-1][j-1][i-1][n][m]=cvoxel66[n][m];
				xlsec66[n][m]+=cvoxel66[n][m]*wgt;
			}
	}

	change_basis(aux6,aux33,xlsec66,xlsec3333,3);

	for(k=0;k<N3;k++)
		for(j=0;j<N2;j++)
			for(i=0;i<N1;i++){
				
				for(n=0;n<6;n++)
					for(m=0;m<6;m++)
						saux[n][m]=cloc[k][j][i][n][m];
				findInverse((double *)saux,6);

				for(n=0;n<6;n++)
					for(m=0;m<6;m++){
						dummy=0.0;
						for(p=0;p<6;p++)
							dummy+=xlsec66[n][p]*saux[p][m];
						taux[n][m]=(i+1/(j+1))*(j+1/(i+1))+dummy;
					}

				findInverse((double *)taux,6);

				for(n=0;n<6;n++)
					for(m=0;m<6;m++){
						dummy=0.0;
						for(p=0;p<6;p++)
							dummy+=saux[n][p]*taux[p][m];
						fsloc[k][j][i][n][m]=dummy;
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

		voigt(cmat66,cmat3333,1);
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
						cmat3333[i][j][k][l]=lambda*identityR2[i][j]*identityR2[k][l]+2.0*mu*identityR4[i][j][k][l];

	}
	
}

void readinput(char filename[100]){
	
	std::string textureFile,propsFile;
	std::string line;
	double aux66[6][6],aux3333[3][3][3][3];
	int n1,n2,n3;
	int i,j;


	std::fstream maininputIn;
	maininputIn.open(filename, std::ios::in);

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

	//
	{
		std::getline(maininputIn,line);
		std::istringstream iss(line);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
	                     		   std::istream_iterator<std::string>{}};
		
		n1=std::stoi(tokens[0]);
		n2=std::stoi(tokens[1]);
		n3=std::stoi(tokens[2]);
	}
	
	//	
	{	
		std::getline(maininputIn,line);
		std::istringstream iss(line);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
	                     		   std::istream_iterator<std::string>{}};
		
		RVEdim[0]=std::stod(tokens[0]);
		RVEdim[1]=std::stod(tokens[1]);
		RVEdim[2]=std::stod(tokens[2]);
	}
	
	//
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

		for(i=0;i<3;i++)
			for (j=0;j<3;j++){
				straingradrate33[i][j]=0.5*(velgrad33[i][j]+velgrad33[j][i]);
				rotationrate33[i][j]=0.5*(velgrad33[i][j]+velgrad33[j][i]);
			}

		for(i=0;i<3;i++)
			IDstraingradrate[i]=straingradrate33[i][i];

		IDstraingradrate[3]=0;
		if(velgrad33[1][2]==1 && velgrad33[3][2]==1)IDstraingradrate[3]=1;
		IDstraingradrate[4]=0;
		if(velgrad33[0][2]==1 && velgrad33[2][0]==1)IDstraingradrate[4]=1;
		IDstraingradrate[5]=0;
		if(velgrad33[0][1]==1 && velgrad33[1][0]==1)IDstraingradrate[5]=1;

		change_basis(straingradrate6,straingradrate33,aux66,aux3333,2);
	}

	//
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

	//
	{
		std::getline(maininputIn,line);
		std::istringstream iss(line);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
	                     		   std::istream_iterator<std::string>{}};
		
		error=std::stod(tokens[0]);
	}

	//
	{
		std::getline(maininputIn,line);
		std::istringstream iss(line);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
	                     		   std::istream_iterator<std::string>{}};
		
		nsteps=std::stoi(tokens[0]);
	}

	//
	{
		std::getline(maininputIn,line);
		std::istringstream iss(line);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
	                     		   std::istream_iterator<std::string>{}};
		
		itermax=std::stoi(tokens[0]);
	}
}
