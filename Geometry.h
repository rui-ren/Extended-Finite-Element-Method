#ifndef GEOMETRY_H
#define GEOMETRY_H
#include<math.h>
#include<vector>
#include<list>
#include<iostream>

using namespace std;

class Geometry
{
public:
	double Crack_Geometry(double& X1, double& X2, double& Y1, double& Y2);
	void Boundary(vector<double> X, vector<double> Y);

protected:
	// The Geometry Parameter of domain;
	vector <double> X;       // The length of the domain
	vector <double> Y;       // The height of the domain

	// The fundamental parameter of the formation;
	double E;		// Young's Modulus
	float Poss;		// Possion's Ratio
	float Pore;		// Pore Pressure
	float Biot;

	// The fundamental parameter of the drilling mud;
	double em;		    // Mud Density
	double Rheology;   // rheology of mud
	double Depth;

	// The fundamental parameter of the fractures; (Linear)
	double Length_Crack;
	double X1; 
	double Y1;
	double X2;
	double Y2;
};

// Crack geometry
double Crack_Geometry(double& X1, double& X2, double& Y1, double& Y2)
{
	double X_Crack = NULL;
	double Crack_theta = (Y2 - Y1) / (X2 - X1);
	double b = Y1 - Crack_theta*X1;
	double Crack = Crack_theta*X_Crack + b;
	return Crack;
};

#endif
