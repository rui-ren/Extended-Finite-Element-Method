
#ifndef _FEMCOMPONENT_H_
#define _FEMCOMPONENT_H_
#include "string.h"
class Domain; 
class FileReader;
class FEMComponent
	/*
	This class is an abstract class, the superclass of all classes that imple-
	ment the components of a finite element mesh 
	contain: elements, nodes, time steps, materials, loads and load-time functions.

	This class defines the two attributes common to all component classes ;
	'number' is primarily used for reading data in the data file. 
	'domain' is used for communicating with other components 
	*/
{
protected:
	int      number;   // read data in the data file
	Domain*  domain;   // used for communicating with other components

public:
	FEMComponent() {}										// constructors
	FEMComponent(int n, Domain* d) { number = n; domain = d; }	 //default constructor
	virtual ~FEMComponent() {}								// destructor

	virtual char*  giveClassName(char*) = 0;
	virtual void   giveKeyword(char*);
	Domain*        giveDomain()const        { return domain; }
	int            giveNumber()const        { return number; }
	double         read(char);
	double         read(char* d)            { return this->read(d, 1); }
	double         read(char*, int);
	int            readIfHas(char*);
	int            readInteger(char* d)     { return this->readInteger(d, 1); }
	int            readInteger(char* d, int i);
	int            readNumberOf(char*);
	void           readString(char* d, char* s)  { this->readString(d, 1, s); }
	void           readString(char*, int, char*);
	int            readWhetherHas(char*);
};
#endif // _FEMCMPNN_H_
