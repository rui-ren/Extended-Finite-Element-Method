// Use the geometryentity to find enriched nodes!
// Calculate enriched nodes!

#ifndef	_GEOMETRYENTITY_H_
#define	_GEOMETRYENTITY_H_

#include"FEM_Component.h"
#include"geometrydescription.h"
#include"geometry_info.h"
#include<stdio.h>
#include<iostream>
#include"geometry_2D.h"
using namespace std;

class Element;   // super class define the geometry of an enrichment item.

/*
Vertex, PiecewiseLinear,PicewiseParabolic, Circle and Ellipse.

Usage:
Enrichment items will use this geometry to find	enriched nodes.
*/

class	GeometryEntity :public FEMComponent
{

public:
	// Constructor
	GeometryEntity(int n, Domain* aDomain) :FEMComponent(n, aDomain){
		geoDescription = 0;
	}

	virtual ~GeometryEntity(){}// Destructor
	/*	Definition of an element
	*  Functions to	define an element
	*/
	GeometryEntity*				  typed();
	GeometryEntity*				  ofType(char*);
	char*	giveClassName(char* s)
	{
		return	strcpy(s, "GeometryEntity");
	}
	int	    giveNumber()
	{
		return	FEMComponent::giveNumber();
	}

	/*
	Do	the interaction between	the geometry entity and	the finite elements
	Use this function to check the geometry entity intersect with finite elements
	*/
	virtual bool interactsWith(Element*);

	/*
	Compute	the signed distance of a given point to the receiver, Heaviside equation
	Virtual	function
	*/

	virtual double computeSignedDistanceOfPoint(rui::Point*){ return NULL; }
	/*
	Read from the input file the type of description of	the receiver
	Possible geometry descriptions are standard, level set and vector level set.

	@return :
	StandardDescription (geoDescription=1)
	LevelSetDescription (geoDescription=2)
	VectorLevelSetDescription (geoDescription=3)
	*/


	GeometryDescription*	giveMyGeoDescription();    //return the Heaviside Value of the function
	/*
	Returns	the position of the point compared to the	receiver
	1	if	point	is	inside (circle) above (line) and	0 otherwise
	Virtual	function, so derived	class	have to implemente this	method
	*/

	virtual int	 givePositionComparedTo(rui::Point*)	{ return	NULL; }

	/*
	Check the intersection	between the	geometry entity and	the triangle
	*/
	virtual bool intersects(const rui::Triangle* t){ return	NULL; }
	virtual bool intersects(const rui::Segment* s){ return	NULL; }
	/*!
	Computes the intersection points between	the receiver with	the triangle
	*/
	virtual vector<rui::Point> intersection(const rui::Triangle* t)
	{
		return vector<rui::Point>();
	}

	virtual vector<rui::Point> intersection(const rui::Segment* s)
	{
		return vector<rui::Point>();
	}

	/*
	Computes the intersection points between the receiver with a given element
	Used to partition split elements into sub-cells
	*/
	virtual vector<rui::Point> intersection(Element* e);
	/*!
	Export geometry to Matlab file for plotting.
	*/
	virtual void  exportToMatlab(std::string&){ ; }   //  export the report to Matlab

protected:

	size_t	   geoDescription;// Description of Geometry entity
};


#endif // _GEOMETRYENTITY_H_
