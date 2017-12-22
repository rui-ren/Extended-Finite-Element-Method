
#ifndef _PIECEWISELINEAR_H_
#define _PIECEWISELINEAR_H_

#include "geometryentity.h"
#include "domain.h"
#include "Geometry_info.h"
#include <list>
#include <vector>
#include <iostream>
#include "geometry_2D.h"
using namespace std;

class Vertex;
class GeometryDescription; 
class CrackInterior;

//! define a straight crack as a polyline.
/*
This class implements a polyline composed by several vertices
A PiecewiseLinear is geometry of a straight crack
*/

class PiecewiseLinear :public GeometryEntity
{
public:

	PiecewiseLinear(int, Domain*);   // Constructor
	~PiecewiseLinear();				// Destructor

	/*List of vertices
	*  Methods to get and give vertices of the polyline
	*/
	Vertex*				 giveVertex(int);
	list<rui::Point*>*   getListOfVertices();     // find vertice points
	list<rui::Point*>*   giveMyListOfVertices();  // store vertice points

	/*
	Find the segment of polyline closest to a given point
	returns two vertices of this segment
	*/
	rui::Segment FindSegmentClosestTo(rui::Point*);   // find the element need to be enriched

	/*
	Compute the signed distance of a given point to
	the receiver. signed distance function.
	*/
	double computeSignedDistanceOfPoint(rui::Point*);

	/*
	Check the position of a point w.r.t the receiver
	@return  1 if p is above the receiver
	@return -1 if p is below the receiver
	@return  0 if p is on the receiver
	*/

	int givePositionComparedTo(rui::Point* p);
	/*!
	Make a vector of segments
	*/

	vector<rui::Segment*>*   makeSegments();
	/*!
	Return vector of segments of the receiver
	*/

	vector<rui::Segment*>*   giveSegments();

	/*!
	Check the intersection of the receiver with a triangle
	*/
	bool   intersects(const rui::Triangle* tri);
	bool   intersects(const rui::Segment* s);

	// returns the extreme segments whose end is the vertex v (often v is the TIP)
	// v is the vertex
	// 
	rui::Segment*   giveSegmentContain(Vertex* v);

	/*
	Compute the intersection of the receiver with a triangle
	tri the triangle
	*/
	vector<rui::Point> intersection(const rui::Triangle* tri);
	vector<rui::Point> intersection(const rui::Segment* s);

	// Defines the bounding box of the reveiver
	/*
	Export geometry to Matlab file for plotting.
	*/
	void  exportToMatlab(std::string&);

	void  insertNewVertexAtBack(Vertex* v);
	void  insertNewVertexInFront(Vertex* v);
	void  insertNewSegmentInFront(rui::Segment *s){ segmentList->insert(segmentList->begin(), s); }
	void  insertNewSegmentAtBack(rui::Segment *s){ segmentList->push_back(s); }

	/*
	Check if a point belongs to the receiver
	*/

	bool  isBelongToMe(rui::Point* p);

private:
	std::list<rui::Point*>*     vertexList;   // list of vertices of the PiecewiseLinear
	std::vector<rui::Segment*>* segmentList;   // segments of the Piecewiselinear
	//rui::Rectangle*             myBoundingBox; // bounding box of the PL
};

#endif // _PIECEWISELINEAR_H_