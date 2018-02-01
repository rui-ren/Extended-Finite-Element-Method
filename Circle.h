#ifndef _CIRCLE_H_
#define _CIRCLE_H_

#include "geometryentity.h"
#include "vertex.h"
#include "flotarry.h"
#include "domain.h"
#include "geometry_2D.h"
#include <iostream>

/*! A Circle is geometry of holes
This class implements a circle with center P and radius R
*/

class Cercle :public GeometryEntity
{
public:

	/*!
	Constructor
	*/
	Cercle(int n, Domain* d) :GeometryEntity(n, d), radius(0), center(0){};
	~Cercle(){ ; }           //!< Destructor

	/*!
	Return the center of the receiver
	*/
	Vertex*    giveCenter();
	/*!
	Return the radius of the receiver
	*/
	double     giveRadius();
	/*!
	Check if a point defined by coord is inside (negative) or
	outside the circle (positive value)
	*/
	double     computeSignedDistanceOfPoint(rui::Point*);
	/*!
	Returns the position of the point compared to the receiver
	@return 1 if point is inside and @return 0 otherwise
	*/
	int        givePositionComparedTo(rui::Point*);
	/*!
	Check the intersection between circle and triangle
	@return true if intersects
	*/
	bool       intersects(const rui::Triangle* t);
	bool       intersects(const rui::Segment* s);
	/*!
	Computes the intersection points of circle and triangle
	@param t the triangle
	@return intersection points stored in vector<Rui::Point>
	*/
	std::vector<rui::Point> intersection(const rui::Triangle* t);
	std::vector<rui::Point> intersection(const rui::Segment* s);
	/*!
	Export geometry to Matlab file for plotting.
	Not yet implemented at this point.\todo
	*/
	void  exportToMatlab(std::string&);


private:
	int      center;
	double   radius;
};

#endif // _CIRCLE_H_
