#include "geometryentity.h"
#include "circle.h"  
#include "piecewiselinear.h"
#include "standarddescription.h" 
#include "vectorlevelsetdescription.h"
#include "levelsetdescription.h"
#include "domain.h"
#include "node.h"
#include "assert.h"
#include <stdlib.h>
#include <stdio.h>

GeometryEntity*  GeometryEntity::ofType(char* aClass)
// ***********************************************************************
// Returns a new geometry item, which has the same number than the receiver,
// but belongs to aClass (PiecewiseLinear, Vertex or Circle).
{
	GeometryEntity *newGeoEntity;

	if (!strncmp(aClass, "PiecewiseLinear", 5))
		newGeoEntity = new PiecewiseLinear(number, domain);
	else if (!strcmp(aClass, "Circle"))
		newGeoEntity = new Cercle(number, domain);
	else
		assert(false);

	return newGeoEntity;
}

GeometryEntity*  GeometryEntity::typed()
// ************************************************************************
// Returns a new geometry item, which has the same number than the receiver,
// but belongs to aClass (Piecewiselinear, Vertex or Circle).
{
	GeometryEntity* newGeoEntity;
	char     type[32];

	this->readString("class", type);
	newGeoEntity = this->ofType(type);

	return newGeoEntity;
}

GeometryDescription* GeometryEntity::giveMyGeoDescription()
// ********************************************************
{
	if (geoDescription == 0)
		geoDescription = this->readInteger("geoDescription");

	GeometryDescription * geoDes;
	switch (geoDescription){
	case 1:
		geoDes = new StandardDescription;
		break;
	case 2:
		geoDes = new LevelSetDescription;
		break;
	}

	return geoDes;
}

bool GeometryEntity::interactsWith(Element* e)
// *********************************************
// Depending on the geometry description 
{
	return this->giveMyGeoDescription()->interactsWith(e, this);
}

std::vector<Mu::Point> GeometryEntity::intersection(Element* e)
// **************************************************************
// From the nodes of element, make segments, then compute the intersection
// between segments and the receiver.
// Remark: \todo How about curved elements,ie. high order elements?
{
	std::vector<Mu::Point>  temp, pts, ret;
	for (size_t i = 0; i < e->giveNumberOfNodes(); i++)
	{
		Mu::Point *p = e->giveNode(i + 1)->makePoint();
		pts.push_back(*p);
		delete p;
	}

	for (size_t i = 0; i < pts.size(); i++)
	{
		Mu::Segment s(pts[i], pts[(i + 1) % pts.size()]);
		if (this->intersects(&s))
		{
			temp = this->intersection(&s);
			ret.insert(ret.end(), temp.begin(), temp.end());
		}
	}

	return ret;
}

