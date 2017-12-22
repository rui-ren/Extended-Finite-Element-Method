// This function use to define the elment which can be classified into:
// standard Description
// vector level set Description
// level set Description

#ifndef _GEOMETRYDESCRIPTION_H_
#define _GEOMETRYDESCRIPTION_H_

class Element; 
class GeometryEntity;

// The description type of geometry entities

class GeometryDescription
{
public:
	GeometryDescription(){}   // Constructor
	virtual ~GeometryDescription(){}  // Destructor
	virtual bool interactsWith(Element*, GeometryEntity*) = 0;
};


#endif //_GEOMETRYDESCRIPTION_H_
