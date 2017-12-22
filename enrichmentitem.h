#ifndef _ENRICHMENTITEM_H_
#define _ENRICHMENTITEM_H_

#include "FEM_Component.h"
#include "domain.h"
#include <vector>

// abstract class
class GeometryEntity; 
class EnrichmentFunction;
class Element; 
class Vertex; 
class TimeStep; 
class EnrichmentDetector;

using namespace std;

// Enrichment items : cracks, holes, material interfaces...
// Derive from FEM Component father class

class EnrichmentItem :public FEMComponent
{

public:
	EnrichmentItem(int, Domain*); // Constructor
	virtual ~EnrichmentItem();   // Destructor

	void                 getGeometry();   // get the initial geometry

	vector< EnrichmentFunction* >* giveEnrFuncVector();

	bool interactsWith(Element*);

	virtual void treatEnrichment();

	EnrichmentDetector*   defineMyEnrDetector();

	EnrichmentDetector*   giveMyEnrDetector();

	GeometryEntity*        giveMyGeo();    // return the geometry

	void                   setListOfInteractedElements(Element *e);

	vector<Element*>*      giveElementsInteractWithMe();

	EnrichmentItem*      typed();

	EnrichmentItem*      ofType(char*);

	char*                giveClassName(char* s)
	{
		return strcpy(s, "EnrichmentItem");
	}
	int                  giveNumber()
	{
		return FEMComponent::giveNumber();
	}

	virtual void     printOutputAt(TimeStep*, FILE*, FILE*){}
	/*!
	If a node is enriched by the Heaviside function H(x) and it does not satisfy
	the area inclusion criteria ( see Dolbow et al. 1999) then this node is not
	enriched by H(x) any more => matrix is not singular.
	Virtual method, only derived class CrackInterior overrides this method.
	*/
	virtual void     resolveLinearDependency(){}
	/*!
	A node should not be enriched by both H(x) and branch functions of the same
	crack. Remove H(x), just enriched by branch functions.
	Pure virtual method, only derived class CrackInterior overrides this method.
	New method coded in.
	*/
	virtual void     resolveConflictsInNodalEnrichment(){}
	/*!
	Update the geometry, virtual function that requires derived class (Crack,
	Hole,...) must implement this method.

	*/
	virtual void     UpdateMyGeometry(){ ; }
	/*!
	Update the enrichment. Just CrackTip has actual implementation since the update is
	just local around the crack tip.
	*/
	virtual void     updateEnrichment(){ ; }
	/*!

	*/
	virtual void     treatMeshGeoInteraction(Element* e){ ; }

	//=================== Methods used for debugging ==============================
	virtual void     printYourSelf(){ ; }

protected:
	GeometryEntity*      myGeometry;			//!< the original geometry
	int                  geometryID;			//!< input the ID of the element
	vector<EnrichmentFunction*>*  myEnrichFns;  //!< the enrichment functions
	EnrichmentDetector*  myEnrDetector;			//!< enrichment detector 
	vector<Element*>*    interactedElements;	//!< a list of elements interacting with me
};


#endif //_ENRICHMENTITEM_H_
