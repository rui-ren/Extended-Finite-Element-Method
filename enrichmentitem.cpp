//*************************************************************************************/
#include "enrichmentitem.h"
#include "geometryentity.h"
#include "FEM_Component.h"
#include "enrichmentdetector.h"
#include "enrichmentfunction.h"
#include "crackinterior.h"
#include "cracktip.h"
#include "materialinterface.h"
#include "hole.h"
#include "domain.h"
#include "element.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <typeinfo>
#include <iostream>
using namespace std;

// constructor EnrichmentItem
EnrichmentItem::EnrichmentItem(int n, Domain* aDomain)   // class domain contain all finite element objects.
: FEMComponent(n, aDomain)								 // default aDomain
{
	// initialize the myGeometry, geometryID, myEnrDetector, myEnrichFns, interactedElements
	myGeometry = NULL;	
	geometryID = 0;	
	myEnrDetector = NULL;
	myEnrichFns = NULL;
	interactedElements = NULL;
}

//destructor EnrichmentItem
EnrichmentItem ::  ~EnrichmentItem()
{
	delete myGeometry;
	delete myEnrichFns;
	delete interactedElements;
	delete myEnrDetector;
}

void EnrichmentItem::getGeometry()
// Reads the information of geometry from the input file
{
	if (geometryID == 0)
		geometryID = this->readInteger("geometry");
	myGeometry = domain->giveGeoEntity(geometryID);
}

GeometryEntity* EnrichmentItem::giveMyGeo()
// Displace the geometry
{
	if (myGeometry == NULL)
		this->getGeometry();   // read the geometry
	return myGeometry;
}

bool EnrichmentItem::interactsWith(Element *elem)
{
	myGeometry = this->giveMyGeo();
	return myGeometry->interactsWith(elem);
}

EnrichmentItem*  EnrichmentItem::ofType(char* aClass)
//check the type of enrichment function
// Returns a new enrichment item, which has the same number than the receiver,
// but belongs to aClass (Crack, or Hole).
{
	EnrichmentItem* newEnrichItem;
	if (!strcmp(aClass, "CrackInterior"))			// check the CrackInterior
		newEnrichItem = new CrackInterior(number, domain);
	else if (!strcmp(aClass, "CrackTip"))			// check the CrackTip
		newEnrichItem = new CrackTip(number, domain);
	else if (!strcmp(aClass, "MaterialInterface"))
		newEnrichItem = new MaterialInterface(number, domain);
	else if (!strcmp(aClass, "Hole"))
		newEnrichItem = new Hole(number, domain);
	else {
		cout << aClass<<endl;
		throw logic_error("Out of the EnrichmentItem");
	}
	return newEnrichItem;
}

EnrichmentItem*  EnrichmentItem::typed()
// Returns a new enrichment item, which has the same number than the receiver,
// but belongs to aClass (Crack, or Hole).
{
	EnrichmentItem* newEnrichItem;
	char     type[32];

	this->readString("class", type);
	newEnrichItem = this->ofType(type);

	return newEnrichItem;
}

vector<EnrichmentFunction*>* EnrichmentItem::giveEnrFuncVector()
// Reads from the input file the enrichment functions used to
//	model this enrichment item
{
	EnrichmentFunction *enrFn;

	if (!myEnrichFns)
	{
		myEnrichFns = new vector<EnrichmentFunction*>;
		size_t size = this->readInteger("EnrichmentFuncs");
		for (size_t i = 0; i < size; i++)
		{
			enrFn = domain->giveEnrichmentFunction(this->readInteger("EnrichmentFuncs", i + 2));
			enrFn->setMyEnrichmentItem(this); // enrFn used to model this enr. item
			myEnrichFns->push_back(enrFn);
		}
	}

	return myEnrichFns;
}

EnrichmentDetector* EnrichmentItem::defineMyEnrDetector()
// Returns the EnrichmentDetector of the receiver. Creates it it does not exist  yet 
// Convention : 
//      enrichScheme = 1 : for CrackTip, enriche all nodes of tip-element
//      enrichScheme = 2 : for CrackTip, enriche all nodes within a ball of radius r
//      enrichScheme = 3 : for CrackInterior, enriche all nodes of split-element
//      enrichScheme = 4 : for biofilms applications. 
{
	size_t enrDetectorID = this->readInteger("enrichScheme");

	EnrichmentDetector* myEnrDetector;

	switch (enrDetectorID)
	{
	case 1:
		myEnrDetector = new ElementAroundPointDetect;
		break;
	case 2:
		myEnrDetector = new BallAroundPointDetect;
		break;
	case 3:
		myEnrDetector = new SplitElementDetect;
		break;
	case 4:
		myEnrDetector = new BandOfWidthDetect;
		break;
	}

	return myEnrDetector;
}

EnrichmentDetector* EnrichmentItem::giveMyEnrDetector()
{
	if (!myEnrDetector)
		myEnrDetector = this->defineMyEnrDetector();

	return myEnrDetector;
}

void EnrichmentItem::setListOfInteractedElements(Element *e)
// ***********************************************************
{
	if (!interactedElements)
		interactedElements = new std::vector<Element*>;

	interactedElements->push_back(e);
}

vector<Element*>* EnrichmentItem::giveElementsInteractWithMe()
{
	return interactedElements;
}

void EnrichmentItem::treatEnrichment()
// *************************************
// get the corresponding enrichment detector
// This enrichment detector how to set enriched nodes ...
{
	EnrichmentDetector* myEnrDetector = this->giveMyEnrDetector();
	myEnrDetector->setEnrichedNodes(this);

	// if call Element::conflicts(Circle*), need to reset Element::checked = false for the next time
	size_t enrDetectorID = this->readInteger("enrichScheme");
	if (enrDetectorID == 2)
	{
		for (size_t i = 0; i < domain->giveNumberOfElements(); i++)
			domain->giveElement(i + 1)->clearChecked();
	}
}

