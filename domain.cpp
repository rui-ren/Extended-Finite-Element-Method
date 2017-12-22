//   file DOMAIN.CPP

#include "domain.h"
#include "element.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "material.h"
#include "load.h"
#include "enrichmentitem.h"
#include "crackinterior.h"
#include "enrichmentfunction.h"
#include "geometryentity.h"
#include "piecewiselinear.h"
#include "loadtime.h"
#include "timinteg.h"
#include "linsyst.h"
#include "newtonraphson.h"
#include "string.h"
#include "freader.h"
#include "delaunay.h"
#include "crackgrowthdirectionlaw.h"
#include "crackgrowthincrementlaw.h"
#include "cracktip.h"
#include "vertex.h"
#include "clock.h"
#include <typeinfo>
#include <fstream> // filestream flux sur les fichiers
#include <iostream>

Domain::Domain()
// Constructor. Creates a new domain.
{
	dataFileName = NULL;
	elementList = new List(0);
	nodeList = new List(0);
	materialList = new List(0);
	loadList = new List(0);
	loadTimeFunctionList = new List(0);
	enrichmentFunctionList = new List(0); // for XFEM
	enrichmentItemList = new List(0); // for XFEM
	geoEntityList = new List(0); // for XFEM
	directionLaw = NULL;                   // for XFEM
	incrementLaw = NULL;                   // for XFEM
	isFEM = false;                         // for XFEM
	isXFEM = false;                        // for XFEM
	timeIntegrationScheme = NULL;
	nlSolver = NULL;
	inputStream = NULL;
	outputStream = NULL;
	numberOfElements = 0;
	numberOfNodes = 0;
	numberOfFreeDofs = 0;
	numberOfEnrichmentFunctions = 0;        // for XFEM
	numberOfEnrichmentItems = 0;			// for XFEM
	numberOfGeoEntities = 0;				// for XFEM
	unknownArray = NULL;
	planeElasticity = PlaneStress;
	numfig = 0;								//number of matlab figures 
	enrichedNodesList = NULL;				// list of enriched nodes, just for quick post-processing
	isMultiMatDomain = false;
	//numberOfInteractedElements = 0 ;
}


Domain::Domain(char* s)
// Constructor. Creates a new domain with data file 's'.
{
	dataFileName = new char[strlen(s) + 1];
	strcpy(dataFileName, s);

	elementList = new List(0);
	nodeList = new List(0);
	materialList = new List(0);
	loadList = new List(0);
	loadTimeFunctionList = new List(0);
	enrichmentFunctionList = new List(0);	// for XFEM
	enrichmentItemList = new List(0);		// for XFEM
	geoEntityList = new List(0);			// for XFEM
	directionLaw = NULL;                   // for XFEM
	incrementLaw = NULL;                   // for XFEM
	isFEM = false;
	isXFEM = false;
	timeIntegrationScheme = NULL;
	nlSolver = NULL;
	inputStream = NULL;
	outputStream = NULL;
	numberOfElements = 0;
	numberOfNodes = 0;
	numberOfEnrichmentFunctions = 0;
	numberOfEnrichmentItems = 0;
	numberOfFreeDofs = 0;
	unknownArray = NULL;
	planeElasticity = PlaneStress;
	enrichedNodesList = NULL;
}

Domain :: ~Domain()
// Destructor.
{
	delete dataFileName;
	delete s00FileName;
	delete s01FileName;
	delete hisFileName;
	delete logFileName;
	delete disFileName;
	delete strFileName;
	delete elementList;
	delete nodeList;
	delete materialList;
	delete enrichmentFunctionList; // XFEM
	delete enrichmentItemList;     // XFEM
	delete directionLaw;
	delete incrementLaw;
	delete nlSolver;
	delete loadList;
	delete loadTimeFunctionList;
	delete timeIntegrationScheme;
	delete unknownArray;
	delete inputStream;
	delete enrichedNodesList;
	if (outputStream)
		fclose(outputStream);
}

Skyline*  Domain::computeTangentStiffnessMatrix()
// Computes the tangent stiffness matrix of the receiver submitted to the
// displacements d.
{
	Element     *elem;
	FloatMatrix *k;
	IntArray    *loc;
	Skyline     *stiffnessMatrix;

	stiffnessMatrix = new Skyline();
	const int noElems = this->giveNumberOfElements();

	for (size_t i = 0; i < noElems; i++)
	{
		elem = this->giveElement(i + 1);
		k = elem->computeLhsAt(timeIntegrationScheme->giveCurrentStep());
		loc = elem->giveLocationArray();
		stiffnessMatrix->assemble(k, loc);
		elem->reinitializeStiffnessMatrix();
		delete k;
	}
	return stiffnessMatrix;
}

Skyline*  Domain::giveInitialStiffnessMatrix()
{
	Element     *elem;
	FloatMatrix *k;
	IntArray    *loc;
	Skyline     *answer;

	answer = new Skyline();

	for (size_t i = 0; i < this->giveNumberOfElements(); i++)
	{
		elem = this->giveElement(i + 1);
		k = elem->GiveStiffnessMatrix();
		loc = elem->giveLocationArray();
		answer->assemble(k, loc);
		delete k;
	}

	return answer;
}

Skyline*  Domain::giveInitialMassMatrix()
{
	Element     *elem;
	FloatMatrix *k;
	IntArray    *loc;
	Skyline     *answer;

	answer = new Skyline();

	for (size_t i = 0; i < this->giveNumberOfElements(); i++)
	{
		elem = this->giveElement(i + 1);
		k = elem->giveMassMatrix();
		loc = elem->giveLocationArray();
		answer->assemble(k, loc);
	}

	return answer;
}


void  Domain::formTheSystemAt(TimeStep* stepN)
// Assembles the system of linear equations, at the current time step.
// No assembly here, this will be made in NLSolver!!! - SC - 25.7.97
{
	this->giveNumberOfElements();
	for (size_t i = 0; i < numberOfElements; i++)
		//this -> giveElement(i) -> assembleYourselfAt(stepN) ;
		this->giveElement(i + 1);

	int nNodes = this->readNumberOf("Node");
	for (size_t i = 0; i < nNodes; i++)
		//this -> giveNode(i) -> assembleYourLoadsAt(stepN) ;
		this->giveNode(i + 1);
}


char*  Domain::giveDataFileName()
// Returns the name of the file containing the data of the problem.
// The other files handling is also done here! - SC 08.97
{
	char s[64];
	char temp[64] = "";
	FILE *hisFile, *s00File, *s01File, *logFile, *strFile, *disFile, *mlbFile;

	if (!dataFileName) {
		printf("***************************************************************\n");
		printf("***                                                         ***\n");
		printf("***  Object Oriented Enriched Finite Elements Library       ***\n");
		printf("***          OpenXFEM++ based on FEMOBJ, 08.2005            ***\n");
		printf("***       With applications to Fracture Mechanics           ***\n");
		printf("*** (c) Nguyen Vinh Phu, Stephane Bordas and Cyrille Dunant ***\n");
		printf("***                                                         ***\n");
		printf("*************************************************************** \n");

		printf("Please enter the name of the data file: \n");
		gets(s);


		dataFileName = new char[strlen(s) + 1];
		strcpy(dataFileName, s);
		strncat(temp, dataFileName, strlen(s) - 3);
		s00FileName = new char[strlen(s) + 1];
		strcpy(s00FileName, temp);
		strcat(s00FileName, "s00");
		s01FileName = new char[strlen(s) + 1];
		strcpy(s01FileName, temp);
		strcat(s01FileName, "s01");
		s00File = fopen(s00FileName, "wb");
		fclose(s00File);
		s01File = fopen(s01FileName, "wb");
		fclose(s01File);
		//dis & str file names for ascii results, initialize files
		disFileName = new char[strlen(s) + 1];
		strcpy(disFileName, temp);
		strcat(disFileName, "dis");
		strFileName = new char[strlen(s) + 1];
		strcpy(strFileName, temp);
		strcat(strFileName, "str");
		disFile = fopen(disFileName, "w");
		fclose(disFile);
		strFile = fopen(strFileName, "w");
		fclose(strFile);
		//his file name for history, write first dummy line
		hisFileName = new char[strlen(s) + 1];
		strcpy(hisFileName, temp);
		strcat(hisFileName, "his");
		hisFile = fopen(hisFileName, "w");
		fprintf(hisFile, "           0           0           0           0\n");
		fclose(hisFile);
		logFileName = new char[strlen(s) + 1];
		strcpy(logFileName, temp);
		strcat(logFileName, "log");
		logFile = fopen(logFileName, "w");
		fclose(logFile);



		///sb2004June17
		///matlab file name for output export and postpro in MATLAB
		mlbFileName = new char[strlen(s) + 1];
		strcpy(mlbFileName, temp);
		strcat(mlbFileName, "m");
		mlbFile = fopen(mlbFileName, "w");
		fclose(mlbFile);

		///RB 2004nov.
		///matlab file name for output figures
		char temp2[64] = ""; // must declare a null string in order to use strncat
		strncat(temp2, s, strlen(s) - 4); // copy temp to mlbFileNameFigs suppressing the last character (here the dot)
		strcat(temp2, "Figs.m");
		mlbFileNameFigs = new char[strlen(s) + 5];
		strcpy(mlbFileNameFigs, temp2);
		mlbFile = fopen(mlbFileNameFigs, "w");
		fclose(mlbFile);

		char temp4[64] = ""; // must declare a null string in order to use strncat

		strncat(temp4, s, strlen(s) - 4); // copy temp to mlbFileNameFigs suppressing the last character (here the dot)
		fclose(mlbFile);
	}

	return dataFileName;
}


Element*  Domain::giveElement(int n)
// Returns the n-th element. Creates this element if it does not exist yet.
{
	Element* elem;

	if (elementList->includes(n))
		elem = (Element*)elementList->at(n);
	else {
		elem = (Element*)Element(n, this).typed();
		elementList->put(n, elem);
	}

	return elem;
}


FileReader*  Domain::giveInputStream()
// Returns an input stream on the data file of the receiver.
{
	if (inputStream)
		return inputStream->reset();

	else {
		if (outputStream) {              // flush output stream, if it exists
			fclose(outputStream);
			outputStream = NULL;
		}
		inputStream = new FileReader(this->giveDataFileName());
		return inputStream;
	}
}


Load*  Domain::giveLoad(int n)
// Returns the n-th load. Creates this load if it does not exist yet.
{
	Load* load;

	if (loadList->includes(n))
		load = (Load*)loadList->at(n);
	else {
		load = (Load*)Load(n, this).typed();
		loadList->put(n, load);
	}

	return load;
}


LoadTimeFunction*  Domain::giveLoadTimeFunction(int n)
// Returns the n-th load-time function. Creates this fuction if it does
// not exist yet.
{
	LoadTimeFunction* ltf;

	if (loadTimeFunctionList->includes(n))
		ltf = (LoadTimeFunction*)loadTimeFunctionList->at(n);
	else {
		ltf = (LoadTimeFunction*)LoadTimeFunction(n, this).typed();
		loadTimeFunctionList->put(n, ltf);
	}

	return ltf;
}


Material*  Domain::giveMaterial(int n)
// Returns the n-th material. Creates this material if it does not exist
// yet.
{
	Material* mat;

	if (materialList->includes(n))
		mat = (Material*)materialList->at(n);
	else {
		mat = Material(n, this).typed();
		materialList->put(n, mat);
	}

	return mat;
}

NLSolver*  Domain::giveNLSolver()
// Returns the nlSolver. Creates this nonlinear solver if it does not exist
// yet.
{

	if (nlSolver)
		return nlSolver;
	else {
		nlSolver = NLSolver(1, this).typed();
	}

	return nlSolver;
}

Node*  Domain::giveNode(int n)
// Returns the n-th node. Creates this node if it does not exist yet.
{
	Node *node;

	if (nodeList->includes(n))
		node = (Node*)nodeList->at(n);
	else {
		node = new Node(n, this);
		nodeList->put(n, node);
	}

	return node;
}


int  Domain::giveNumberOfElements()
// Returns the number of elements the problem consists of.
{
	if (!numberOfElements)
		numberOfElements = this->readNumberOf("Element");
	return numberOfElements;
}


FILE*  Domain::giveOutputStream()
// Returns an output stream on the data file of the receiver.
{
	if (!outputStream) {
		if (inputStream) {                // flush input stream, if it exists
			delete inputStream;
			inputStream = NULL;
		}
		outputStream = fopen(dataFileName, "a");
	}

	return outputStream;
}


TimeIntegrationScheme*  Domain::giveTimeIntegrationScheme()
// Returns the time integration algorithm. Creates it if it does not
// exist yet.
{
	TimeIntegrationScheme* scheme;

	if (timeIntegrationScheme)
		return  timeIntegrationScheme;
	else {
		scheme = TimeIntegrationScheme(1, this).typed();
		timeIntegrationScheme = scheme;
		return scheme;
	}
}


void  Domain::instanciateYourself()
// Creates all objects mentioned in the data file.
// Exception : time step 2 and subsequent ones are not instanciated.
{
	int i, n;

#  ifdef VERBOSE
	printf("Reading all data from input file \n");
#  endif

	this->giveTimeIntegrationScheme();

	//NLSolver
	this->giveNLSolver()->instanciateYourself();
	this->giveNLSolver()->giveLinearSystem()->carveYourselfFor(this);

	n = this->readNumberOf("Node");
	nodeList->growTo(n);
	for (i = 1; i <= n; i++)
		this->giveNode(i)->instanciateYourself();

	n = this->giveNumberOfElements();
	elementList->growTo(n);
	for (i = 1; i <= n; i++)
		this->giveElement(i)->instanciateYourself();

	n = this->readNumberOf("Material");
	materialList->growTo(n);
	for (i = 1; i <= n; i++)
		this->giveMaterial(i)->instanciateYourself();

	n = this->readNumberOf("Load");
	loadList->growTo(n);
	for (i = 1; i <= n; i++)
		this->giveLoad(i)->instanciateYourself();

	n = this->readNumberOf("LoadTimeFunction");
	loadTimeFunctionList->growTo(n);
	for (i = 1; i <= n; i++)
		this->giveLoadTimeFunction(i)->instanciateYourself();

}


int  Domain::readNumberOf(char* type)
// Gets from the data file the number of objects of type 'type' (e.g.
// Element, or Node) that the receiver possesses.
{
	char value[8];

	this->giveInputStream()->read(type, value);
	return atoi(value);
}


void  Domain::solveYourself()
// Solves the problem described by the receiver.
{

	isFEM = true;
	TimeStep* currentStep;

	this->giveTimeIntegrationScheme();
	while (currentStep = timeIntegrationScheme->giveNextStep())
		this->solveYourselfAt(currentStep);

}


void  Domain::solveYourselfAt(TimeStep* stepN)
// Solves the problem at the current time step.
{

	std::cout << " ****************************************************" << std::endl;
	std::cout << " ***        SOLVING THE EQUATION SYSTEM           ***" << std::endl;
	std::cout << " ****************************************************" << std::endl;

	if (unknownArray) {
		delete unknownArray;
	}

	unknownArray = this->giveNLSolver()->Solve();

	std::cout << " ****************************************************" << std::endl;
	std::cout << " ***               POST-PROCESSING                ***" << std::endl;
	std::cout << " ****************************************************" << std::endl;

	this->terminate(stepN);

}

FloatArray*  Domain::ComputeFunctionalAt(FloatArray* dxacc)
{
	Element    *elem;
	Node       *node;
	FloatArray *f;
	IntArray   *loc;

	FloatArray *answer = new FloatArray(this->giveNumberOfFreeDofs());

	for (size_t i = 0; i < this->giveNumberOfElements(); i++)
	{
		elem = this->giveElement(i + 1);
		f = elem->computeRhsAt(timeIntegrationScheme->giveCurrentStep(), dxacc);
		loc = elem->giveLocationArray();
		answer->assemble(f, loc);
		delete f;
	}

	for (size_t i = 0; i < this->giveNumberOfNodes(); i++)
	{
		node = this->giveNode(i + 1);
		f = node->ComputeLoadVectorAt(timeIntegrationScheme->giveCurrentStep());
		if (f) {
			loc = node->giveLocationArray();
			answer->assemble(f, loc);
			delete f;
		}
	}

	return answer;
}

FloatArray*  Domain::giveInitialLoadVector()
{
	Element    *elem;
	Node       *node;
	FloatArray *answer, *f;
	IntArray   *loc;

	answer = new FloatArray(this->giveNumberOfFreeDofs());

	int nElem = this->giveNumberOfElements();
	for (size_t i = 0; i < nElem; i++)
	{
		elem = this->giveElement(i + 1);
		f = elem->ComputeLoadVectorAt(timeIntegrationScheme->giveCurrentStep());
		loc = elem->giveLocationArray();
		answer->assemble(f, loc);
		delete f;
	}

	int nNodes = this->giveNumberOfNodes();
	for (size_t i = 0; i < nNodes; i++)
	{
		node = this->giveNode(i + 1);
		f = node->ComputeLoadVectorAt(timeIntegrationScheme->giveCurrentStep());
		if (f) {
			loc = node->giveLocationArray();
			answer->assemble(f, loc);
			delete f;
		}
	}

	return answer;
}

FloatArray*  Domain::giveInitialDisplacementVector()
{
	int pos;
	double val;
	FloatArray *aFloatArray;
	aFloatArray = new FloatArray(this->giveNumberOfFreeDofs());

	int nNodes = this->giveNumberOfNodes();
	for (size_t i = 1; i <= nNodes; i++)
	{
		for (size_t j = 1; j <= 2; j++)
		{
			if (this->giveNode(i)->giveDof(j)->hasBc() == 0)
			{
				val = this->giveNode(i)->giveDof(j)->giveUnknown('d', timeIntegrationScheme->giveCurrentStep());
				pos = this->giveNode(i)->giveDof(j)->giveEquationNumber();
				aFloatArray->at(pos) = val;
			}

		}
	}

	return aFloatArray;
}


Skyline*  Domain::ComputeJacobian()
// Evaluates the jacobian matrix of the functional
{
	return this->computeTangentStiffnessMatrix();
}

FloatArray*  Domain::GivePastUnknownArray()
// Returns an array with the unknowns of the previous
// step. If this is the first step, returns an array
// full of zeroes.
// Modification made by NVP 2005-09-05 for XFEM !!!
// This modification is bad !!!
{
	int nNodes, pos;
	double val;
	FloatArray *aFloatArray;
	aFloatArray = new FloatArray(this->giveNumberOfFreeDofs());

	//SC - 12.98
	int anInitialTimeStep = 0; //for Newmark
	if (timeIntegrationScheme->isStatic()) anInitialTimeStep = 1; //for Static

	if (timeIntegrationScheme->giveCurrentStep()->giveNumber() == anInitialTimeStep || this->isXFEMorFEM())
		return aFloatArray;
	else {
		nNodes = this->giveNumberOfNodes();
		for (size_t i = 1; i <= nNodes; i++)
		{
			for (size_t j = 1; j <= 2; j++)
			{
				if (this->giveNode(i)->giveDof(j)->hasBc() == 0) {
					pos = this->giveNode(i)->giveDof(j)->giveEquationNumber();
					val = 0;
					if (timeIntegrationScheme->isStatic()) {
						//d(n-1)
						val = this->giveNode(i)->giveDof(j)->givePastUnknown('d', timeIntegrationScheme->giveCurrentStep());
					}
					else if (timeIntegrationScheme->isNewmark()) {
						//predictor for D
						val = this->giveNode(i)->giveDof(j)->giveUnknown('D', timeIntegrationScheme->giveCurrentStep());
					}
					aFloatArray->at(pos) = val;
				}
			}
		}
	}

	return aFloatArray;
}

FloatArray*  Domain::GiveInitialGuess()
// Returns the starting point for the Newton-Raphson algorithm.
{
	return this->GivePastUnknownArray();
}
/*
int  Domain :: giveNumberOfFreeDofs ()
// Returns the number of unknowns of the receiver.
{
int numberOfFreeDofs = 0;

for (size_t i = 0 ; i < this->giveNumberOfNodes(); i++)
{
for (size_t j = 1 ; j <= this->giveNode(i+1)->giveNumberOfDofs() ; j++)
{
if (this->giveNode(i+1)->giveDof(j)->hasBc() == 0)
numberOfFreeDofs++;
}
}

return numberOfFreeDofs;
} */

int  Domain::giveNumberOfFreeDofs()
// ***********************************
// Returns the number of unknowns of the receiver.
// Modified version 2005-09-08.
// Why we must to compute this quantity several times ??? Rewrite so that it is computed
// only once time for each step.
// Attention with this, at the end of each step, remember to reset numberOfFreeDofs to ZERO
// But, I did :)
{
	if (numberOfFreeDofs == 0)
	{
		for (size_t i = 0; i < this->giveNumberOfNodes(); i++)
			for (size_t j = 0; j < this->giveNode(i + 1)->giveNumberOfDofs(); j++)
			{
				if (this->giveNode(i + 1)->giveDof(j + 1)->hasBc() == 0)
					numberOfFreeDofs++;
			}
	}
	// DEBUG ONLY

	return numberOfFreeDofs;
}

int  Domain::giveNumberOfNodes()
// Returns the number of nodes of the mesh.
{
	if (!numberOfNodes)
		numberOfNodes = this->readNumberOf("Node");
	return numberOfNodes;
}

EnrichmentItem*  Domain::giveEnrichmentItem(int n)
// ************************************************
// Returns the n-th enr. item. Creates this enr. item if it does not exist yet.
{
	EnrichmentItem* enrItem;

	if (enrichmentItemList->includes(n))
		enrItem = (EnrichmentItem*)enrichmentItemList->at(n);
	else {
		enrItem = (EnrichmentItem*)EnrichmentItem(n, this).typed();
		enrichmentItemList->put(n, enrItem);
	}

	return enrItem;
}

EnrichmentFunction*  Domain::giveEnrichmentFunction(int n)
// **********************************************************
// Returns the n-th enrichment function. Creates this function if it
// does not exist yet.
{
	EnrichmentFunction* enrFunction;

	if (enrichmentFunctionList->includes(n))
		enrFunction = (EnrichmentFunction*)enrichmentFunctionList->at(n);
	else {
		enrFunction = (EnrichmentFunction*)EnrichmentFunction(this, n).typed();
		enrichmentFunctionList->put(n, enrFunction);
	}

	return enrFunction;
}

GeometryEntity*  Domain::giveGeoEntity(int n)
// **********************************************
// Returns the n-th geometry entity. Creates this geo. entity
// if it does not exist yet.
{
	GeometryEntity *geo;

	if (geoEntityList->includes(n))
		geo = (GeometryEntity*)geoEntityList->at(n);
	else {
		geo = (GeometryEntity*)GeometryEntity(n, this).typed();
		geoEntityList->put(n, geo);
	}

	return geo;
}


size_t  Domain::giveNumberOfEnrichmentItems()
// ********************************************
// Returns the number of enrichment items in the domain.
{
	if (numberOfEnrichmentItems == 0)
		numberOfEnrichmentItems = this->readNumberOf("EnrichmentItem");
	return numberOfEnrichmentItems;
}

size_t  Domain::giveNumberOfEnrichmentFunctions()
// ************************************************
// Returns the number of enrichment functions in the domain.
{
	if (numberOfEnrichmentFunctions == 0)
		numberOfEnrichmentFunctions = this->readNumberOf("EnrichmentFunction");
	return numberOfEnrichmentFunctions;
}

size_t  Domain::giveNumberOfGeoEntities()
// ****************************************
// Returns the number of geometry entities in the domain.
{
	if (numberOfGeoEntities == 0)
		numberOfGeoEntities = this->readNumberOf("GeometryEntity");
	return numberOfGeoEntities;
}

FieldType  Domain::giveFieldType()
// *******************************
// Returns plane strain or plane stress
{
	//  HOW TO IMPLEMENT THIS ? JUST READ FROM DATA FILE !!!
	planeElasticity = PlaneStress;

	return planeElasticity;
}

CrackGrowthDirectionLaw*  Domain::giveGrowthDirectionLaw()
// **********************************************************
{
	if (directionLaw)
		return directionLaw;
	else
	{
		directionLaw = CrackGrowthDirectionLaw(1, this).typed();
	}

	return directionLaw;
}

CrackGrowthIncrementLaw*  Domain::giveGrowthIncrementLaw()
// **********************************************************
{

	if (incrementLaw)
		return incrementLaw;
	else
	{
		incrementLaw = CrackGrowthIncrementLaw(1, this).typed();
	}

	return incrementLaw;
}


void  Domain::terminate(TimeStep* stepN)
// ****************************************
// Performs all operations (printings and updates) for terminating time
// step stepN.
{
	TimeStep* nextStep;
	Element*  elem;
	FILE      *s00File, *s01File, *disFile, *strFile, *hisFile;

	int convergenceStatus = this->giveNLSolver()->giveConvergenceStatus();

	/* No longer necessary since all info have been exported to Matlab M files

	//his file
	hisFile = fopen(hisFileName, "a");
	fprintf(hisFile,"    0%5d%5d    .00000    1  %10.5f%5d\n",this->giveNLSolver()->giveNumberOfIterations(),
	-this->giveNLSolver()->giveConvergenceStatus(),(float)(stepN->giveNumber()),(stepN->giveNumber()));
	fclose(hisFile);


	//nodes

	disFile = fopen(disFileName, "a");
	fprintf (disFile,"\nSolution %d, # iterations = %d, convergence status = %d\n",stepN->giveNumber(),

	this->giveNLSolver()->giveNumberOfIterations(),this->giveNLSolver()->giveConvergenceStatus());
	s00File = fopen(s00FileName, "ab");

	int nNodes = nodeList -> giveSize() ;
	for (size_t i = 0 ; i < nNodes ; i++)
	this -> giveNode(i+1) -> printOutputAt(stepN, disFile, s00File) ;
	fclose(s00File);
	fclose(disFile);

	*/

	//saving for matlab
	string comment("%");
	string st(mlbFileName);
	st += "_";
	char buffer[50];
	//sprintf(stepN->giveNumber(),buffer,10); // convert the step number into a char*
	st.append(buffer);
	st += "  Convergence Status =";
	//sprintf(this->giveNLSolver()->giveConvergenceStatus(),buffer,10); // convert the convergence status into a char*
	st.append(buffer);

	comment += st + "\n";

	if (stepN->giveNumber() == 1)
	{
		this->exportMatlabHeader(mlbFileName);
	}

	if (isFEM && this->giveNLSolver()->giveConvergenceStatus() == 1) // means step has converged)
	{
		this->exportMeshToMatlab(mlbFileName);
		this->exportDeformedMeshToMatlab(mlbFileName);
		this->exportNodeResultsToMatlab(mlbFileName);
		this->exportGaussPointsToMatlab(mlbFileName);
		this->exportConnectivityToMatlab(mlbFileName);
		this->exportMatlabPlotCommands(mlbFileNameFigs);
		this->exportElementResultsToMatlab(mlbFileName);
	}
	//this->exportElementResultsToMatlab(mlbFileName);

	//elements
	//strFile = fopen(strFileName, "a");
	//fprintf (strFile,"\nSolution %d, # iterations = %d, convergence status = %d\n",stepN->giveNumber(),

	//this->giveNLSolver()->giveNumberOfIterations(),this->giveNLSolver()->giveConvergenceStatus());

	//s01File = fopen(s01FileName, "ab");

	//fclose(s01File);
	//fclose(strFile);

	if (isXFEM)
	{
		// with XFEM, the underlying mesh is fixed, then just export the mesh one time !!!
		if (stepN->giveNumber() == 1)
		{
			std::cout << " Export mesh to Matlab ... " << endl;
			this->exportMeshToMatlab(mlbFileName);
			this->exportConnectivityToMatlab(mlbFileName);
			this->exportMatlabPlotCommands(mlbFileNameFigs);
		}

		this->exportEnrichedNodesToMatlab(mlbFileName);
		this->exportCrackGeoToMatlab(mlbFileName);

		if (this->readNumberOf("ExportDisplacement") == 1)
		{
			std::cout << " Export deformed mesh to Matlab ... " << endl;
			this->exportDeformedMeshToMatlab(mlbFileName);
		}
		if (this->readNumberOf("ExportStress") == 1)
		{
			std::cout << " Export stresses at Gps to Matlab ... " << endl;
			this->exportGaussPointsToMatlab(mlbFileName);
			this->exportElementResultsToMatlab(mlbFileName);
		}

		this->checkMultiMatDomain();

		// STRESS INTENSITY FACTORS COMPUTATION
		std::cout << " STRESS INTENSITY FACTORS COMPUTATION ... " << std::endl;
		EnrichmentItem  *enrItem;
		strFile = fopen(strFileName, "a");
		s01File = fopen(s01FileName, "ab");

		for (size_t i = 0; i < numberOfEnrichmentItems; i++)
		{
			enrItem = this->giveEnrichmentItem(i + 1);
			enrItem->printOutputAt(stepN, strFile, s01File);
		}

		fclose(s01File);
		fclose(strFile);

		// UPDATE THE ENRICHMENT WHEN CRACKS GROW
		// no need update for last step
		if (stepN->isNotTheLastStep())
		{
			std::cout << " ENRICHMENT UPDATE ... " << std::endl;
			for (size_t i = 0; i < numberOfEnrichmentItems; i++)
			{
				enrItem = this->giveEnrichmentItem(i + 1);
				enrItem->updateEnrichment();
			}
		}

		// recompute Gauss points for updated elements ...
		// Compute at the end of each step, makes the code modification small.
		if (stepN->isNotTheLastStep())
		{
			for (size_t i = 0; i < numberOfElements; i++)
			{
				Element* e = this->giveElement(i + 1);
				if (e->isUpdatedElement())
					e->computeGaussPoints();
			}
		}
		// ----------------------  debug only --------------------
		this->exportCrackGeoToMatlab(mlbFileName);
		this->exportEnrichedNodesToMatlab(mlbFileName);
		// --------------------------------------------------------

		// set checked = false for all elements. Used for Element::conflicts(Circle*).
		for (size_t i = 0; i < numberOfElements; i++)
			this->giveElement(i + 1)->clearChecked();

	} // end of if(isXFEM)

	  // UPDATE ELEMENTS
	for (size_t i = 0; i < numberOfElements; i++)
	{
		elem = this->giveElement(i + 1);
		//elem -> printOutputAt(stepN, strFile, s01File) ;
		elem->updateYourself(); // erase components of Gauss Points !!!
	}
	// UPDATE NODES

	for (size_t i = 0; i < nodeList->giveSize(); i++)
		this->giveNode(i + 1)->updateYourself();

	if (stepN->isNotTheLastStep())
	{
		nextStep = new
			TimeStep(stepN->giveNumber() + 1, timeIntegrationScheme);
		delete nextStep;
	}

	//checking if the last step was a diverged one!
	if (this->giveNLSolver()->giveConvergenceStatus() == 0)
	{
		std::cout << " The step was diverged " << std::endl;
		assert(false);
	}

	timeIntegrationScheme->updateYourself();
	this->giveNLSolver()->updateYourself();

	numberOfFreeDofs = 0; // reset the number of free dofs
}

void Domain::solveFractureMechanicsProblem()
// *******************************************
// Solves the problem described by the receiver.
{
	isXFEM = true;

	TimeStep* currentStep;
	this->giveTimeIntegrationScheme();

	while (currentStep = timeIntegrationScheme->giveNextStep())
		this->solveFractureMechanicsProblemAt(currentStep);
}

void  Domain::solveFractureMechanicsProblemAt(TimeStep* stepN)
// ************************************************************
// Solves the problem at the current time step for fracture mechanics
// 1. First step :
//      + mesh-geo interaction
//      + set enriched nodes
//      + resolve conflicts in enrichment ( nodes can not be enriched by both H(x) and branch functions)
//      + resolve linear dependency for step enriched nodes
//      + solve the system of equation : Ku = f
// 2. From step 2 ( quasi-static crack growth)
//      + update the geometry of cracks
//      + update the enrichment
//      + solve the updated system of equation : Ku = f
{
	if (stepN->giveNumber() == 1) // Step 1
	{

		// do the mesh geo interaction
		std::cout << " ****************************************************" << std::endl;
		std::cout << " ***      DOING THE MESH GEOMETRY INTERACTION     ***" << std::endl;
		std::cout << " ****************************************************" << std::endl;
		this->treatMeshGeoInteractionPhase1();
		this->treatMeshGeoInteractionPhase2();

		// set enriched nodes
		std::cout << " ****************************************************" << std::endl;
		std::cout << " ***            SET ENRICHED NODES                ***" << std::endl;
		std::cout << " ****************************************************" << std::endl;
		this->treatEnrichment();
		this->treatEnrichmentForBranchedCracks(); // 12-11-2005

												  // resolve conflicts in enrichment of nodes
		std::cout << " ****************************************************" << std::endl;
		std::cout << " ***       RESOLVE CONFLICTS IN ENRICHMENT        ***" << std::endl;
		std::cout << " ****************************************************" << std::endl;

		this->resolveConflictsInEnrichment(); // new !!! 2005-09-11

											  // ONLY NEED FOR INTERFACE CRACK PROBLEMS !!!
											  // so that the Delaunay triangulation code work
											  // the CrackInterior and the MaterialInterface coincide => colinear points !!!
		for (size_t i = 0; i < numberOfElements; i++)
		{
			Element* e = this->giveElement(i + 1);
			if (e->isSplitElement())
				e->resolveConflictsInEnrItems();
		}

		// resolve linear dependencies for enrichment
		std::cout << " ****************************************************" << std::endl;
		std::cout << " ***   RESOLVE LINEAR DEPENDENCIES IN ENRICHMENT  ***" << std::endl;
		std::cout << " ****************************************************" << std::endl;
		this->resolveLinearDependencyForEnrichment();
	}

	// DEBUG ONLY FOR CHECK MATERIALS

	FILE       *mlbFile;
	char        mlbFileName[15] = "checkMaterial"; // name of the file
	strcat(mlbFileName, ".dat");          // extension .m ( Matlab M file)
	mlbFile = fopen(mlbFileName, "w");

	Element* elt;
	Material *mat;
	std::vector<size_t> matIDs;
	std::cout << "CHECKING MATERIALS" << endl;

	for (size_t i = 0; i < numberOfElements; i++)
	{
		elt = this->giveElement(i + 1);
		//matIDs = elt->giveMatIDs();
		mat = elt->giveMaterial();
		fprintf(mlbFile, "%d ", elt->giveNumber());
		if (elt->containsMultiMats() == false)
			fprintf(mlbFile, "%d ", mat->giveNumber());
		else
			fprintf(mlbFile, " 2 Vat Lieu !!!");
		fprintf(mlbFile, "\n");
	}

	fclose(mlbFile);
	// END DEBUG FOR CHECK MATERIALS*/

	/*// CHECK THE COMPUTATION OF THE ABSOLUTE OF SIGNED DISTANCE

	FILE       *mlbFile;
	char        mlbFileName[15] = "signed" ; // name of the file
	strcat(mlbFileName,".m");          // extension .m ( Matlab M file)
	mlbFile = fopen(mlbFileName, "w");

	EnrichmentItem *enr = this->giveEnrichmentItem(1);

	std::vector<Element*> *elems = enr->giveElementsInteractWithMe();

	std::vector<EnrichmentFunction*> *step = enr->giveEnrFuncVector();
	(*step)[0]->findActiveEnrichmentItem(enr);

	fprintf (mlbFile,"splitelems = [ ");
	Node *node ;

	for (size_t i = 0 ; i < elems->size() ; i++)
	{
	size_t numnode = (*elems)[i]->giveNumberOfNodes();
	for (size_t n = 0 ; n < numnode ; n++)
	{
	node = (*elems)[i]->giveNode(n+1);
	fprintf(mlbFile,"%d ",node->giveNumber());
	if(n == numnode-1)
	fprintf(mlbFile,"\n");
	}
	}

	fprintf (mlbFile," ];\n ");

	fprintf (mlbFile,"node = [ \n");
	double val1, val2;

	for (size_t i = 0 ; i < this->giveNumberOfNodes() ; i++)
	{
	val1 = this->giveNode(i+1)->giveCoordinate(1);
	val2 = this->giveNode(i+1)->giveCoordinate(2);
	fprintf (mlbFile,"%15.14e %15.14e \n",val1, val2);
	}

	fprintf (mlbFile," ];\n ");

	fprintf (mlbFile,"stepfunc = [ ");
	double val ;

	for (size_t i = 0 ; i < elems->size() ; i++)
	{
	size_t numnode = (*elems)[i]->giveNumberOfNodes();
	for (size_t n = 0 ; n < numnode ; n++)
	{
	node = (*elems)[i]->giveNode(n+1);
	val = (*step)[0]->EvaluateYourSelfAt((*elems)[i],node);
	fprintf(mlbFile,"%3.1e ",val);
	if(n == numnode-1)
	fprintf(mlbFile,"\n");
	}
	}

	fprintf (mlbFile," ];\n ");
	fprintf (mlbFile,"figure;\n");
	fprintf (mlbFile,"hold on\n");
	fprintf (mlbFile,"plot_field(node,splitelems,'Q4',stepfunc);\n");
	fprintf (mlbFile,"plot_mesh(node,splitelems,'Q4','k-');\n");
	fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
	fprintf (mlbFile,"  crFig = plot(crack(1,1:2:size(crack,2)-1),crack(1,2:2:size(crack,2)),'k');\n");
	fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
	fprintf (mlbFile,"end\n");

	fclose (mlbFile) ;

	// CHECK THE COMPUTATION OF THE ABSOLUTE OF SIGNED DISTANCE*/


	if (unknownArray)      // delete unknowns of previous step
	{
		delete unknownArray;
	}

	// solving the system of equations
	std::cout << " ****************************************************" << std::endl;
	std::cout << " ***        SOLVING THE EQUATION SYSTEM           ***" << std::endl;
	std::cout << " ****************************************************" << std::endl;

	unknownArray = this->giveNLSolver()->Solve();

	std::cout << std::endl;
	std::cout << " Well, finish the solution of system of equation       " << std::endl;
	std::cout << " Oh, we're getting results !!!  " << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	/*// DEBUG ONLY FOR LEVEL SET
	for(size_t i = 0 ; i < numberOfElements ; i++)
	{
	Element* e = this->giveElement(i+1);
	if(e->isEnriched())
	{
	std::cout<< e->giveNumber()<<endl;
	for(size_t j = 0 ; j < e->giveNumberOfNodes() ; j++)
	std::cout<< e->giveNode(j+1)->giveLevelSet(this->giveEnrichmentItem(1))<<endl;
	std::cout<<endl;
	}
	}*/

	std::cout << " ****************************************************" << std::endl;
	std::cout << " ***               POST-PROCESSING                ***" << std::endl;
	std::cout << " ****************************************************" << std::endl;

	this->terminate(stepN);

	std::cout << "                     THE END                         " << std::endl;
	std::cout << " ****************************************************" << std::endl;
	std::cout << std::endl;
}

void Domain::treatMeshGeoInteractionPhase1()
// *******************************************
// doing the phase 1 mesh-geometry interaction to find out
// elements interacted with the discontinuities
// Phase 1 : find elements intersecting with crack interiors not with tips
// Just for the first step. For crack growth, using updated geometry and ...
// 2005-09-07
{
	Element *e;
	for (size_t i = 0; i < this->giveNumberOfElements(); i++)
	{
		e = this->giveElement(i + 1);
		e->treatGeoMeshInteraction();
	}

	// CRACK NUCLEATION
	// if the crack completely falls inside ONE element
	// then need more check to find tip element
	/*if(numberOfInteractedElements == 0)
	{
	std::cout << " CRACK NUCLEATION !!! " << endl;
	EnrichmentItem *enrItem ;
	for (size_t i = 0 ; i < this -> giveNumberOfEnrichmentItems() ; i++)
	{
	enrItem = this -> giveEnrichmentItem(i+1);
	if( typeid(*enrItem) == typeid(CrackTip) )
	{
	CrackTip *tip = (dynamic_cast<CrackTip*>(enrItem));
	for (size_t i = 0 ; i < this -> giveNumberOfElements() ; i++)
	{
	e = this -> giveElement(i+1);
	if ( tip->interactsWith(e) )
	{
	e->isEnrichedWith(tip) ;
	tip->setListOfInteractedElements(e);
	e->setEnrichmentForMyNodes(tip);
	}
	}
	}
	if( typeid(*enrItem) == typeid(CrackInterior) )
	{
	CrackInterior *cr = (dynamic_cast<CrackInterior*>(enrItem));
	cr->giveMyTips();
	}
	}
	}*/
}

void Domain::treatMeshGeoInteractionPhase2()
// *******************************************
// doing the phase 2 mesh-geometry interaction to find out
// elements interacted with the discontinuities
// 2005-09-07
{

	EnrichmentItem *enrItem;
	for (size_t i = 0; i < this->giveNumberOfEnrichmentItems(); i++)
	{
		enrItem = this->giveEnrichmentItem(i + 1);
		if (typeid(*enrItem) == typeid(CrackInterior))
			(dynamic_cast<CrackInterior*>(enrItem))->treatMeshGeoInteractionForMyTips();
	}

}

void Domain::treatEnrichment()
// *****************************
// set enriched nodes , done after doing the mesh-geo interaction
{
	for (size_t i = 0; i < this->giveNumberOfEnrichmentItems(); i++)
		this->giveEnrichmentItem(i + 1)->treatEnrichment();
}

void Domain::treatEnrichmentForBranchedCracks()
// **********************************************
// Loop on CrackInterior objects, if it is a branched crack,
// compute the intersection point between itself and its master crack, say point A
// Detecting element containing point A, nodes of this element will be enriched by
// the junction J(x).
// 12-11-2005
{
	EnrichmentItem *enrItem;
	for (size_t i = 0; i < this->giveNumberOfEnrichmentItems(); i++)
	{
		enrItem = this->giveEnrichmentItem(i + 1);
		if (typeid(*enrItem) == typeid(CrackInterior))
		{
			if ((dynamic_cast<CrackInterior*>(enrItem))->IsaBranchedCrack())
			{
				(dynamic_cast<CrackInterior*>(enrItem))->treatJunctionEnrichment();
			}
		}
	}
}

void Domain::resolveConflictsInEnrichment()
// ******************************************
// A node should not be enriched by both H(x) and branch functions of the SAME
// CRACK. Remove H(x), just enriched by branch functions.
// New method coded in 2005-09-11.
{
	for (size_t i = 0; i < this->giveNumberOfEnrichmentItems(); i++)
		this->giveEnrichmentItem(i + 1)->resolveConflictsInNodalEnrichment();
}

void Domain::resolveLinearDependencyForEnrichment()
// **************************************************
// Remove node in H(x) enriched node list if this node is not satisfied
// area inclusion creterion. See Daux et al. 2000 for details.
// 2006-06-06
{
	for (size_t i = 0; i < this->giveNumberOfEnrichmentItems(); i++)
		this->giveEnrichmentItem(i + 1)->resolveLinearDependency();
}


std::map<Node*, vector<Element*> > Domain::buildNodalSupports()
// *************************************************************
// build the support for all nodes in the domain
// Idea of Raphael Bonnaz
{
	for (size_t i = 0; i < this->giveNumberOfElements(); i++)
	{
		Element *elem = this->giveElement(i + 1);
		for (size_t j = 0; j < elem->giveNumberOfNodes(); j++)
		{
			Node *aNode = elem->giveNode(j + 1);
			nodalSupport[aNode].push_back(elem);
		}
	}
	return nodalSupport;
}

std::map<Node*, vector<Element*> > Domain::giveNodalSupports()
// ************************************************************
{
	if (nodalSupport.size() == 0)
	{
		nodalSupport = this->buildNodalSupports();
	}
	return nodalSupport;
}

void Domain::exportMatlabHeader(char* filename)
{
	FILE  	 *mlbFile;

	mlbFile = fopen(filename, "w");
	fprintf(mlbFile, "clear all;\n");
	fprintf(mlbFile, "close all;\n");
	fclose(mlbFile);

}

void Domain::exportMatlabPlotCommandsUndeformedMesh(char* filename)
// *******************************************************************
{
	FILE       *mlbFile;
	mlbFile = fopen(filename, "a");
	fprintf(mlbFile, "figure;\n");
	fprintf(mlbFile, "plot(node(:,1),node(:,2),'o');\n");
	fprintf(mlbFile, "axis equal;\n");
	fclose(mlbFile);
}

void Domain::exportMatlabPlotCommands(char* filename)
// ****************************************************
// Original version : Stephane Bordas and Raphael Bonnaz
// Modified version for XFEM, von Mises stress contour plot : Vinh Phu Nguyen 2005-07-20
{
	FILE       *mlbFile;
	mlbFile = fopen(filename, "a");

	// PLOT DISPLACEMENT RESULTS AT EACH STEP
	size_t numSteps = this->giveTimeIntegrationScheme()->giveNumberOfSteps();
	size_t i = this->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber();
	numfig++;

	// PLOT UNDEFORMED MESH
	// used convention for Matlab routine plot_mesh()
	// elementType = 'T3'(3 node triangle elements),'Q4'(4 node quadrilateral elements)...

	fprintf(mlbFile, "switch size(element,2)\n");
	fprintf(mlbFile, "  case 3 \n");
	fprintf(mlbFile, "  elementType = 'T3';\n");
	fprintf(mlbFile, "  case 6\n");
	fprintf(mlbFile, "  elementType = 'T6';\n");
	fprintf(mlbFile, "  case 4 \n");
	fprintf(mlbFile, "  elementType = 'Q4';\n");
	fprintf(mlbFile, "  case 8 \n");
	fprintf(mlbFile, "  elementType = 'Q8';\n");
	fprintf(mlbFile, "  otherwise \n");
	fprintf(mlbFile, "  error('Grr. Unknown element type');\n");
	fprintf(mlbFile, "end\n\n");

	// PLOT UNDERFORMED CONFIGURATION + CRACKS + ENRICHED NODES
	fprintf(mlbFile, "figure;\n");
	fprintf(mlbFile, "set(gcf,'Color','white');\n");
	fprintf(mlbFile, "plot_mesh(node,element,elementType,'b-');\n");

	// plot crack geometry, enriched nodes
	if (isXFEM)
	{
		// crack geometry
		fprintf(mlbFile, "hold on\n");
		fprintf(mlbFile, "for i = 1 : size(crack,1)\n");
		fprintf(mlbFile, "  crFig = plot(crack(i,1:2:size(crack,2)),crack(i,2:2:size(crack,2)),'r*-');\n");
		fprintf(mlbFile, "  set(crFig,'LineWidth',1.5)\n");
		fprintf(mlbFile, "end\n");

		// enriched nodes
		fprintf(mlbFile, "h1 = plot(StepEnrichedNodes(:,1),StepEnrichedNodes(:,2),'ro');\n");
		fprintf(mlbFile, "set(h1,'MarkerSize',10);\n");
		fprintf(mlbFile, "h2 = plot(TipEnrichedNodes(:,1),TipEnrichedNodes(:,2),'bs');\n");
		fprintf(mlbFile, "set(h2,'MarkerSize',10);\n");
		//fprintf (mlbFile,"h3 = plot(JunctionEnrichedNodes(:,1),TipEnrichedNodes(:,2),'cd');\n");
		//fprintf (mlbFile,"set(h3,'MarkerSize',10);\n");
		fprintf(mlbFile, "h3 = plot(InterfaceEnrichedNodes(:,1),InterfaceEnrichedNodes(:,2),'r*');\n");
		fprintf(mlbFile, "set(h3,'MarkerSize',10);\n");
		fprintf(mlbFile, "legend([h1,h2,h3,],'Interior','Tip','Interface')\n");

		fprintf(mlbFile, "plot(XGP(:,1),XGP(:,2),'b*');\n"); // check Gauss points
	}
	fprintf(mlbFile, "axis off ;\n");
	// use Ben Hinkle's facilities to export Matlab plots to EPS files.
	fprintf(mlbFile, "ops = struct('bounds','tight');\n");
	fprintf(mlbFile, "exportfig(gcf,'undeformed.eps',ops,'color','rgb');\n");
	fprintf(mlbFile, "\n");


	// PLOT DEFORMED CONFIGURATION
	fprintf(mlbFile, "figure;\n");
	fprintf(mlbFile, "set(gcf,'Color','white');\n");
	fprintf(mlbFile, "maxX = max(node(:,1));\n");
	fprintf(mlbFile, "minX = min(node(:,1));\n");
	fprintf(mlbFile, "L = maxX - minX;\n"); // valid for only rectangle domain !!!
	fprintf(mlbFile, "dispNorm=L/max(sqrt(U_%d(:,1).^2+U_%d(:,2).^2)); \n", i, i);
	fprintf(mlbFile, "scaleFact=0.1*dispNorm;\n");
	fprintf(mlbFile, "plot_mesh(node+scaleFact*[U_%d(:,1) U_%d(:,2)],element,elementType,'k-');\n", i, i);
	fprintf(mlbFile, "axis off ;\n");
	fprintf(mlbFile, "title('Deformed configuration');\n");
	fprintf(mlbFile, "exportfig(gcf,'deformed.eps',ops);\n"); // export to EPS file
	fprintf(mlbFile, "\n");

	// DISPLACEMENT CONTOUR PLOTS ON DEFORMED MESH

	//UX
	numfig++;
	fprintf(mlbFile, "figure;\n");
	fprintf(mlbFile, "plot_field(node+scaleFact*[U_%d(:,1) U_%d(:,2)],element,elementType,U_%d(:,1));\n", i, i, i);
	//fprintf(mlbFile,"hold on\n");
	//fprintf(mlbFile,"plot_mesh(X+scaleFact*[U_%d(:,1) U_%d(:,2)],element,elementType,'g-');\n",i,i);
	fprintf(mlbFile, "colorbar;\n");
	fprintf(mlbFile, "axis off ;\n");
	fprintf(mlbFile, "title('X-displacement');\n");
	fprintf(mlbFile, "exportfig(gcf,'ux.eps',ops,'color','rgb');\n"); // export to EPS file

																	  //UY
	numfig++;
	fprintf(mlbFile, "figure;\n");
	fprintf(mlbFile, "plot_field(node+scaleFact*[U_%d(:,1) U_%d(:,2)],element,elementType,U_%d(:,2));\n", i, i, i);
	//fprintf(mlbFile,"hold on\n");
	//fprintf(mlbFile,"plot_mesh(X+scaleFact*[U_%d(:,1) U_%d(:,2)],node,elementType,'g-');\n",i,i);
	fprintf(mlbFile, "colorbar;\n");
	fprintf(mlbFile, "axis off ;\n");
	fprintf(mlbFile, "title('Y-displacement');\n");
	fprintf(mlbFile, "exportfig(gcf,'uy.eps',ops,'color','rgb');\n"); // export to EPS file
	fprintf(mlbFile, "\n");

	// PLOT STRESS CONTOURS AT GAUSS POINTS
	// From stresses at Gauss points, triangulating them using Delaunay function
	// using the tecplotout of Chessa to export to Tecplot software : STRESS.DAT .

	fprintf(mlbFile, "tri = DELAUNAY(XGP(:,1),XGP(:,2));\n");
	fprintf(mlbFile, "varname{1}='X';\n");
	fprintf(mlbFile, "varname{2}='Y';\n");
	fprintf(mlbFile, "varname{3}='sigma_xx';\n");
	fprintf(mlbFile, "varname{4}='sigma_yy';\n");
	fprintf(mlbFile, "varname{5}='sigma_xy';\n");
	fprintf(mlbFile, "varname{6}='sigma_vonMises';\n");
	fprintf(mlbFile, "le = size(sigma_1,2);\n");
	fprintf(mlbFile, "sigma=[sigma_1(:,1:5:le) sigma_1(:,2:5:le) sigma_1(:,3:5:le) sigma_1(:,5:5:le)];\n");
	fprintf(mlbFile, "value=[XGP(:,1) XGP(:,2) sigma];\n");
	fprintf(mlbFile, "tecplotout(tri,'T3',XGP,varname,value,'stress.dat');\n\n");
	int stepNumber =
		this->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber();
	for (int compo = 1; compo <= 5; compo++)
	{
		numfig++;
		fprintf(mlbFile, "figure;\n");

		if (compo == 1)
		{
			fprintf(mlbFile, "plot_field(XGP,tri,'T3',sigma_%d(:,1:5:size(sigma_%d,2)));\n", stepNumber, stepNumber);
			fprintf(mlbFile, "title('\\sigma_{xx}'); \n");
			fprintf(mlbFile, "stressFile = 'sigmaXX.eps' ; \n");
		}
		if (compo == 2)
		{
			fprintf(mlbFile, "plot_field(XGP,tri,'T3',sigma_%d(:,2:5:size(sigma_%d,2)));\n", stepNumber, stepNumber);
			fprintf(mlbFile, "title('\\sigma_{yy}'); \n");
			fprintf(mlbFile, "stressFile = 'sigmaYY.eps' ; \n");
		}
		if (compo == 3)
		{
			fprintf(mlbFile, "plot_field(XGP,tri,'T3',sigma_%d(:,3:5:size(sigma_%d,2)));\n", stepNumber, stepNumber);
			fprintf(mlbFile, "title('\\sigma_{xy}'); \n");
			fprintf(mlbFile, "stressFile = 'sigmaXY.eps' ; \n");
		}
		if (compo == 4) {
			fprintf(mlbFile, "title('\\sigma_{zz}'); \n");
			fprintf(mlbFile, "stressFile = 'sigmaZZ.eps' ; \n");
		}
		if (compo == 5)
		{
			fprintf(mlbFile, "plot_field(XGP,tri,'T3',sigma_%d(:,5:5:size(sigma_%d,2)));\n", stepNumber, stepNumber);
			fprintf(mlbFile, "title('Equivalent von Mises stress'); \n");
			fprintf(mlbFile, "stressFile = 'vonMises.eps' ; \n");
		}

		fprintf(mlbFile, "colorbar;\n");
		fprintf(mlbFile, "axis off ;\n");
		fprintf(mlbFile, "exportfig(gcf,stressFile,ops,'color','rgb');\n\n"); // export to EPS file
	}

	fclose(mlbFile);
}

void Domain::exportMeshToMatlab(char* filename)
// **********************************************
// Write the Matlab file containing the coord of nodes of the mesh
// with format :
//      X = [ x_node1 y_node1
//            x_node2 y_node2
//            ...
//            x_noden y_noden ]
{
	FILE       *mlbFile;
	mlbFile = fopen(filename, "a");
	int nNodes = nodeList->giveSize();
	fprintf(mlbFile, "node = [ \n");
	double val1 = 0.0, val2 = 0.0;

	for (size_t i = 0; i < nNodes; i++)
	{
		val1 = this->giveNode(i + 1)->giveCoordinate(1);
		val2 = this->giveNode(i + 1)->giveCoordinate(2);
		fprintf(mlbFile, "%15.14e %15.14e \n", val1, val2);
	}

	fprintf(mlbFile, " ];\n ");

	fclose(mlbFile);
}

void Domain::exportEnrichedNodesToMatlab(char* filename)
// *******************************************************
// Export enriched nodes to Matlab for plotting ...
{
	FILE    *mlbFile;
	mlbFile = fopen(filename, "a");

	// Step-enriched nodes

	fprintf(mlbFile, "StepEnrichedNodes = [ \n");
	double x, y;
	for (size_t i = 0; i < enrichedNodesList->size(); i++)
	{
		if ((*enrichedNodesList)[i]->isStepEnriched())
		{
			x = (*enrichedNodesList)[i]->giveCoordinate(1);
			y = (*enrichedNodesList)[i]->giveCoordinate(2);
			fprintf(mlbFile, "%15.14e %15.14e \n", x, y);
		}
	}
	fprintf(mlbFile, " ];\n ");

	// Tip-enriched nodes

	fprintf(mlbFile, "TipEnrichedNodes = [ \n");

	for (size_t i = 0; i < enrichedNodesList->size(); i++)
	{
		if ((*enrichedNodesList)[i]->isTipEnriched())
		{
			x = (*enrichedNodesList)[i]->giveCoordinate(1);
			y = (*enrichedNodesList)[i]->giveCoordinate(2);
			fprintf(mlbFile, "%15.14e %15.14e \n", x, y);
		}
	}
	fprintf(mlbFile, " ];\n ");

	// Junction-enriched nodes

	fprintf(mlbFile, "JunctionEnrichedNodes = [ \n");

	for (size_t i = 0; i < enrichedNodesList->size(); i++)
	{
		if ((*enrichedNodesList)[i]->isJunctionEnriched())
		{
			x = (*enrichedNodesList)[i]->giveCoordinate(1);
			y = (*enrichedNodesList)[i]->giveCoordinate(2);
			fprintf(mlbFile, "%15.14e %15.14e \n", x, y);
		}
	}
	fprintf(mlbFile, " ];\n ");

	// Interface-enriched nodes

	fprintf(mlbFile, "InterfaceEnrichedNodes = [ \n");

	for (size_t i = 0; i < enrichedNodesList->size(); i++)
	{
		if ((*enrichedNodesList)[i]->isInterfaceEnriched())
		{
			x = (*enrichedNodesList)[i]->giveCoordinate(1);
			y = (*enrichedNodesList)[i]->giveCoordinate(2);
			fprintf(mlbFile, "%15.14e %15.14e \n", x, y);
		}
	}
	fprintf(mlbFile, " ];\n ");

	fclose(mlbFile);
}

void Domain::exportCrackGeoToMatlab(char* filename)
// **************************************************
// Write the Matlab file contaning the crack geometry
// with format :
// crack = [x1 y1 x2 y2 .... xn yn  % of the first crack
//         ...]

{
	string crack("crack");
	string bracket(" = [\n");
	char stepnum[32];

	string theString = crack + bracket; //string to store the geo of crack

	GeometryEntity *geoEntity;
	for (size_t i = 0; i < this->giveNumberOfGeoEntities(); i++)
	{
		geoEntity = this->giveGeoEntity(i + 1);
		geoEntity->exportToMatlab(theString);

	}

	string bracket2("];\n \n");
	theString += bracket2;

	ofstream destination(filename, ios::out | ios::app);
	destination << theString;
	destination.close();
}

void Domain::exportDeformedMeshToMatlab(char* filename)
// *******************************************************
{
	FILE       *mlbFile;
	TimeStep* timeStep =
		this->giveTimeIntegrationScheme()->giveCurrentStep();
	int stepNumber =
		this->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber();
	mlbFile = fopen(filename, "a");

	// DISPLACEMENT FIELD
	fprintf(mlbFile, "U_%d = [ \n", stepNumber);
	double d1 = 0.0, d2 = 0.0;

	for (size_t i = 0; i < nodeList->giveSize(); i++)
	{
		d1 = this->giveNode(i + 1)->giveDof(1)->giveUnknown('d', timeStep); // u_i
		d2 = this->giveNode(i + 1)->giveDof(2)->giveUnknown('d', timeStep); // v_i
		fprintf(mlbFile, "%15.14e %15.14e \n", d1, d2);
	}
	fprintf(mlbFile, " ];\n ");

	fclose(mlbFile);
}

void Domain::exportNodeResultsToMatlab(char* filename)
// *****************************************************
{
	TimeStep* timeStep =
		this->giveTimeIntegrationScheme()->giveCurrentStep();
	int stepNumber =
		this->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber();
	FILE       *mlbFile;

	mlbFile = fopen(filename, "a");

	int nNodes = nodeList->giveSize();
	fprintf(mlbFile, "u1_%d = [ \n", stepNumber);
	double val = 0.0;

	for (size_t i = 0; i < nNodes; i++)
	{
		val = this->giveNode(i + 1)->giveDof(1)->giveUnknown('d', timeStep);
		fprintf(mlbFile, "%15.14e \n", val);
	}

	fprintf(mlbFile, " ];\n ");

	fprintf(mlbFile, "u2_%d = [ \n", stepNumber);
	for (size_t i = 0; i < nNodes; i++)
	{
		val = this->giveNode(i + 1)->giveDof(2)->giveUnknown('d', timeStep);
		fprintf(mlbFile, "%15.14e \n", val);
	}

	fprintf(mlbFile, " ];\n ");
	fclose(mlbFile);
}

void Domain::exportElementResultsToMatlab(char* filename)
// ********************************************************
{

	TimeStep* timeStep =
		this->giveTimeIntegrationScheme()->giveCurrentStep();
	int stepNumber =
		this->giveTimeIntegrationScheme()->giveCurrentStep()->giveNumber();

	// ------------------------------------------------------------------------------
	//STRESSES at Gauss points
	Element* elt; //elt used in loop
	string sigma("sigma_");
	string bracket(" = [ \n");
	char stepnum[32];
	//_itoa(stepNumber,stepnum,10);
	sigma += stepnum;
	string theString = sigma + bracket; //string to store the stresses at the elements

	for (size_t i = 0; i < elementList->giveSize(); i++)
	{
		elt = this->giveElement(i + 1);
		elt->exportStressResultsToMatlab(theString);
	}
	string bracket2("];\n \n");
	theString += bracket2;
	// ------------------------------------------------------------------------------


	//STRAINS
	string epsilon("epsilon_");
	epsilon += stepnum;
	theString += (epsilon + bracket); //string to store the stresses at the elements
	for (size_t i = 0; i < elementList->giveSize(); i++)
	{
		elt = this->giveElement(i + 1);
		elt->computeStrain(timeStep);
		elt->exportStrainResultsToMatlab(theString);
	}

	theString += bracket2;

	ofstream destination(filename, ios::out | ios::app);
	destination << theString;
	destination.close();
}

void Domain::exportGaussPointsToMatlab(char* filename)
// *****************************************************
// Exports the global coordinates of Gauss points.
{
	Element* elt;       //elt used in loop
	string xgp("XGP");
	string bracket(" = [ \n");
	string theString = xgp + bracket; //string to store the GPs of the elements

	for (size_t i = 0; i < this->giveNumberOfElements(); i++) {
		elt = this->giveElement(i + 1);
		elt->exportGaussPointsToMatlab(theString);
	}
	string bracket2("];\n \n");
	theString += bracket2;
	ofstream destination(filename, ios::out | ios::app);
	destination << theString;
	destination.close();
}


void Domain::exportConnectivityToMatlab(char* filename)
// ******************************************************
// connectivity = [1  2  34      // connectivity of element 1(T3 element)
//                 12 23 10      // connectivity of element 2(T3 element)
//                 ...]
{
	FILE       *mlbFile;
	mlbFile = fopen(filename, "a");

	fprintf(mlbFile, "nbelt = %d \n", this->giveNumberOfElements());
	fprintf(mlbFile, "element = [ \n");

	Element* elt;
	int numnode;
	Node* node;

	for (size_t e = 0; e < this->giveNumberOfElements(); e++)
	{
		elt = this->giveElement(e + 1);
		numnode = elt->giveNumberOfNodes();
		for (size_t n = 1; n <= numnode; n++)
		{
			node = elt->giveNode(n);
			fprintf(mlbFile, "%d ", node->giveNumber());
			if (n == numnode)
				fprintf(mlbFile, "\n");
		}
	}

	fprintf(mlbFile, "];\n\n");
	fclose(mlbFile);
}

void Domain::setEnrichedNodesList(Node* aNode)
// *********************************************
{
	if (enrichedNodesList == NULL)
		enrichedNodesList = new std::vector<Node*>;

	if (find(enrichedNodesList->begin(), enrichedNodesList->end(), aNode)
		== enrichedNodesList->end())
		enrichedNodesList->push_back(aNode);

}
void Domain::removeNodeFromEnrichedNodesList(Node* aNode)
// ********************************************************
{
	remove(enrichedNodesList->begin(), enrichedNodesList->end(), aNode);
}

bool Domain::isXFEMorFEM()
// *************************
{
	return (isXFEM) ? true : false;
}

bool Domain::isMultiMaterialDomain()
// *********************************
{
	return isMultiMatDomain;
}

void Domain::checkMultiMatDomain()
// ********************************
{
	// is multi material problem ?

	size_t nMats = materialList->giveSize();
	isMultiMatDomain = (nMats != 1) ? true : false;
}
