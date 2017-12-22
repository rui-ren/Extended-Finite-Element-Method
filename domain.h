//   ********************
//   *** CLASS DOMAIN ***
//   ********************


#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "list.h"
#include "newtonraphson.h"
#include "auxiliaryfield.h"
#include <stdio.h>
#include <vector>
#include <map>

class Element; class Node; class Material; class TimeIntegrationScheme;
class TimeStep; class Load; class LoadTimeFunction; class LinearSystem;
class FileReader; class Skyline;	class NewtonRaphson; class NLSolver;
// XFEM classes
class EnrichmentItem; class GeometryEntity; class EnrichmentFunction;
class CrackGrowthDirectionLaw; class CrackGrowthIncrementLaw;

//! Domain contains all of finite element objects.
class Domain
	/*!
	This class implements the domain of a finite element problem.\n
	DESCRIPTION :\n
	The domain stores its components (elements, nodes, materials, loads, load-
	time functions, enrichmentitems, enrichmentfunction,...) in lists : 'elementList',
	'nodeList', etc. Its component 'timeIntegrationScheme' needs no list, since it is unique. The domain
	also possesses the system of linear algebraic equations of the problem
	('linearSystem') and the NonLinarSolver ('nlSolver').
	The domain possesses two streams 'inputStream' and 'outputStream' that
	the components use in order to read their input and write their output in
	the data file.\n
	TASKS :

	- Managing its components and the linear system.
	The domain maintains lists of its components. The domain is responsible
	for creating these components (methods 'giveElement','giveLinearSystem',
	etc). Asking the domain is the only way for a component (e.g., an ele-
	ment) to have access to another component (e.g., its material), or to
	the data file, or to the linear system.
	- Returning the input/output streams. Since these two streams act on the
	same file (the data file), the domain is responsible for always closing
	one when the other one is to be used.
	*/
{
private:
	char*       				dataFileName;
	char*       				s00FileName;
	char*       				s01FileName;
	char*       				hisFileName;
	char*       				logFileName;
	char*       				disFileName;
	char*       				strFileName;
	List*        			elementList;
	List*        			nodeList;
	List*        			materialList;
	List*        			loadList;
	List*        			enrichmentFunctionList; //!< List of enrichment functions
	List*        			enrichmentItemList;     //!< List of enrichment item:cracks, holes...
	List*						geoEntityList;         //!< List of geometry entities of discontinnuous 
	CrackGrowthDirectionLaw* directionLaw; //!< Law of crack growth direction
	CrackGrowthIncrementLaw* incrementLaw; //!< Law of crack growth increment
	bool                  isFEM; //!< marker, standard FEM problems
	bool                  isXFEM; //!< marker, XFEM problems
	List*						loadTimeFunctionList;
	TimeIntegrationScheme*	timeIntegrationScheme;
	NLSolver*					nlSolver;
	int						   numberOfElements;//!< number of elements
	size_t					   numberOfEnrichmentFunctions;//!< number of total enrichment functions
	size_t	   			   numberOfEnrichmentItems;//!< number of enrichment items
	size_t                numberOfGeoEntities;        //!< number of geometry entities
	int						   numberOfFreeDofs;
	int						   numberOfNodes;//!< number of nodes
	FileReader*  			inputStream;
	FILE*      				outputStream;
	FloatArray*				unknownArray;
	FieldType             planeElasticity; //!< plane strain or plane stress
	std::map<Node*, std::vector<Element*> > nodalSupport; //!< the supports of all nodes in the domain
	std::vector<Node*>*   enrichedNodesList; //!< store enriched nodes just for efficiently post processing.( avoid looping on all nodes).
	bool                  isMultiMatDomain; //!< multi material domains, 2005-09-21 
	// -------------------------- post processing ------------------------------
	int            numfig; //!< number of matlab figures exported sb-2004June29
	char*          mlbFileName;//!< Matlab export file SB June17-2004
	char*          mlbFileNameFigs;//!< Matlab export file for the figures
	char*          mlbFileNameSummary;//!< Matlab export file for the summary : contains only the error for each step

public:

	Domain();                             //!< constructors
	Domain(char*);                        //!< constructors
	~Domain();                            //!< destructor

	/** @name Solving
	*  Methods for solving ...
	*/
	//@{
	void               solveYourself();
	void               solveYourselfAt(TimeStep*);
	void               formTheSystemAt(TimeStep*);
	int                giveNumberOfElements();
	void               terminate(TimeStep*);
	//@}

	/** @name Solving fracture mechanics problem
	*  Solvable problems :
	*       static unique crack, quasi-static crack growth
	*       static interfacial crack
	*       structure with holes & cracks
	*/
	//@{
	void               solveFractureMechanicsProblem();
	void               solveFractureMechanicsProblemAt(TimeStep* stepN);
	//@}

	/** @name Management of the mesh components
	*  Methods used to handle the mesh of a FE problem
	*/
	//@{
	Element*            giveElement(int);//!< return the ith element
	NLSolver*		       giveNLSolver();
	Load*               giveLoad(int);
	LoadTimeFunction*   giveLoadTimeFunction(int);
	Material*           giveMaterial(int);//!< return the ith material
	Node*               giveNode(int);//!< return the ith node
	EnrichmentItem*     giveEnrichmentItem(int);//!< return the ith enrichment item
	EnrichmentFunction* giveEnrichmentFunction(int);//!< return the ith enrichment function
	GeometryEntity*     giveGeoEntity(int);   //!< returns the geometry entity  
	FieldType           giveFieldType(); //!< plane strain or plane stress
	CrackGrowthDirectionLaw* giveGrowthDirectionLaw(); //!< return the crack growth direction law
	CrackGrowthIncrementLaw* giveGrowthIncrementLaw(); //!< return the crack growth increment law
	TimeIntegrationScheme*   giveTimeIntegrationScheme();
	void                     instanciateYourself();
	//@}

	/** @name Input / Ouput
	*  Methods for reading input file and return output files
	*/
	//@{
	char*              giveDataFileName();
	FileReader*        giveInputStream();
	FILE*              giveOutputStream();
	int                readNumberOf(char*);
	char*              giveS00FileName() { return s00FileName; };
	char*              giveS01FileName() { return s01FileName; };
	char*              giveHisFileName() { return hisFileName; };
	char*              giveLogFileName() { return logFileName; };
	char*              giveDisFileName() { return logFileName; };
	char*              giveStrFileName() { return logFileName; };
	//@}

	/** @name Functional Functions
	*  Functions to compute the functional
	*/
	//@{
	FloatArray*		GiveInitialGuess();
	FloatArray*     GivePastUnknownArray();
	FloatArray*		ComputeFunctionalAt(FloatArray*);
	Skyline*			ComputeJacobian();
	Skyline*        computeTangentStiffnessMatrix();
	//@}

	/** @name Global numbering Functions
	*  global numbering
	*/
	//@{
	int				giveNumberOfNodes();
	int				giveNumberOfFreeDofs();
	//@}

	//solution array which is stored here!!!
	FloatArray* giveUnknownArray() { return unknownArray; }

	//to compute initial acceleration
	Skyline*     giveInitialStiffnessMatrix();
	Skyline*     giveInitialMassMatrix();
	FloatArray*  giveInitialLoadVector();
	FloatArray*  giveInitialDisplacementVector();

	size_t       giveNumberOfEnrichmentItems(); //!< returns total number of enrichment items
	size_t       giveNumberOfEnrichmentFunctions(); //!< returns total number of enrichment functions
	size_t       giveNumberOfGeoEntities();     //!< returns total number of geometry entities
	void         increaseNumberOfEnrichmentItems(){ numberOfEnrichmentItems++; }
	void         increaseNumberOfEnrichmentFunctions(){ numberOfEnrichmentFunctions++; }

	/*!
	Find out elements interacting with the discontinuities.
	Modified in 2005-09-07 as following:
	The geo-mesh interaction procedure is divided into two steps :
	1. Find out elements interacting with all enrichment items except CrackTips.
	2. Find out elements interacting with CrackTips based on set of elements interacting
	with their associated CrackInteriors.
	For example, considering two two-tips cracks, i.e., we have 6 EnrichmentItems
	With the previous implementation,i.e., loop on all elements(1000 eles) for each EnrichmentItems
	then we must do 6*1000 computations.
	Now, with this improved way, we just do about 2*1000 + ... computations.
	Furthermore, this way allows the tips to be on element edge, the crack can be align with
	element edges.
	*/
	void         treatMeshGeoInteractionPhase1();
	void         treatMeshGeoInteractionPhase2();
	/*!
	Set enriched nodes for the whole domain.
	Do this at each time step once all geometries have been updated at the end
	of the previous step.
	*/
	void         treatEnrichment();
	/*!
	Resolve conflicts in enrichment of all nodes in the mesh
	Improved version : 2005-09-11
	*/
	void         resolveConflictsInEnrichment();
	/*!
	Treat the linear dependency.
	Nodes not be satisfied the area inclusion criteria will be removed from
	the list of H(x) enriched nodes. See N.Moes et al. (1999).
	*/
	void         resolveLinearDependencyForEnrichment();
	/*!
	Branched cracks problem
	Procedure to find nodes enriched by junction J(x).
	Start: 12-11-2005.
	*/
	void         treatEnrichmentForBranchedCracks();

	/*!
	Build the supports for all nodes in the mesh
	*/
	std::map<Node*, std::vector<Element*> >  buildNodalSupports();
	std::map<Node*, std::vector<Element*> >  giveNodalSupports();

	/** @name Post processing functions
	*  Export data to Matlab for post processing
	*  SB2004-06-17
	*/
	//@{
	void exportMatlabPlotCommandsUndeformedMesh(char* filename);
	void exportMatlabPlotCommands(char* filename);
	void exportMatlabHeader(char* filename);
	void exportMeshToMatlab(char* filename);
	void exportCrackGeoToMatlab(char* filename);         // NVP 2005-07-18
	void exportEnrichedNodesToMatlab(char* filename);    // NVP 2005-09-02
	void exportDeformedMeshToMatlab(char* filename);
	void exportNodeResultsToMatlab(char* filename);
	void exportElementResultsToMatlab(char* filename);
	void exportGaussPointsToMatlab(char* filename);
	void exportConnectivityToMatlab(char* filename);
	//@}

	/*!
	Insert new enriched \c Node aNode into "enrichedNodesList" of the domain.
	Whenever a node is enriched, this method is called.
	*/
	void  setEnrichedNodesList(Node* aNode);
	void  removeNodeFromEnrichedNodesList(Node* aNode);
	/*!
	Returns true if XFEM and false otherwise.
	2005-09-05
	*/
	bool  isXFEMorFEM();
	/*!
	Return true if the receiver composed of multi materials
	Implemented 2005-09-21 to solve edge crack growth in bimaterial plate
	*/
	bool  isMultiMaterialDomain();
	void  checkMultiMatDomain(); // called after list of materials defined !!!
};


#endif // _DOMAIN_H_
