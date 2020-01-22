#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <G4Polyhedra.hh>
#include <G4Trd.hh>
#include <G4Cons.hh>
#include <G4Torus.hh>
#include <G4UnionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <globals.hh>

#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using std::vector;
using std::stringstream;
using std::max;

#include "Xenon100LXeSensitiveDetector.hh"
#include "Xenon100PmtSensitiveDetector.hh"

#include "Xenon100DetectorConstruction.hh"

map<G4String, G4double> Xenon100DetectorConstruction::m_hGeometryParameters;

Xenon100DetectorConstruction::Xenon100DetectorConstruction()
{
	m_pRotationXPlus225 = new G4RotationMatrix();
	m_pRotationXPlus225->rotateX(22.5*deg);

	m_pRotationXMinus225 = new G4RotationMatrix();
	m_pRotationXMinus225->rotateX(-22.5*deg);

	m_pRotationXPlus45 = new G4RotationMatrix();
	m_pRotationXPlus45->rotateX(45.*deg);

	m_pRotationXMinus45 = new G4RotationMatrix();
	m_pRotationXMinus45->rotateX(-45.*deg);

	m_pRotationXPlus90 = new G4RotationMatrix();
	m_pRotationXPlus90->rotateX(90.*deg);

	m_pRotationXMinus90 = new G4RotationMatrix();
	m_pRotationXMinus90->rotateX(-90.*deg);

	m_pRotationYPlus90 = new G4RotationMatrix();
	m_pRotationYPlus90->rotateY(90.*deg);

	m_pRotationYMinus90 = new G4RotationMatrix();
	m_pRotationYMinus90->rotateY(-90.*deg);

	m_pRotationXPlus90YMinus90 = new G4RotationMatrix();
	m_pRotationXPlus90YMinus90->rotateX(90.*deg);
	m_pRotationXPlus90YMinus90->rotateY(-90.*deg);

	m_pRotationX180 = new G4RotationMatrix();
	m_pRotationX180->rotateX(180.*deg);

}

Xenon100DetectorConstruction::~Xenon100DetectorConstruction()
{
}

G4VPhysicalVolume*
Xenon100DetectorConstruction::Construct()
{	
//	DefineColours();
	
	DefineMaterials();

	DefineGeometryParameters();

	ConstructLaboratory();

    ConstructShield();

	ConstructXenon();

	ConstructBell();

	ConstructFieldCage();

	ConstructTPC();

	ConstructPmtArrays();

	ConstructCryostat();

	PrintGeometryInformation();

	return m_pLabPhysicalVolume;
}

G4double
Xenon100DetectorConstruction::GetGeometryParameter(const char *szParameter)
{
	return m_hGeometryParameters[szParameter];
}

//void
//Xenon100DetectorConstruction::DefineColours()
//{
//	//////////////////////////////////////////////////////////
//	// DEFINE COLOURS							//			//
//	G4Colour hLeadColour(0.5, 0.5, 0.5);		// grey		//
//	G4Colour hPolyColour(1.0,	1.0, 1.0);		// white	//
//	G4Colour hCopperColour(1.0, 1.0, 0.0);		// yellow	//
//	G4Colour hSteelColour(0.0, 0.0, 1.0);		// xlblue	//
//	G4Colour hTeflonColour(1.0, 0.0, 1.0);		// magenta	//
//	G4Colour hCirlexColour(0.7, 0.4, 0.1);		// brown	//
//	G4Colour hLXeColour(0.5, 0.0, 0.0);			// xlred	//
//	G4Colour hGXeColour(0.0, 0.75, 0.0);		// lgreen	//
//	G4Colour hPMTColour(1.0, 0.0, 0.0);			// red		//
//	G4Colour hQuartzColour(0.0, 0.5, 0.0);		// xlgreen	//
//	G4Colour hGridColour(0.5, 0.5, 0.5);		// grey		//
//	G4Colour hPhotoCathodeColour(1.0, 0.5, 0.0);// orange	//
//	//////////////////////////////////////////////////////////
//}

void
Xenon100DetectorConstruction::DefineMaterials()
{

	//================================== elements ===================================
	G4Element *Xe = new G4Element("Xenon",     "Xe", 54., 131.293	*g/mole);
	G4Element *H  = new G4Element("Hydrogen",  "H",  1.,  1.0079	*g/mole);
	G4Element *C  = new G4Element("Carbon",    "C",  6.,  12.011	*g/mole);
	G4Element *N  = new G4Element("Nitrogen",  "N",  7.,  14.007	*g/mole);
	G4Element *O  = new G4Element("Oxygen",    "O",  8.,  15.999	*g/mole);
	G4Element *F  = new G4Element("Fluorine",  "F",  9.,  18.998	*g/mole);
	G4Element *Na = new G4Element("Sodium",    "Na", 11., 22.9897	*g/mole);
	G4Element *Al = new G4Element("Aluminium", "Al", 13., 26.982	*g/mole);
	G4Element *Si = new G4Element("Silicon",   "Si", 14., 28.086	*g/mole);
	G4Element *P  = new G4Element("Fosforus",  "P",  15., 30.9738	*g/mole);
	G4Element *S  = new G4Element("Sulphur",   "S",  16., 32.065	*g/mole);
	G4Element *Ti = new G4Element("Titanium",  "Ti", 22., 47.867	*g/mole);
	G4Element *Cr = new G4Element("Chromium",  "Cr", 24., 51.996	*g/mole);
	G4Element *Mn = new G4Element("Manganese", "Mn", 25., 54.938	*g/mole);
	G4Element *Fe = new G4Element("Iron",      "Fe", 26., 55.85		*g/mole);
	G4Element *Ni = new G4Element("Nickel",    "Ni", 28., 58.693	*g/mole);
	G4Element *Cu = new G4Element("Copper",    "Cu", 29., 63.546	*g/mole);
	G4Element *Mo = new G4Element("Molybdenum","Mo", 42., 95.94		*g/mole);
	G4Element *Pb = new G4Element("Lead",      "Pb", 82., 207.2		*g/mole);
	
	//================================== materials ================================== 
	G4NistManager* pNistManager = G4NistManager::Instance();

	//------------------------------------- air -------------------------------------
	pNistManager->FindOrBuildMaterial("G4_AIR");
	// density = 0.00120479 g/cm3
	//------------------------------------- water -------------------------------------
	pNistManager->FindOrBuildMaterial("G4_WATER");

	//----------------------------------- vacuum ------------------------------------
	G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
	Vacuum->AddElement(N, 0.755);
	Vacuum->AddElement(O, 0.245);
	
	//-------------------------------- liquid xenon ---------------------------------
	G4Material *LXe = new G4Material("LXe", 2.9172*g/cm3, 1, kStateLiquid, 168.15*kelvin, 1.5*atmosphere);
	LXe->AddElement(Xe, 1);

	const G4int iNbEntries = 3;

	G4double pdLXePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdLXeScintillation[]    = {0.1,     1.0,     0.1};
	G4double pdLXeRefractiveIndex[]  = {1.63,    1.61,    1.58};
	G4double pdLXeAbsorbtionLength[] = {100.*cm, 100.*cm, 100.*cm};
	G4double pdLXeScatteringLength[] = {30.*cm,  30.*cm,  30.*cm};

	const G4double dFieldQuenchingFactor = 1.00;

	G4MaterialPropertiesTable *pLXePropertiesTable = new G4MaterialPropertiesTable();

	pLXePropertiesTable->AddProperty("FASTCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("SLOWCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("RINDEX", pdLXePhotonMomentum, pdLXeRefractiveIndex, iNbEntries);
	pLXePropertiesTable->AddProperty("ABSLENGTH", pdLXePhotonMomentum, pdLXeAbsorbtionLength, iNbEntries);
	pLXePropertiesTable->AddProperty("RAYLEIGH", pdLXePhotonMomentum, pdLXeScatteringLength, iNbEntries);

	pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", dFieldQuenchingFactor*1./(21.6*eV));
	pLXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
	pLXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pLXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pLXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);

//    LXe->SetMaterialPropertiesTable(pLXePropertiesTable);

	//-------------------------------- gaseous xenon --------------------------------
	G4Material *GXe = new G4Material("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere);
	GXe->AddElement(Xe, 1);

	G4double pdGXePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdGXeScintillation[]    = {0.1,     1.0,     0.1};
	G4double pdGXeRefractiveIndex[]  = {1.00,    1.00,    1.00};
	G4double pdGXeAbsorbtionLength[] = {100*m,   100*m,   100*m};
	G4double pdGXeScatteringLength[] = {100*m,   100*m,   100*m};

	G4MaterialPropertiesTable *pGXePropertiesTable = new G4MaterialPropertiesTable();

	pGXePropertiesTable->AddProperty("FASTCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("SLOWCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("RINDEX", pdGXePhotonMomentum, pdGXeRefractiveIndex, iNbEntries);
	pGXePropertiesTable->AddProperty("ABSLENGTH", pdGXePhotonMomentum, pdGXeAbsorbtionLength, iNbEntries);
	pGXePropertiesTable->AddProperty("RAYLEIGH", pdGXePhotonMomentum, pdGXeScatteringLength, iNbEntries);

	pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1./(21.6*eV));
	pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 10000./(keV));
	pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
	pGXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pGXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pGXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);

//    GXe->SetMaterialPropertiesTable(pGXePropertiesTable);	

	//----------------------------------- quartz ------------------------------------
	// ref: http://www.sciner.com/Opticsland/FS.htm
	G4Material *Quartz = new G4Material("Quartz", 2.201*g/cm3, 2, kStateSolid, 168.15*kelvin, 1.5*atmosphere);
	Quartz->AddElement(Si, 1);
	Quartz->AddElement(O, 2);

	G4double pdQuartzPhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdQuartzRefractiveIndex[]  = {1.50,    1.56,    1.60};
	G4double pdQuartzAbsorbtionLength[] = {30*m,    30*m,    30*m};

	G4MaterialPropertiesTable *pQuartzPropertiesTable = new G4MaterialPropertiesTable();

	pQuartzPropertiesTable->AddProperty("RINDEX", pdQuartzPhotonMomentum, pdQuartzRefractiveIndex, iNbEntries);
	pQuartzPropertiesTable->AddProperty("ABSLENGTH", pdQuartzPhotonMomentum, pdQuartzAbsorbtionLength, iNbEntries);

	Quartz->SetMaterialPropertiesTable(pQuartzPropertiesTable);

	//------------------------------- stainless steel -------------------------------
	G4Material *SS304LSteel = new G4Material("SS304LSteel", 8.00*g/cm3, 10, kStateSolid);
	SS304LSteel->AddElement(C,  0.0008);
	SS304LSteel->AddElement(Si, 0.01);
	SS304LSteel->AddElement(Mn, 0.02);
	SS304LSteel->AddElement(P,  0.00045);
	SS304LSteel->AddElement(S,  0.0003);
	SS304LSteel->AddElement(Ni, 0.12);
	SS304LSteel->AddElement(Cr, 0.17);
	SS304LSteel->AddElement(Mo, 0.025);
	SS304LSteel->AddElement(Ti, 0.004);
	SS304LSteel->AddElement(Fe, 0.64945);
	

	G4double pdSteelPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdSteelReflectivity[]   = {0.15,    0.2,     0.15};
	G4MaterialPropertiesTable *pSteelPropertiesTable = new G4MaterialPropertiesTable();

	pSteelPropertiesTable->AddProperty("REFLECTIVITY", pdSteelPhotonMomentum, pdSteelReflectivity, iNbEntries);
	SS304LSteel->SetMaterialPropertiesTable(pSteelPropertiesTable);

	//---------------------------------- aluminium ----------------------------------
	pNistManager->FindOrBuildMaterial("G4_Al");

	//---------------------------- photocathode aluminium ---------------------------
	G4Material *PhotoCathodeAluminium = new G4Material("PhotoCathodeAluminium", 8.00*g/cm3, 1, kStateSolid);
	PhotoCathodeAluminium->AddElement(Al, 1);

	G4double pdPhotoCathodePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdPhotoCathodeRefractiveIndex[]  = {1.50,    1.56,    1.60};
	G4double pdPhotoCathodeAbsorbtionLength[] = {1.*nm,   1.*nm,   1.*nm};

	G4MaterialPropertiesTable *pPhotoCathodePropertiesTable = new G4MaterialPropertiesTable();

	pPhotoCathodePropertiesTable->AddProperty("RINDEX", pdPhotoCathodePhotonMomentum, pdPhotoCathodeRefractiveIndex, iNbEntries);
	pPhotoCathodePropertiesTable->AddProperty("ABSLENGTH", pdPhotoCathodePhotonMomentum, pdPhotoCathodeAbsorbtionLength, iNbEntries);

	PhotoCathodeAluminium->SetMaterialPropertiesTable(pPhotoCathodePropertiesTable);

	//-------------------------------------------------------------------------------
	// cirlex
	G4Material *Cirlex = new G4Material("Cirlex", 1.43*g/cm3, 4, kStateSolid);
	Cirlex->AddElement(C, 22);
	Cirlex->AddElement(H, 10);
	Cirlex->AddElement(N, 2);
	Cirlex->AddElement(O, 5);

	//-------------------------------------------------------------------------------
	// ceramics
	G4Material *Ceramics = new G4Material("Ceramics", 1.00*g/cm3, 4, kStateSolid);
	Ceramics->AddElement(Na, 1);
	Ceramics->AddElement(Al, 1);
	Ceramics->AddElement(Si, 1);
	Ceramics->AddElement(O, 2);

	//----------------------------- grid mesh aluminium------------------------------
	G4Material *GridMeshAluminium = new G4Material("GridMeshAluminium", 0.174*g/cm3, 1, kStateSolid);
	GridMeshAluminium->AddElement(Al, 1);
	
	G4double pdGridMeshPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double *pdGridMeshRefractiveIndex = pdLXeRefractiveIndex;
	G4double dAbsortionLength = 2.10*mm; // exp(-GridMeshThickness/2.10) = 0.94
	G4double pdGridMeshAbsortionLength[] = {dAbsortionLength, dAbsortionLength, dAbsortionLength};
	
	G4MaterialPropertiesTable *pGridMeshPropertiesTable = new G4MaterialPropertiesTable();

	pGridMeshPropertiesTable->AddProperty("RINDEX", pdGridMeshPhotonMomentum, pdGridMeshRefractiveIndex, iNbEntries);
	pGridMeshPropertiesTable->AddProperty("ABSLENGTH", pdGridMeshPhotonMomentum, pdGridMeshAbsortionLength, iNbEntries);
	GridMeshAluminium->SetMaterialPropertiesTable(pGridMeshPropertiesTable);

	//------------------------------------ teflon -----------------------------------
	G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
	Teflon->AddElement(C, 0.240183);
	Teflon->AddElement(F, 0.759817);

	G4double pdTeflonPhotonMomentum[]  = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdTeflonRefractiveIndex[] = {1.63,    1.61,    1.58};
	G4double pdTeflonReflectivity[]    = {0.95,    0.95,    0.95};
	G4double pdTeflonSpecularLobe[]    = {0.01,    0.01,    0.01};
	G4double pdTeflonSpecularSpike[]   = {0.01,    0.01,    0.01};
	G4double pdTeflonBackscatter[]     = {0.01,    0.01,    0.01};
	G4double pdTeflonEfficiency[]      = {1.0,     1.0,     1.0};
	G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

	pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, iNbEntries);
	pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, iNbEntries);
	pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, iNbEntries);
	pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, iNbEntries);
	Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

	//------------------------------------- lead ------------------------------------
	G4Material *Lead = new G4Material("Lead", 11.34*g/cm3, 1);
	Lead->AddElement(Pb, 1); 

	//--------------------------------- polyethylene --------------------------------
	G4Material *Polyethylene = new G4Material("Polyethylene", 0.94*g/cm3, 2, kStateSolid);
	Polyethylene->AddElement(C, 2);
	Polyethylene->AddElement(H, 4);	

	//------------------------------------ copper -----------------------------------
	G4Material *Copper = new G4Material("Copper", 8.92*g/cm3, 1);
	Copper->AddElement(Cu, 1); 

	G4double pdCopperPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdCopperReflectivity[]   = {0.15,    0.2,     0.15};
	G4MaterialPropertiesTable *pCopperPropertiesTable = new G4MaterialPropertiesTable();

	pCopperPropertiesTable->AddProperty("REFLECTIVITY", pdCopperPhotonMomentum, pdCopperReflectivity, iNbEntries);
	Copper->SetMaterialPropertiesTable(pCopperPropertiesTable);
}

void
Xenon100DetectorConstruction::DefineGeometryParameters()
{
	// general
	m_hGeometryParameters["CryostatLidInnerWallDomeTopToVetoLiquidLevel"] = 70.*mm;

	m_hGeometryParameters["PmtToGrid"]             = 10.3*mm;
	m_hGeometryParameters["GridToAnode"]           = 5.08*mm;
	m_hGeometryParameters["AnodeToLiquidLevel"]    = 1.58*mm;
	m_hGeometryParameters["LiquidLevelToGrid"]     = 3.5*mm;
	m_hGeometryParameters["DriftLength"]           = 304.8*mm;
	m_hGeometryParameters["CathodeToPmt"]          = 12.5*mm;

	// xenon
	m_hGeometryParameters["GXeHeight"] = 87.5*mm;

	// bell
	m_hGeometryParameters["BellTopThickness"]  = 2.5*mm;
	m_hGeometryParameters["BellSideThickness"] = 1.5*mm;
	m_hGeometryParameters["BellRadius"]        = 188.275*mm;	
	m_hGeometryParameters["BellHeight"]        = 95.*mm;	

	m_hGeometryParameters["TopPmtTeflonHolderRadius"]         = 182.00*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderThickness"]      = 31.12*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderCutDepth"]       = 27.94*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderCutWidth"]       = 26.72*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderBottomCutWidth"] = 21.00*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderBottomCutAngle"] = 42.00*mm;	
	m_hGeometryParameters["TopPmtTeflonHolderToBell"]         = 38.48*mm;	

	m_hGeometryParameters["BellPmtSupportRingThickness"] = 2.5*mm;	
	m_hGeometryParameters["BellPmtSupportRingRadius"]    = 163.3*mm;	
	m_hGeometryParameters["BellPmtSupportRingHeight"]    = 9.5*mm;	

	m_hGeometryParameters["TopVetoAngleThickness"]   = 2.5*mm;	
	m_hGeometryParameters["TopVetoAngleLowerHeight"] = 12.59*mm;	
	m_hGeometryParameters["TopVetoAngleUpperHeight"] = 29.84*mm;	
	m_hGeometryParameters["TopVetoAngleLength"]      = 29.60*mm;	
	m_hGeometryParameters["TopVetoAngleWidth"]       = 27.94*mm;	

	m_hGeometryParameters["UpperSideVetoAngleThickness"] = 2.5*mm;	
	m_hGeometryParameters["UpperSideVetoAngleHeight"]    = 50.60*mm;	
	m_hGeometryParameters["UpperSideVetoAngleLength"]    = 28.27*mm;	
	m_hGeometryParameters["UpperSideVetoAngleWidth"]     = 27.94*mm;	
	
	// grids
	m_hGeometryParameters["TopGridRingRadius"]     		= 161.3*mm;
	m_hGeometryParameters["AnodeGridRingRadius"]   		= 161.3*mm;
	m_hGeometryParameters["BottomGridRingRadius"]  		= 161.3*mm; 
	m_hGeometryParameters["CathodeGridRingRadius"] 		= 152.4*mm;
	m_hGeometryParameters["GridRingWidth"]         		= 6.35*mm;
	m_hGeometryParameters["GridRingThickness"]     		= 1.5*mm;
	m_hGeometryParameters["GridMeshThickness"]     		= 0.13*mm;
	m_hGeometryParameters["GridHolderRadius"]       	= 157.6*mm;
	m_hGeometryParameters["GridHolderWidth"]        	= 20.0*mm;
	m_hGeometryParameters["GridHolderCutHeight"]    	= 2.7*mm;
    m_hGeometryParameters["CathodeToScreenMesh"] 		= 10.*mm;
	m_hGeometryParameters["CopperScreenThickness"] 		= 0.5*mm;	
	m_hGeometryParameters["CopperScreenOuterRadius"] 	= 186.*mm;	
	m_hGeometryParameters["TopPlateThickness"] 			= 3.*mm;	
	m_hGeometryParameters["TopPlateRadius"]    			= 188.275*mm;	
	m_hGeometryParameters["TopPlateWidth"]     			= 35.875*mm;	

	// tpc
	m_hGeometryParameters["TopPlateThickness"] 			= 3.*mm;	
	m_hGeometryParameters["TopPlateRadius"]    			= 188.275*mm;	
	m_hGeometryParameters["TopPlateWidth"]     			= 35.875*mm;	

	m_hGeometryParameters["TeflonPanelRadius"]    		= 152.4*mm;
	m_hGeometryParameters["TeflonPanelThickness"] 		= 6.35*mm;
	m_hGeometryParameters["TeflonPanelHeight"]    		= 297.18*mm;

	m_hGeometryParameters["TeflonRodsRadius"]  			= 4.76*mm;
	m_hGeometryParameters["TeflonRodsHeight"]  			= 297.18*mm;
	m_hGeometryParameters["TeflonRodsOffsetX"] 			= 178.59*mm;

	m_hGeometryParameters["BottomPlateThickness"]    	= 3.*mm;	
	m_hGeometryParameters["BottomPlateRadius"]       	= 168.28*mm;	
	m_hGeometryParameters["BottomPlateWidth"]        	= 15.88*mm;	
	m_hGeometryParameters["BottomPlateAnchorLength"] 	= 15.07*mm;	
	m_hGeometryParameters["BottomPlateAnchorWidth"]  	= 9.53*mm;	

	m_hGeometryParameters["UpperBasePlateThickness"]    = 3.*mm;	
	m_hGeometryParameters["UpperBasePlateRadius"]       = 160.02*mm;	
	m_hGeometryParameters["UpperBasePlateWidth"]        = 7.62*mm;	
	m_hGeometryParameters["UpperBasePlateAnchorLength"] = 12.66*mm;	
	m_hGeometryParameters["UpperBasePlateAnchorWidth"]  = 9.53*mm;	

	m_hGeometryParameters["LowerTeflonPanelRadius"]    	= 152.4*mm;
	m_hGeometryParameters["LowerTeflonPanelThickness"] 	= 6.35*mm;
	//m_hGeometryParameters["LowerTeflonPanelHeight"]    	= 40.13*mm;
 	m_hGeometryParameters["LowerTeflonPanelHeight"]    = 49.81*mm; // from Elizabeth's code, to fit with the copper screen?
	m_hGeometryParameters["LowerBasePlateThickness"]    = 2.5*mm;	
	m_hGeometryParameters["LowerBasePlateRadius"]       = 160.34*mm;	
	m_hGeometryParameters["LowerBasePlateSupportRingRadius"] = 152.82*mm;	

	m_hGeometryParameters["BottomPmtPlateThickness"] = 2.5*mm;	
	m_hGeometryParameters["BottomPmtPlateRadius"]    = 149.23*mm;	

	m_hGeometryParameters["BellSupportCylinderHeight"]    = 120.*mm;	
	m_hGeometryParameters["BellSupportCylinderRadius"]    = 4.*mm;	
	m_hGeometryParameters["BellSupportCylinderOffsetX"]   = 66.*mm;	
	m_hGeometryParameters["BellSupportCylinderOffsetY"]   = 114.4*mm;	

	m_hGeometryParameters["BottomVetoAngleThickness"] = 2.5*mm;	
	m_hGeometryParameters["BottomVetoAngleHeight"]    = 33.10*mm;	
	m_hGeometryParameters["BottomVetoAngleLength"]    = 31.75*mm;	
	m_hGeometryParameters["BottomVetoAngleWidth"]     = 27.94*mm;	

	m_hGeometryParameters["LowerSideVetoAngleThickness"]   = 2.5*mm;	
	m_hGeometryParameters["LowerSideVetoAngleHeight"]      = 31.75*mm;	
	m_hGeometryParameters["LowerSideVetoAngleShortLength"] = 14.33*mm;	
	m_hGeometryParameters["LowerSideVetoAngleLongLength"]  = 31.75*mm;	
	m_hGeometryParameters["LowerSideVetoAngleWidth"]       = 27.94*mm;	

	m_hGeometryParameters["TeflonSideVetoLiningRadius"]    = 203.2*mm;
	m_hGeometryParameters["TeflonSideVetoLiningThickness"] = 3.175*mm;

	m_hGeometryParameters["TeflonLedBlockWidth"]    	= 20.*mm;
	m_hGeometryParameters["TeflonLedBlockHeight"]   	= 10.*mm;
	m_hGeometryParameters["TeflonLedBlockDepth"]    	= 25.*mm;
	m_hGeometryParameters["TeflonLedBlockOffsetX"] 		= -68.6*mm;
	m_hGeometryParameters["TeflonLedBlockOffsetY"] 		= 123.5*mm;
	m_hGeometryParameters["TeflonLedBlockRotationZ"] 	= -35.*deg;

	m_hGeometryParameters["TeflonBottomReflectorThickness"] 	= 1.5 *mm;
	m_hGeometryParameters["TeflonBottomReflectorTopRadius"] 	= GetGeometryParameter("BottomPmtPlateRadius");
	m_hGeometryParameters["TeflonBottomReflectorBottomRadius"] 	= 150. *mm;
	m_hGeometryParameters["BottomReflectorHolderRadius"] 		= 3.2 *mm;
	m_hGeometryParameters["BottomReflectorHolderHeight"] 		= 30. *mm;
    m_hGeometryParameters["BottomReflectorHolderOffsetX"] 		= GetGeometryParameter("TeflonBottomReflectorTopRadius")-5.*mm;

	// pmts
	m_hGeometryParameters["PmtWidth"]                 = 25.4*mm;
	m_hGeometryParameters["PmtSpacing"]               = 2.03*mm;
	m_hGeometryParameters["PmtWindowWidth"]           = 25.00*mm;
	m_hGeometryParameters["PmtWindowThickness"]       = 1.50*mm;
	m_hGeometryParameters["PmtCasingWidth"]           = 25.40*mm;
	m_hGeometryParameters["PmtCasingHeight"]          = 27.00*mm;
	m_hGeometryParameters["PmtCasingThickness"]       = 0.50*mm;
	m_hGeometryParameters["PmtPhotoCathodeWidth"]     = 22.00*mm;
	m_hGeometryParameters["PmtPhotoCathodeThickness"] = 0.50*mm;
	m_hGeometryParameters["PmtBaseThickness"]         = 1.50*mm;
	m_hGeometryParameters["PmtToPmtBase"]             = 3.00*mm;
	m_hGeometryParameters["PmtBaseSpacerHeight"]      = 6.00*mm;

	m_hGeometryParameters["NbTopPmtsSixthRing"]  = 30;
	m_hGeometryParameters["NbTopPmtsFifthRing"]  = 24;
	m_hGeometryParameters["NbTopPmtsFourthRing"] = 20;
	m_hGeometryParameters["NbTopPmtsThirdRing"]  = 14;
	m_hGeometryParameters["NbTopPmtsSecondRing"] = 8;
	m_hGeometryParameters["NbTopPmtsFirstRing"]  = 2;
	m_hGeometryParameters["NbTopPmts"] = 98;

	m_hGeometryParameters["NbBottomPmtsFirstRow"]   = 4;
	m_hGeometryParameters["NbBottomPmtsSecondRow"]  = 7;
	m_hGeometryParameters["NbBottomPmtsThirdRow"]   = 9;
	m_hGeometryParameters["NbBottomPmtsFourthRow"]  = 10;
	m_hGeometryParameters["NbBottomPmtsFifthRow"]   = 10;
	m_hGeometryParameters["NbBottomPmtsSixthRow"]   = 10;
	m_hGeometryParameters["NbBottomPmtsSeventhRow"] = 10;
	m_hGeometryParameters["NbBottomPmtsEighthRow"]  = 9;
	m_hGeometryParameters["NbBottomPmtsNinthRow"]   = 7;
	m_hGeometryParameters["NbBottomPmtsTenthRow"]   = 4;
	m_hGeometryParameters["NbBottomPmts"] = 80;

	m_hGeometryParameters["NbTopVetoPmts"]       = 32;
	m_hGeometryParameters["NbBottomVetoPmts"]    = 32;

	m_hGeometryParameters["PmtSixthRingRadius"]  = 166.84*mm;
	m_hGeometryParameters["PmtFifthRingRadius"]  = 136.53*mm;
	m_hGeometryParameters["PmtFourthRingRadius"] = 106.2*mm;
	m_hGeometryParameters["PmtThirdRingRadius"]  = 75.87*mm;

	m_hGeometryParameters["VetoPmtToCryostatInnerWallDomeBottom"] = 120.*mm;

	// cryostat
	m_hGeometryParameters["CryostatVesselFlangeRadius"]        = 247.65*mm;
	m_hGeometryParameters["CryostatVesselFlangeThickness"]     = 25.4*mm;
	m_hGeometryParameters["CryostatVesselInnerWallRadius"]     = 203.2*mm;
	m_hGeometryParameters["CryostatVesselInnerWallHeight"]     = 715.772*mm;
	m_hGeometryParameters["CryostatVesselInnerWallThickness"]  = 1.5*mm;
	m_hGeometryParameters["CryostatVesselInnerWallDomeHeight"] = 76.10*mm;
	m_hGeometryParameters["CryostatVesselOuterWallRadius"]     = 228.6*mm;
	m_hGeometryParameters["CryostatVesselOuterWallThickness"]  = 1.5*mm;

	m_hGeometryParameters["CryostatLidFlangeRadius"]        = 247.65*mm;
	m_hGeometryParameters["CryostatLidFlangeThickness"]     = 25.4*mm;
	m_hGeometryParameters["CryostatLidInnerWallRadius"]     = 203.2*mm;
	m_hGeometryParameters["CryostatLidInnerWallHeight"]     = 163.83*mm;
	m_hGeometryParameters["CryostatLidInnerWallThickness"]  = 1.5*mm;
	m_hGeometryParameters["CryostatLidInnerWallDomeHeight"] = 31.29*mm;
	m_hGeometryParameters["CryostatLidOuterWallRadius"]     = 184.15*mm;
	m_hGeometryParameters["CryostatLidOuterWallThickness"]  = 1.5*mm;

	m_hGeometryParameters["CryostatLidInTubeRadius"]            	= 38.1*mm;
	//m_hGeometryParameters["CryostatLidInTubeHeight"]            	= 245.6*mm;
	m_hGeometryParameters["CryostatLidInTubeHeight"]            	= 210*mm; // not real dimension, to make it fit
	m_hGeometryParameters["CryostatLidInTubeHeightInside"]      	= 129.4*mm;
	m_hGeometryParameters["CryostatLidInTubeThickness"]         	= 1.5*mm;
	m_hGeometryParameters["CryostatLidInTubeFlangeRadius"]      	= 58.674*mm;
	m_hGeometryParameters["CryostatLidInTubeFlangeThickness"]   	= 31.75*mm;
	m_hGeometryParameters["CryostatLidInTubeOffsetX"]           	= 88.9*mm;
	m_hGeometryParameters["CryostatLidInTubeElbowLength"]       	= 123.57*mm;
	m_hGeometryParameters["CryostatLidInTubeTopLength"]         	= 75.59*mm;
	m_hGeometryParameters["CryostatLidInTubeTopFlangeRadius"]     	= 50.8*mm;
	m_hGeometryParameters["CryostatLidInTubeTopFlangeThickness"]   	= 25.4*mm;
	m_hGeometryParameters["CryostatLidInTubeJacketRadius"]      	= 50.8*mm;
	m_hGeometryParameters["CryostatLidInTubeJacketHeight"]      	= 179.73*mm;
	m_hGeometryParameters["CryostatLidInTubeOuterJacketRadius"] 	= 63.5*mm;
	m_hGeometryParameters["CryostatLidInTubeOuterJacketHeight"] 	= 129.9*mm;

	m_hGeometryParameters["CryostatLidOutTubeRadius"]            	= 38.1*mm;
	//m_hGeometryParameters["CryostatLidOutTubeHeight"]            	= 245.6*mm;
	m_hGeometryParameters["CryostatLidOutTubeHeight"]            	= 210.0*mm; // not real, to make it fit
	m_hGeometryParameters["CryostatLidOutTubeThickness"]         	= 1.5*mm;
	m_hGeometryParameters["CryostatLidOutTubeOffsetX"]           	= 88.9*mm;
	m_hGeometryParameters["CryostatLidOutTubeElbowLength"]       	= 123.57*mm;
	m_hGeometryParameters["CryostatLidOutTubeTopLength"]         	= 75.59*mm;
	m_hGeometryParameters["CryostatLidOutTubeTopFlangeRadius"]      = 50.8*mm;
	m_hGeometryParameters["CryostatLidOutTubeTopFlangeThickness"]   = 25.4*mm;
	m_hGeometryParameters["CryostatLidOutTubeTopToShieldLength"]    = 100.*mm;
	m_hGeometryParameters["CryostatLidOutTubeJacketRadius"]      	= 50.8*mm;
	m_hGeometryParameters["CryostatLidOutTubeJacketHeight"]      	= 154.6*mm;

	m_hGeometryParameters["ShieldHeight"]                	= 1950.*mm;
	m_hGeometryParameters["ShieldWidth"]                 	= 1700.*mm;
	m_hGeometryParameters["ShieldDepth"]                 	= 1700.*mm;
	m_hGeometryParameters["ShieldPolishLeadThickness"]   	= 150.*mm;
	m_hGeometryParameters["ShieldFrenchLeadThickness"]   	= 50.*mm;
	m_hGeometryParameters["ShieldPolyethyleneThickness"] 	= 200.*mm;
	
	m_hGeometryParameters["ShieldCopperThicknessThick"] 	= 50.*mm;
	m_hGeometryParameters["ShieldCopperThicknessThin"] 		= 5.*mm;

	m_hGeometryParameters["PolyethyleneSlabThickness"] 	= 250.*mm;
	
	//m_hGeometryParameters["ShieldInnerLeadThickness"] = 3.*mm;
	
	
	// derived quantities
	m_hGeometryParameters["LXeHeight"] = GetGeometryParameter("CryostatVesselInnerWallHeight")
		-GetGeometryParameter("CryostatLidInnerWallHeight")+GetGeometryParameter("CryostatVesselFlangeThickness")
		+GetGeometryParameter("CryostatLidInnerWallDomeHeight")+GetGeometryParameter("CryostatVesselInnerWallDomeHeight");

	m_hGeometryParameters["CryostatVesselInnerWallDomeRadius"] = (pow(GetGeometryParameter("CryostatVesselInnerWallRadius"),2)+pow(GetGeometryParameter("CryostatVesselInnerWallDomeHeight"),2))/(2*GetGeometryParameter("CryostatVesselInnerWallDomeHeight"));

	m_hGeometryParameters["CryostatLidInnerWallDomeRadius"] = (pow(GetGeometryParameter("CryostatVesselInnerWallRadius"),2)+pow(GetGeometryParameter("CryostatLidInnerWallDomeHeight"),2))/(2*GetGeometryParameter("CryostatLidInnerWallDomeHeight"));	

	m_hGeometryParameters["CryostatLidInTubeOffsetDifference"] = sqrt(pow(GetGeometryParameter("CryostatLidInnerWallDomeRadius"),2)-pow(GetGeometryParameter("CryostatLidInTubeOffsetX"),2))-(GetGeometryParameter("CryostatLidInnerWallDomeRadius")-GetGeometryParameter("CryostatLidInnerWallDomeHeight"));

	m_hGeometryParameters["CryostatLidInnerWallDomeTopToLiquidLevel"] = GetGeometryParameter("CryostatLidInnerWallDomeHeight")
		-GetGeometryParameter("CryostatLidInTubeOffsetDifference")+GetGeometryParameter("CryostatLidInTubeHeightInside")
		+GetGeometryParameter("BellTopThickness")+GetGeometryParameter("GXeHeight");

	m_hGeometryParameters["CryostatVesselFlangeToZero"] = GetGeometryParameter("CryostatLidInnerWallHeight")
		-GetGeometryParameter("CryostatLidInnerWallDomeHeight")
		+GetGeometryParameter("CryostatLidInnerWallDomeTopToLiquidLevel")+GetGeometryParameter("LiquidLevelToGrid");

	m_hGeometryParameters["CavityDepth"] = 	GetGeometryParameter("ShieldDepth")*0.5
											- GetGeometryParameter("ShieldPolishLeadThickness")
											- GetGeometryParameter("ShieldFrenchLeadThickness")
											- GetGeometryParameter("ShieldPolyThickness")
											- GetGeometryParameter("ShieldCopperThicknessThick");

	// verifications
//    assert(GetGeometryParameter("TeflonPanelHeight") == GetGeometryParameter("DriftLength")-GetGeometryParameter("GridRingThickness"));
}

void
Xenon100DetectorConstruction::ConstructLaboratory()
{
	//---------------------------------- laboratory ---------------------------------
	const G4double dLabHalfX = 5.0*m;
	const G4double dLabHalfY = 5.0*m;
	const G4double dLabHalfZ = 5.0*m;

	G4Box *pLabBox = new G4Box("LabBox", dLabHalfX, dLabHalfY, dLabHalfZ);
	
	G4Material *Air = G4Material::GetMaterial("G4_AIR");

	m_pLabLogicalVolume = new G4LogicalVolume(pLabBox, Air, "LabLogicalVolume", 0, 0, 0);

	m_pLabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(), m_pLabLogicalVolume, "Lab", 0, false, 0);

	m_pLabLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	m_pMotherLogicalVolume = m_pLabLogicalVolume;
}

void
Xenon100DetectorConstruction::ConstructShield()
{
	
	//=================================== shield ====================================
	const G4double dShieldHalfZ = 0.5*GetGeometryParameter("ShieldHeight");
	const G4double dShieldHalfX = 0.5*GetGeometryParameter("ShieldWidth");
	const G4double dShieldHalfY = 0.5*GetGeometryParameter("ShieldDepth");

	const G4double dPolishLeadThickness = GetGeometryParameter("ShieldPolishLeadThickness");
	const G4double dFrenchLeadThickness = GetGeometryParameter("ShieldFrenchLeadThickness");
	const G4double dInnerLeadThickness = GetGeometryParameter("ShieldInnerLeadThickness");
	const G4double dPolyethyleneThickness = GetGeometryParameter("ShieldPolyethyleneThickness");
	const G4double dCopperThicknessThick = GetGeometryParameter("ShieldCopperThicknessThick");
	const G4double dCopperThicknessThin = GetGeometryParameter("ShieldCopperThicknessThin");

	const G4double dPolishLeadHalfZ = dShieldHalfZ;
	const G4double dPolishLeadHalfX = dShieldHalfX;
	const G4double dPolishLeadHalfY = dShieldHalfY;
	const G4double dPolishLeadOffsetZ = 0.;

	const G4double dFrenchLeadHalfZ = dPolishLeadHalfZ - dPolishLeadThickness;
	const G4double dFrenchLeadHalfX = dPolishLeadHalfX - dPolishLeadThickness;
	const G4double dFrenchLeadHalfY = dPolishLeadHalfY - dPolishLeadThickness;
	const G4double dFrenchLeadOffsetZ = dPolishLeadOffsetZ;
	
	const G4double dPolyethyleneHalfZ = dFrenchLeadHalfZ - dFrenchLeadThickness;
	const G4double dPolyethyleneHalfX = dFrenchLeadHalfX - dFrenchLeadThickness;
	const G4double dPolyethyleneHalfY = dFrenchLeadHalfY - dFrenchLeadThickness;
	const G4double dPolyethyleneOffsetZ = dFrenchLeadOffsetZ;

	const G4double dCopperHalfZ = dPolyethyleneHalfZ - dPolyethyleneThickness;
	const G4double dCopperHalfX = dPolyethyleneHalfX - dPolyethyleneThickness;
	const G4double dCopperHalfY	= dPolyethyleneHalfY - dPolyethyleneThickness;
	const G4double dCopperOffsetZ = dPolyethyleneOffsetZ;

	const G4double dCavityHalfZ = dCopperHalfZ - 0.5*dCopperThicknessThick - 0.5*dCopperThicknessThin;
	const G4double dCavityHalfX = dCopperHalfX - dCopperThicknessThick;
	const G4double dCavityHalfY = dCopperHalfY  - dCopperThicknessThick;
	const G4double dCavityOffsetZ = - dCopperThicknessThick*0.5;

	const G4double dPolySlab_HalfX = dPolishLeadHalfX+20*cm;
	const G4double dPolySlab_HalfY = dPolishLeadHalfY;
	const G4double dPolySlab_HalfZ = 0.5*GetGeometryParameter("PolyethyleneSlabThickness");
	//const G4double dPolySlabOffsetX = dPolishLeadOffsetX;
	//const G4double dPolySlabOffsetY = dPolishLeadOffsetY;
	//const G4double dPolySlabOffsetZ = dPolishLeadOffsetZ-dShieldHalfZ-dPolySlab_HalfZ;
	const G4double dPolySlabOffsetX = 0.0;
	const G4double dPolySlabOffsetY = 0.0;
	const G4double dPolySlabOffsetZ = 0.0-dShieldHalfZ-dPolySlab_HalfZ;

/*	const G4double dPolyethyleneBottomPlateHalfZ = 25.0 *mm;
	const G4double dPolyethyleneBottomPlateOffsetZ = -dFrenchLeadHalfZ - dPolyethyleneBottomPlateHalfZ;
	
	const G4double dInnerLeadHalfZ = dPolyethyleneHalfZ + dInnerLeadThickness;
	const G4double dInnerLeadHalfX = dPolyethyleneHalfX + dInnerLeadThickness;
	const G4double dInnerLeadHalfY = dPolyethyleneHalfY + dInnerLeadThickness;
	const G4double dInnerLeadOffsetZ = dPolyethyleneOffsetZ;
*/	

	const G4double Water_thick 		= 0.2 *m;
	const G4double WaterTop_thick	= 0.4 *m;	

	const G4double WaterTop_HalfX	= dPolishLeadHalfX;
	const G4double WaterTop_HalfY	= dPolishLeadHalfY;
	const G4double WaterTop_HalfZ	= WaterTop_thick*0.5;
	
	const G4double WaterBack_HalfX	= dPolishLeadHalfX;
	const G4double WaterBack_HalfY	= Water_thick*0.5;
	const G4double WaterBack_HalfZ	= dPolishLeadHalfZ;

	const G4double WaterLeft_HalfX	= Water_thick*0.5;
	const G4double WaterLeft_HalfY	= dPolishLeadHalfY;
	const G4double WaterLeft_HalfZ	= dPolishLeadHalfZ;

	const G4double WaterRight_HalfX	= Water_thick*0.5;
	const G4double WaterRight_HalfY	= dPolishLeadHalfY;
	const G4double WaterRight_HalfZ	= dPolishLeadHalfZ;
	
	const G4double WaterTop_x = 0.0;
	const G4double WaterTop_y = 0.0;
	const G4double WaterTop_z = dShieldHalfZ + WaterTop_thick*0.5;

	const G4double WaterBack_x = 0.0;
	const G4double WaterBack_y = dShieldHalfY + Water_thick*0.5;
	const G4double WaterBack_z = 0.0;

	const G4double WaterLeft_x = -dShieldHalfX - Water_thick*0.5;
	const G4double WaterLeft_y = 0.0;
	const G4double WaterLeft_z = 0.0;

	const G4double WaterRight_x = dShieldHalfX + Water_thick*0.5;
	const G4double WaterRight_y = 0.0;
	const G4double WaterRight_z = 0.0;


	//====================== cryostat support bars =================================
	const G4double SBar_thick		= 50 *mm;
	const G4double SBar_height		= 60 *mm;
	const G4double SBarSide_length	= 730 *mm;
	const G4double SBarFront_width	= 610 *mm;
	
	const G4double dSBarSideLeftOffsetX		= -255 *mm - SBar_thick*0.5;
	const G4double dSBarSideRightOffsetX	= 255 *mm + SBar_thick*0.5;
	const G4double dSBarSideOffsetY			= -dCavityHalfY + SBarSide_length/2;
	const G4double dSBarFrontOffsetY		= -dCavityHalfY + SBarSide_length + SBar_thick/2;

	const G4double dSBarOffsetZ		= 362 - 800 - SBar_height/2; //top of the top flange - distance to the top of the bar
	
	//--------------------------------- polish lead ---------------------------------
	G4Box *pShieldPolishLeadBox = new G4Box("ShieldPolishLeadBox", dPolishLeadHalfX, dPolishLeadHalfY, dPolishLeadHalfZ);
	G4Material *Lead = G4Material::GetMaterial("Lead");
	m_pShieldPolishLeadLogicalVolume = new G4LogicalVolume(pShieldPolishLeadBox, Lead, "ShieldPolishLeadLogicalVolume", 0, 0, 0);
	m_pShieldPolishLeadPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPolishLeadOffsetZ), 
		m_pShieldPolishLeadLogicalVolume, "ShieldPolishLead", m_pLabLogicalVolume, false, 0);

	//--------------------------------- french lead ---------------------------------
	G4Box *pShieldFrenchLeadBox = new G4Box("ShieldFrenchLeadBox", dFrenchLeadHalfX, dFrenchLeadHalfY, dFrenchLeadHalfZ);
	m_pShieldFrenchLeadLogicalVolume = new G4LogicalVolume(pShieldFrenchLeadBox, Lead, "ShieldFrenchLeadLogicalVolume", 0, 0, 0);
	m_pShieldFrenchLeadPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dFrenchLeadOffsetZ), 
		m_pShieldFrenchLeadLogicalVolume, "ShieldFrenchLead", m_pShieldPolishLeadLogicalVolume, false, 0);
	
/*	//--------------------------------- inner lead ----------------------------------
	G4Box *pShieldInnerLeadBox = new G4Box("ShieldInnerLeadBox", dInnerLeadHalfX, dInnerLeadHalfY, dInnerLeadHalfZ);
	m_pShieldInnerLeadLogicalVolume = new G4LogicalVolume(pShieldInnerLeadBox, Lead, "ShieldInnerLeadLogicalVolume", 0, 0, 0);
	m_pShieldInnerLeadPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dInnerLeadOffsetZ), 
		m_pShieldInnerLeadLogicalVolume, "ShieldInnerLead", m_pShieldFrenchLeadLogicalVolume, false, 0);
*/	
	//-------------------------------- polyethylene ---------------------------------
	G4Box *pShieldPolyethyleneBox = new G4Box("ShieldPolyethyleneBox", dPolyethyleneHalfX, dPolyethyleneHalfY, dPolyethyleneHalfZ);
	G4Material *Polyethylene = G4Material::GetMaterial("Polyethylene");
	m_pShieldPolyethyleneLogicalVolume = new G4LogicalVolume(pShieldPolyethyleneBox, Polyethylene, "ShieldPolyethyleneLogicalVolume", 0, 0, 0);
	m_pShieldPolyethylenePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPolyethyleneOffsetZ), 
		m_pShieldPolyethyleneLogicalVolume, "ShieldPolyethylene", m_pShieldFrenchLeadLogicalVolume, false, 0);

	//---------------------------- polyethylene ground slab -------------------------------------
	G4Box* pPolyethyleneSlab = new G4Box("PolyethyleneSlabBox",dPolySlab_HalfX,dPolySlab_HalfY,dPolySlab_HalfZ);
	m_pPolyethyleneSlabLogicalVolume = new G4LogicalVolume(pPolyethyleneSlab, Polyethylene, "PolyethyleneSlabLogicalVolume", 0, 0, 0);
	m_pPolyethyleneSlabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dPolySlabOffsetX, dPolySlabOffsetY, dPolySlabOffsetZ),
		m_pPolyethyleneSlabLogicalVolume, "PolyethyleneSlabPhysicalVolume", m_pLabLogicalVolume, false, 0);
	
	//---------------------------------- copper -------------------------------------
	G4Box *pShieldCopperBox = new G4Box("ShieldCopperBox", dCopperHalfX, dCopperHalfY, dCopperHalfZ);	
	G4Material *Copper = G4Material::GetMaterial("Copper");
	m_pShieldCopperLogicalVolume = new G4LogicalVolume(pShieldCopperBox, Copper, "ShieldCopperLogicalVolume", 0, 0, 0);
	m_pShieldCopperPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dCopperOffsetZ), 
		m_pShieldCopperLogicalVolume, "ShieldCopper", m_pShieldPolyethyleneLogicalVolume, false, 0);

	//---------------------------------- water tanks --------------------------------
	G4Material *Water = G4Material::GetMaterial("G4_WATER");

	G4Box *pWaterTopBox 	= new G4Box("WaterTopBox", 	WaterTop_HalfX, WaterTop_HalfY, WaterTop_HalfZ);
	G4Box *pWaterBackBox 	= new G4Box("WaterBackBox", WaterBack_HalfX, WaterBack_HalfY, WaterBack_HalfZ);
	G4Box *pWaterLeftBox 	= new G4Box("WaterLeftBox", WaterLeft_HalfX, WaterLeft_HalfY, WaterLeft_HalfZ);
	G4Box *pWaterRightBox 	= new G4Box("WaterRightBox", WaterRight_HalfX, WaterRight_HalfY, WaterRight_HalfZ);
	
	m_pWaterTopLogicalVolume 	= new G4LogicalVolume(pWaterTopBox, 	Water, "WaterTopLogicalVolume", 0, 0, 0);
	m_pWaterBackLogicalVolume 	= new G4LogicalVolume(pWaterBackBox, 	Water, "WaterBackLogicalVolume", 0, 0, 0);
	m_pWaterLeftLogicalVolume 	= new G4LogicalVolume(pWaterLeftBox, 	Water, "WaterLeftLogicalVolume", 0, 0, 0);
	m_pWaterRightLogicalVolume 	= new G4LogicalVolume(pWaterRightBox, 	Water, "WaterRightLogicalVolume", 0, 0, 0);
	
	m_pWaterTopPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(WaterTop_x, WaterTop_y, WaterTop_z), m_pWaterTopLogicalVolume, "WaterTop", m_pLabLogicalVolume, false, 0);
	m_pWaterBackPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(WaterBack_x, WaterBack_y, WaterBack_z), m_pWaterBackLogicalVolume, "WaterBack", m_pLabLogicalVolume, false, 0);
	m_pWaterLeftPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(WaterLeft_x, WaterLeft_y, WaterLeft_z), m_pWaterLeftLogicalVolume, "WaterLeft", m_pLabLogicalVolume, false, 0);
	m_pWaterRightPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(WaterRight_x, WaterRight_y, WaterRight_z), m_pWaterRightLogicalVolume, "WaterRight", m_pLabLogicalVolume, false, 0);
	
	//----------------------------------- cavity ------------------------------------
	G4Box *pShieldCavityBox = new G4Box("ShieldCavityBox", dCavityHalfX, dCavityHalfY, dCavityHalfZ);
	G4Material *Air = G4Material::GetMaterial("G4_AIR");
	m_pShieldCavityLogicalVolume = new G4LogicalVolume(pShieldCavityBox, Air, "ShieldCavityLogicalVolume", 0, 0, 0);
	m_pShieldCavityPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dCavityOffsetZ),
		m_pShieldCavityLogicalVolume, "ShieldCavity", m_pShieldCopperLogicalVolume, false, 0);
	m_pMotherLogicalVolume = m_pShieldCavityLogicalVolume;
	
/*	//-------------------------------- polyethylene, bottom plate 5 cm thick
		G4Box *pShieldPolyethyleneBottomPlate = new G4Box("ShieldPolyethyleneBottomPlate", dPolyethyleneHalfX, dPolyethyleneHalfY, dPolyethyleneBottomPlateHalfZ);
	m_pShieldPolyethyleneBottomPlateLogicalVolume = new G4LogicalVolume(pShieldPolyethyleneBottomPlate, Polyethylene, "ShieldPolyethyleneBottomPlateLogicalVolume", 0, 0, 0);
	m_pShieldPolyethyleneBottomPlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPolyethyleneBottomPlateOffsetZ), 
		m_pShieldPolyethyleneBottomPlateLogicalVolume, "ShieldPolyethyleneBottomPlate", m_pShieldCavityLeadLogicalVolume, false, 0);
*/	
	//----------------------------------- cryostat support bars ---------------------
	G4Box* pSBarSideBox		= new G4Box("SBarSideBox",	SBar_thick*0.5, SBarSide_length*0.5, SBar_height*0.5);
	G4Box* pSBarFrontBox	= new G4Box("SBarFrontBox",	SBarFront_width*0.5, SBar_thick*0.5, SBar_height*0.5);
	
	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
	
	m_pSBarSideLeftLogicalVolume	= new G4LogicalVolume(pSBarSideBox, SS304LSteel, "SBarSideLeftLogicalVolume", 0, 0, 0);
	m_pSBarSideRightLogicalVolume	= new G4LogicalVolume(pSBarSideBox, SS304LSteel, "SBarSideRightLogicalVolume", 0, 0, 0);
	m_pSBarFrontLogicalVolume		= new G4LogicalVolume(pSBarFrontBox, SS304LSteel, "SBarFrontLogicalVolume", 0, 0, 0);
	
	m_pSBarSideLeftPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dSBarSideLeftOffsetX, dSBarSideOffsetY, dSBarOffsetZ), 
		m_pSBarSideLeftLogicalVolume, "SBarSideLeft", m_pShieldCavityLogicalVolume, false, 0);
		
	m_pSBarSideRightPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dSBarSideRightOffsetX, dSBarSideOffsetY, dSBarOffsetZ), 
		m_pSBarSideRightLogicalVolume, "SBarSideRight", m_pShieldCavityLogicalVolume, false, 0);
	
	m_pSBarFrontPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0.0, dSBarFrontOffsetY, dSBarOffsetZ), 
		m_pSBarFrontLogicalVolume, "SBarFront", m_pShieldCavityLogicalVolume, false, 0);

	//----------------------------------- lead brick for the source with the copper holder ----------
	G4double Lead_cutoff 	= 1.3 *cm;
	G4double Copper_cutoff	= 0.6 *cm;
	G4double Lead_offset 	= 23.51 *cm;

	G4double LeadBrick_width	= 10.0 *cm;
	G4double LeadBrick_length	= 16.0 *cm - Lead_cutoff;
	G4double LeadBrick_height	= 10.0 *cm;

	G4double LeadBrickTor_oD	= 1.5 *cm;
	G4double LeadBrickTor_iD	= 0.0 *cm;
	G4double LeadBrickTor_R		= 30.63 *cm;
	//G4double LeadBrickTor_R		= 30.63 *cm;

	G4double LongHolder_width 	= 1.91 *cm;
	G4double LongHolder_length 	= 5.08 *cm;
	G4double LongHolder_height 	= 25.67 *cm;
	
	G4double LeadBrickCollimatorHole_width 		= 0.5 *cm;
	G4double LeadBrickCollimatorHole_height		= 0.5 *cm;
	G4double LeadBrickCollimatorHole_angle		= 7 *deg;
	G4double LeadBrickCollimatorHole_length		= 15.14 *cm;
	G4double LeadBrickCollimatorHole_offsetY	= 2.89 *cm;

	// connecting piece
	G4double ShortHolder1_width		= 5.08 *cm;
	G4double ShortHolder1_length	= LongHolder_length;
	G4double ShortHolder1_height	= 1.91 *cm;
	// side piece
	G4double ShortHolder2_width		= 1.91 *cm;
	G4double ShortHolder2_length	= 5.08 *cm;
	G4double ShortHolder2_height	= 7.94 *cm;
	// upper piece
	G4double ShortHolder3_width		= 5.08 *cm;
	G4double ShortHolder3_length	= 16.0 *cm -Copper_cutoff;
	G4double ShortHolder3_height	= 1.91 *cm;

	G4RotationMatrix Rot0;
	Rot0.rotateX (0. *deg);
	G4ThreeVector LeadBrickTor_V 	= G4ThreeVector(LeadBrickTor_R, 0, 0);
	G4Transform3D LeadBrickTor_T 	= G4Transform3D(Rot0, LeadBrickTor_V);
	
	
	G4RotationMatrix Rot7;
	Rot7.rotateZ (353. *deg);
	G4ThreeVector LeadBrickCollimatorHole_V	= G4ThreeVector(LeadBrick_length*0.5, LeadBrickCollimatorHole_offsetY, 0);
	G4Transform3D LeadBrickCollimatorHole_T = G4Transform3D(Rot7, LeadBrickCollimatorHole_V);
	

	G4Box* pLeadBrickBox				= new G4Box("LeadBrickBox",	LeadBrick_length*0.5, LeadBrick_width*0.5, LeadBrick_height*0.5);
	G4Torus* pLeadBrickTor				= new G4Torus("LeadBrickTor", 0.5*LeadBrickTor_iD, 0.5*LeadBrickTor_oD, LeadBrickTor_R, 0, 330. *deg);
	G4Box* pLeadBrickCollimatorHoleBox 	= new G4Box("LeadBrickCollimatorHoleBox", LeadBrickCollimatorHole_length*0.5, LeadBrickCollimatorHole_width*0.5, LeadBrickCollimatorHole_height*0.5);

	G4SubtractionSolid* pLeadBrickBox1	= new G4SubtractionSolid("LeadBrickBox1",	pLeadBrickBox,	pLeadBrickTor, LeadBrickTor_T);
	G4SubtractionSolid* pLeadBrickBox2	= new G4SubtractionSolid("LeadBrickBox2",	pLeadBrickBox1,	pLeadBrickCollimatorHoleBox, LeadBrickCollimatorHole_T);

	G4Box* pLongHolderBox	= new G4Box("LongHolderBox",LongHolder_width*0.5, LongHolder_length*0.5, LongHolder_height*0.5);
	G4Box* pShortHolder1Box	= new G4Box("ShortHolder1Box",ShortHolder1_width*0.5, ShortHolder1_length*0.5, ShortHolder1_height*0.5);
	G4Box* pShortHolder2Box	= new G4Box("ShortHolder2Box",ShortHolder2_width*0.5, ShortHolder2_length*0.5, ShortHolder2_height*0.5);
	G4Box* pShortHolder3Box	= new G4Box("ShortHolder3Box",ShortHolder3_length*0.5, ShortHolder3_width*0.5, ShortHolder3_height*0.5);

	m_pLeadBrickLogicalVolume		= new G4LogicalVolume(pLeadBrickBox2, 	Lead, 	"LeadBrickLogicalVolume", 0, 0, 0);
	m_pLongHolderLogicalVolume		= new G4LogicalVolume(pLongHolderBox, 	Copper, "LongHolderLogicalVolume", 0, 0, 0);
	m_pShortHolder1LogicalVolume	= new G4LogicalVolume(pShortHolder1Box, Copper, "ShortHolder1LogicalVolume", 0, 0, 0);
	m_pShortHolder2LogicalVolume	= new G4LogicalVolume(pShortHolder2Box, Copper, "ShortHolder2LogicalVolume", 0, 0, 0);
	m_pShortHolder3LogicalVolume	= new G4LogicalVolume(pShortHolder3Box, Copper, "ShortHolder2LogicalVolume", 0, 0, 0);

	G4double dLongHolderOffsetX	= dSBarSideLeftOffsetX - SBar_thick*0.5 - LongHolder_width*0.5;
	G4double dLongHolderOffsetY	= 0;
	G4double dLongHolderOffsetZ	= dSBarOffsetZ - SBar_height*0.5 + LongHolder_height*0.5;

	m_pLongHolderPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dLongHolderOffsetX,dLongHolderOffsetY,dLongHolderOffsetZ), 
		m_pLongHolderLogicalVolume, "LongHolder", m_pShieldCavityLogicalVolume, false, 0);

	G4double dShortHolder1OffsetX	= dLongHolderOffsetX + LongHolder_width*0.5 + ShortHolder1_width*0.5;
	G4double dShortHolder1OffsetY	= 0;
	G4double dShortHolder1OffsetZ	= dLongHolderOffsetZ - LongHolder_height*0.5 + ShortHolder2_height - ShortHolder1_height*0.5;

	m_pShortHolder1PhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dShortHolder1OffsetX,dShortHolder1OffsetY,dShortHolder1OffsetZ), 
		m_pShortHolder1LogicalVolume, "ShortHolder1", m_pShieldCavityLogicalVolume, false, 0);

	G4double dShortHolder2OffsetX	= dShortHolder1OffsetX + ShortHolder1_width*0.5 + ShortHolder2_width*0.5;
	G4double dShortHolder2OffsetY	= 0;
	G4double dShortHolder2OffsetZ	= dLongHolderOffsetZ - LongHolder_height*0.5 + ShortHolder2_height*0.5;

	m_pShortHolder2PhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dShortHolder2OffsetX,dShortHolder2OffsetY,dShortHolder2OffsetZ), 
		m_pShortHolder2LogicalVolume, "ShortHolder2", m_pShieldCavityLogicalVolume, false, 0);

	G4double dShortHolder3OffsetX	= dLongHolderOffsetX + Copper_cutoff;
	G4double dShortHolder3OffsetY	= 0;
	G4double dShortHolder3OffsetZ	= dLongHolderOffsetZ + LongHolder_height*0.5 + ShortHolder3_height*0.5;

	m_pShortHolder3PhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dShortHolder3OffsetX,dShortHolder3OffsetY,dShortHolder3OffsetZ), 
		m_pShortHolder3LogicalVolume, "ShortHolder3", m_pShieldCavityLogicalVolume, false, 0);
	
	G4double dLeadBrickOffsetX	= dShortHolder3OffsetX + ShortHolder3_length*0.5 - LeadBrick_length*0.5;
	G4double dLeadBrickOffsetY	= 0;
	G4double dLeadBrickOffsetZ	= dShortHolder3OffsetZ + ShortHolder3_height*0.5 + LeadBrick_height*0.5;

	m_pLeadBrickPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dLeadBrickOffsetX,dLeadBrickOffsetY,dLeadBrickOffsetZ), 
		m_pLeadBrickLogicalVolume, "LeadBrick", m_pShieldCavityLogicalVolume, false, 0);
	
	//----------------------------------- copper pipe to insert the source ----------
	G4double SPipe_oD			= 10.0 *mm;
	G4double SPipe_wall			= 1.0 *mm;
	G4double SPipe_iD			= SPipe_oD - 2*SPipe_wall;
	G4double SPipeLine_length	= 400. *mm;
	G4double SPipe_R			= 30.6 *cm; // 30.63 originally
	
	G4RotationMatrix RotM90x;
	G4RotationMatrix RotP90x;

	RotP90x.rotateX (90. *deg);
	RotM90x.rotateX(-90. *deg);
	
	G4Tubs* SPipeLineR	= new G4Tubs("SPipeLineR", 0.5*SPipe_iD, 0.5*SPipe_oD, 0.5*SPipeLine_length, 0. * deg, 360. * deg);
	G4Torus* SPipeTor	= new G4Torus("SPipeTor", 0.5*SPipe_iD, 0.5*SPipe_oD, SPipe_R, 0, 330. *deg);
	
	G4ThreeVector SPipeLineR_V = G4ThreeVector(SPipe_R,		-SPipeLine_length*0.5, 0);
	
	G4Transform3D SPipeLineR_T = G4Transform3D(RotP90x, SPipeLineR_V);
	
	G4UnionSolid* SourcePipe2	= new G4UnionSolid("SourcePipe2",	SPipeTor,	SPipeLineR, SPipeLineR_T);
	
	G4double dSourcePipeOffsetZ	= dLeadBrickOffsetZ;
	G4cout <<"Z coordinate of the copper pipe = "<< dLeadBrickOffsetZ << G4endl;
	
	m_pSourcePipeLogicalVolume	= new G4LogicalVolume(SourcePipe2,	Copper, "SourcePipeLogicalVolume");
	m_pSourcePipePhysicalVolume	= new G4PVPlacement(0,G4ThreeVector(0,0,dSourcePipeOffsetZ), 
		m_pSourcePipeLogicalVolume, "SourcePipe", m_pShieldCavityLogicalVolume, false, 0);



/*	//----------------------------------- wire (Th228 source) ---------------------
	G4double Wire_iD		= 0.0 *mm;
	G4double Wire_oD		= 1.8 *mm;
	G4double Wire_length	= 400. *mm;
	G4double Wire_R			= 30.6 *cm; // 30.63 originally
		
	G4Torus* WireTor		= new G4Torus("WireTor", 0.5*Wire_iD, 0.5*Wire_oD, Wire_R, 0, 360. *deg);
	
	G4double dWireOffsetZ	= dLeadBrickOffsetZ;
	
	m_pWireLogicalVolume	= new G4LogicalVolume(WireTor,	SS304LSteel, "WireLogicalVolume");
	m_pWirePhysicalVolume	= new G4PVPlacement(0,G4ThreeVector(0,0,dWireOffsetZ), 
		m_pWireLogicalVolume, "Wire", m_pShieldCavityLogicalVolume, false, 0);
*/				

/*	//----------------------------------- lead brick for the source (simplified)---------------------
	G4double LeadBrick_height	= 100. *mm;
	G4double LeadBrick_width	= 100. *mm;
	G4double LeadBrick_thick	= 50. *mm;

	G4Box* pLeadBrickBox		= new G4Box("LeadBrickBox",	LeadBrick_thick*0.5, LeadBrick_width*0.5, LeadBrick_height*0.5);
	
	m_pLeadBrickLogicalVolume	= new G4LogicalVolume(pLeadBrickBox, Lead, "LeadBrickLogicalVolume", 0, 0, 0);
	
	m_pLeadBrickPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(250,0,-163), 
	m_pLeadBrickLogicalVolume, "LeadBrick", m_pShieldCavityLogicalVolume, false, 0);
*/		

/*	//----------------------------------- small lead brick collimator for veto mapping -------------
	G4double Collimator_height	= 32.0 *mm;
	G4double Collimator_width	= 28.0 *mm;
	G4double Collimator_depth	= 50.0 *mm;
	G4double CollimatorHole_oD	= 3.7 *mm;

	G4Box* 	pCollimatorBox		= new G4Box("CollimatorBox",	Collimator_width*0.5, Collimator_depth*0.5, Collimator_height*0.5);
	G4Tubs* pCollimatorHole		= new G4Tubs("CollimatorHole",	0.0, CollimatorHole_oD*0.5, Collimator_width*0.5, 0. *deg, 360. *deg);

	G4RotationMatrix RotM90y;
	G4RotationMatrix RotP90y;

	RotP90y.rotateY (90. *deg);
	RotM90y.rotateY(-90. *deg);

	G4ThreeVector CollimatorHole_V = G4ThreeVector(0, 0, 0);	
	G4Transform3D CollimatorHole_T = G4Transform3D(RotP90y, CollimatorHole_V);

	G4SubtractionSolid *pCollimatorSolid = new G4SubtractionSolid("CollimatorSolid",pCollimatorBox, pCollimatorHole, CollimatorHole_T);

	m_pCollimatorLogicalVolume	= new G4LogicalVolume(pCollimatorSolid, Lead, "CollimatorLogicalVolume", 0, 0, 0);
	
	G4double CollimatorOffsetX 		= 231.0 + Collimator_width*0.5;
	G4double CollimatorPosition1	= 77.0;
	G4double CollimatorPosition2	= -0.8;
	G4double CollimatorPosition3	= -103.0;
	G4double CollimatorPosition4	= -198.0;
	G4double CollimatorPosition5	= -288.0;
	G4double CollimatorPosition6	= -363.0;
	
	G4double CollimatorOffsetZ		= CollimatorPosition2;
	
	m_pCollimatorPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(CollimatorOffsetX,0,CollimatorOffsetZ), 
	m_pCollimatorLogicalVolume, "Collimator", m_pShieldCavityLogicalVolume, false, 0);
*/		


/*	//----------------------------------- small capsule with the source (POSITION: lead brick collimator) ----------
	G4double SourceOrb_oD	= 5. *mm;

	G4Orb* pSourceOrb	= new G4Orb("SourceOrb", SourceOrb_oD*0.5);
	
	m_pCalibrationSourceLogicalVolume	= new G4LogicalVolume(pSourceOrb, Air, "SourceLogicalVolume", 0, 0, 0);
	
	G4double dSourceOffsetX = -30.3 *cm;
	G4double dSourceOffsetY = 3.7 *cm;
	G4double dSourceOffsetZ = dLeadBrickOffsetZ;
	
	m_pCalibrationSourcePhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dSourceOffsetX,dSourceOffsetY,dSourceOffsetZ), 
	m_pCalibrationSourceLogicalVolume, "CalibrationSource", m_pShieldCavityLogicalVolume, false, 0);
*/
	//----------------------------------- small capsule with the source (POSITION: lead brick Center-5deg) ----------
	G4double SourceOrb_oD	= 3. *mm;

	G4Orb* pSourceOrb	= new G4Orb("SourceOrb", SourceOrb_oD*0.5);
	
	m_pCalibrationSourceLogicalVolume	= new G4LogicalVolume(pSourceOrb, Air, "SourceLogicalVolume", 0, 0, 0);
	
	G4double dSourceOffsetX = -30.19 *cm;
	G4double dSourceOffsetY = -2.64 *cm;
	G4double dSourceOffsetZ = dLeadBrickOffsetZ;
	
	m_pCalibrationSourcePhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dSourceOffsetX,dSourceOffsetY,dSourceOffsetZ), 
	m_pCalibrationSourceLogicalVolume, "CalibrationSource", m_pShieldCavityLogicalVolume, false, 0);


/*	//----------------------------------- small capsule with the source (POSITION: RED) ------------
	G4double SourceOrb_oD	= 5. *mm;

	G4Orb* pSourceOrb	= new G4Orb("SourceOrb", SourceOrb_oD*0.5);
	
	m_pCalibrationSourceLogicalVolume	= new G4LogicalVolume(pSourceOrb, Air, "SourceLogicalVolume", 0, 0, 0);
	
	G4double dSourceOffsetX = 28.6 *cm;
	G4double dSourceOffsetY = 10.1 *cm;
	G4double dSourceOffsetZ = dLeadBrickOffsetZ;
	
	m_pCalibrationSourcePhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dSourceOffsetX,dSourceOffsetY,dSourceOffsetZ), 
	m_pCalibrationSourceLogicalVolume, "CalibrationSource", m_pShieldCavityLogicalVolume, false, 0);
*/


/*	//----------------------------------- small capsule with the source (for VETO MAPPING) ------------
	G4double SourceOrb_oD	= 1.0 *mm;

	G4Orb* pSourceOrb	= new G4Orb("SourceOrb", SourceOrb_oD*0.5);
	
	m_pCalibrationSourceLogicalVolume	= new G4LogicalVolume(pSourceOrb, Air, "SourceLogicalVolume", 0, 0, 0);
	
	G4double dSourceOffsetX = CollimatorOffsetX + SourceOrb_oD/2 + 2.0*mm;
	G4double dSourceOffsetY = 0.0 *mm;
	G4double dSourceOffsetZ = CollimatorOffsetZ;
	
	m_pCalibrationSourcePhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dSourceOffsetX,dSourceOffsetY,dSourceOffsetZ), 
	m_pCalibrationSourceLogicalVolume, "CalibrationSource", m_pShieldCavityLogicalVolume, false, 0);
*/

	//---------------------------------- attributes ---------------------------------	
	//G4Colour hPolishLeadColor(0.086, 0.498, 0.788, 0.1);
	//G4Colour hFrenchLeadColor(0.561, 0.078, 0.580, 0.1);
	//G4Colour hInnerLeadColor(0.561, 0.078, 0.580, 0.1);
	//G4Colour hPolyethyleneColor(1.000, 1.000, 1.000, 0.1);
	//G4Colour hCopperColor(1.0, 1.0, 0.0);
	//G4Colour hCavityColor(0.000, 0.000, 0.000, 0.0);
	//G4Colour hLeadColour(0.5, 0.5, 0.5);		// grey
	G4Colour hLeadColour(0.0, 0.0, 0.0);		// grey
	G4Colour hPolyColour(1.0,	1.0, 1.0);		// white
	//G4Colour hCopperColour(1.0, 1.0, 0.0);		// yellow
	G4Colour hCopperColour(0.0, 0.0, 0.0);		// yellow
	//G4Colour hWaterColour(0.0, 0.0, 1.0);		// xlblue
	G4Colour hWaterColour(0.0, 0.0, 0.0);		// xlblue

	//G4Colour hSBarsColour(0.0, 0.0, 1.0);		// xlblue
	//G4Colour hSBarsColour(51, 51, 51);		// dark grey
	G4Colour hSBarsColour(0, 0, 0);			// black

	G4Colour hSourceColour(0.5, 0.0, 0.0);		// xlred

	
	G4VisAttributes *pShieldPolishLeadVisAtt = new G4VisAttributes(hLeadColour);
	pShieldPolishLeadVisAtt->SetVisibility(true);
	m_pShieldPolishLeadLogicalVolume->SetVisAttributes(pShieldPolishLeadVisAtt);//invisible

	G4VisAttributes *pShieldFrenchLeadVisAtt = new G4VisAttributes(hLeadColour);
	pShieldFrenchLeadVisAtt->SetVisibility(true);
	m_pShieldFrenchLeadLogicalVolume->SetVisAttributes(pShieldFrenchLeadVisAtt);//invisible
	
	//G4VisAttributes *pShieldInnerLeadVisAtt = new G4VisAttributes(hLeadColour);
	//pShieldInnerLeadVisAtt->SetVisibility(true);
	//m_pShieldInnerLeadLogicalVolume->SetVisAttributes(pShieldInnerLeadVisAtt);//invisible

	G4VisAttributes *pShieldPolyethyleneVisAtt = new G4VisAttributes(hPolyColour);
	pShieldPolyethyleneVisAtt->SetVisibility(true);
	m_pShieldPolyethyleneLogicalVolume->SetVisAttributes(pShieldPolyethyleneVisAtt);//invisible
	
	G4VisAttributes *pShieldCopperVisAtt = new G4VisAttributes(hCopperColour);
	pShieldCopperVisAtt->SetVisibility(true);
	m_pShieldCopperLogicalVolume->SetVisAttributes(pShieldCopperVisAtt);//invisible

	G4VisAttributes *pShieldWaterVisAtt = new G4VisAttributes(hWaterColour);
	pShieldWaterVisAtt->SetVisibility(true);
	m_pWaterTopLogicalVolume->SetVisAttributes(pShieldWaterVisAtt);//invisible
	m_pWaterBackLogicalVolume->SetVisAttributes(pShieldWaterVisAtt);//invisible
	m_pWaterLeftLogicalVolume->SetVisAttributes(pShieldWaterVisAtt);//invisible
	m_pWaterRightLogicalVolume->SetVisAttributes(pShieldWaterVisAtt);//invisible
	
	G4VisAttributes *pShieldCavityVisAtt = new G4VisAttributes(hPolyColour);
	pShieldCavityVisAtt->SetVisibility(true);
	m_pShieldCavityLogicalVolume->SetVisAttributes(pShieldCavityVisAtt); //invisible

	G4VisAttributes *pSBarsVisAtt = new G4VisAttributes(hSBarsColour);
	pSBarsVisAtt->SetVisibility(true);
	m_pSBarSideLeftLogicalVolume->SetVisAttributes(pSBarsVisAtt);
	m_pSBarSideRightLogicalVolume->SetVisAttributes(pSBarsVisAtt);
	m_pSBarFrontLogicalVolume->SetVisAttributes(pSBarsVisAtt);

	G4VisAttributes *pSourcePipeVisAtt = new G4VisAttributes(hCopperColour);
	pSourcePipeVisAtt->SetVisibility(true);
	//m_pLongHolderLogicalVolume->SetVisAttributes(pSourcePipeVisAtt);
	//m_pShortHolder1LogicalVolume->SetVisAttributes(pSourcePipeVisAtt);
	//m_pShortHolder2LogicalVolume->SetVisAttributes(pSourcePipeVisAtt);
	//m_pShortHolder3LogicalVolume->SetVisAttributes(pSourcePipeVisAtt);
	m_pSourcePipeLogicalVolume->SetVisAttributes(pSourcePipeVisAtt);

	G4VisAttributes *pLeadBrickVisAtt = new G4VisAttributes(hLeadColour);
	pLeadBrickVisAtt->SetVisibility(true);
	//m_pLeadBrickLogicalVolume->SetVisAttributes(pLeadBrickVisAtt);
	//m_pCollimatorLogicalVolume->SetVisAttributes(pLeadBrickVisAtt);

	G4VisAttributes *pSourceVisAtt = new G4VisAttributes(hSourceColour);
	pSourceVisAtt->SetVisibility(true);
	//m_pCalibrationSourceLogicalVolume->SetVisAttributes(pSourceVisAtt);
	//m_pWireLogicalVolume->SetVisAttributes(pSourceVisAtt);

	G4VisAttributes *pInvisibleVisAtt = new G4VisAttributes();
	pInvisibleVisAtt->SetVisibility(false);
	//m_pShieldPolishLeadLogicalVolume	->SetVisAttributes(pInvisibleVisAtt);
	//m_pShieldFrenchLeadLogicalVolume	->SetVisAttributes(pInvisibleVisAtt);
	//m_pShieldPolyethyleneLogicalVolume->SetVisAttributes(pInvisibleVisAtt);
	//m_pShieldCopperLogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	//m_pShieldCavityLogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	//m_pWaterTopLogicalVolume			->SetVisAttributes(pInvisibleVisAtt);
	//m_pWaterBackLogicalVolume			->SetVisAttributes(pInvisibleVisAtt);
	//m_pWaterLeftLogicalVolume			->SetVisAttributes(pInvisibleVisAtt);
	//m_pWaterRightLogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	m_pLongHolderLogicalVolume			->SetVisAttributes(pInvisibleVisAtt);
	m_pShortHolder1LogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	m_pShortHolder2LogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	m_pShortHolder3LogicalVolume		->SetVisAttributes(pInvisibleVisAtt);
	m_pLeadBrickLogicalVolume			->SetVisAttributes(pInvisibleVisAtt);
	m_pCalibrationSourceLogicalVolume	->SetVisAttributes(pInvisibleVisAtt);

}

void
Xenon100DetectorConstruction::ConstructXenon()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< xenon >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dVesselInnerWallHeight = GetGeometryParameter("CryostatVesselInnerWallHeight");
	const G4double dVesselInnerWallRadius = GetGeometryParameter("CryostatVesselInnerWallRadius");
	const G4double dVesselInnerWallDomeHeight = GetGeometryParameter("CryostatVesselInnerWallDomeHeight");
	const G4double dVesselInnerWallDomeRadius = GetGeometryParameter("CryostatVesselInnerWallDomeRadius");

	const G4double dVesselFlangeToZero = GetGeometryParameter("CryostatVesselFlangeToZero");
	const G4double dVesselFlangeThickness = GetGeometryParameter("CryostatVesselFlangeThickness");

	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInnerWallDomeRadius = GetGeometryParameter("CryostatLidInnerWallDomeRadius");

	const G4double dLidInnerWallDomeTopToVetoLiquidLevel = GetGeometryParameter("CryostatLidInnerWallDomeTopToVetoLiquidLevel");
	const G4double dLidInnerWallDomeTopToLiquidLevel = GetGeometryParameter("CryostatLidInnerWallDomeTopToLiquidLevel");

	const G4double dLidInTubeRadius = GetGeometryParameter("CryostatLidInTubeRadius");
	const G4double dLidInTubeThickness = GetGeometryParameter("CryostatLidInTubeThickness");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");
	const G4double dLidInTubeOffsetX = GetGeometryParameter("CryostatLidInTubeOffsetX");

	const G4double dBellRadius = GetGeometryParameter("BellRadius");
	const G4double dBellTopThickness = GetGeometryParameter("BellTopThickness");
	const G4double dBellSideThickness = GetGeometryParameter("BellSideThickness");

	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
	const G4double dGXeHeight = GetGeometryParameter("GXeHeight");

	G4Material *LXe = G4Material::GetMaterial("LXe");
	G4Material *GXe = G4Material::GetMaterial("GXe");

	//================================ liquid xenon =================================

	//---------------------- tube with upper and lower bulges -----------------------
	const G4double dLXeRadius = dVesselInnerWallRadius;
	const G4double dLXeHalfZ = 0.5*dLXeHeight;

	const G4double dLXeUpperBulgeRadius = dLidInnerWallDomeRadius;
	const G4double dLXeUpperBulgeHeight = dLidInnerWallDomeHeight;
	const G4double dLXeUpperBulgeCutOffsetZ = dLXeHalfZ-dLXeUpperBulgeRadius;

	const G4double dLXeLowerBulgeRadius = dVesselInnerWallDomeRadius;
	const G4double dLXeLowerBulgeHeight = dVesselInnerWallDomeHeight;
	const G4double dLXeLowerBulgeCutOffsetZ = -dLXeHalfZ+dLXeLowerBulgeRadius;

	const G4double dLXeOffsetZ = dVesselFlangeToZero-dVesselFlangeThickness-dVesselInnerWallHeight-dVesselInnerWallDomeHeight+dLXeHalfZ;

	G4Tubs *pLXeTubs = new G4Tubs("LXeTubs", 0.*cm, dLXeRadius, dLXeHalfZ, 0.*deg, 360.*deg);

	G4Sphere *pLXeUpperBulgeCut = new G4Sphere("LXeUpperBulgeCut", dLXeUpperBulgeRadius,
		dLXeUpperBulgeRadius+2*dLXeUpperBulgeHeight, 0.*deg, 360.*deg, 0.*deg, 90.*deg);

	G4SubtractionSolid *pLXeTubeWithCut1 = new G4SubtractionSolid("LXeTubeWithCut1",
		pLXeTubs, pLXeUpperBulgeCut, 0, G4ThreeVector(0., 0., dLXeUpperBulgeCutOffsetZ));

	G4Sphere *pLXeLowerBulgeCut = new G4Sphere("LXeUpperBulgeCut", dLXeLowerBulgeRadius,
		dLXeLowerBulgeRadius+2*dLXeLowerBulgeHeight, 0.*deg, 360.*deg, 90.*deg, 180.*deg);

	G4SubtractionSolid *pLXeTubeWithCut2 = new G4SubtractionSolid("LXeTubeWithCut1",
		pLXeTubeWithCut1, pLXeLowerBulgeCut, 0, G4ThreeVector(0., 0., dLXeLowerBulgeCutOffsetZ));

	m_pLXeLogicalVolume = new G4LogicalVolume(pLXeTubeWithCut2, LXe, "LXeVolume", 0, 0, 0);

	m_pLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
		m_pLXeLogicalVolume, "LXe", m_pMotherLogicalVolume, false, 0);
	
	//================================ gaseous xenon ================================

	//----------------------------- veto gaseous xenon ------------------------------
	const G4double dVetoGXeRadius = dVesselInnerWallRadius;
	const G4double dVetoGXeHalfZ = 0.5*dLidInnerWallDomeTopToVetoLiquidLevel;
	const G4double dVetoGXeOffsetZ = dLXeHalfZ-dVetoGXeHalfZ;

	const G4double dVetoGXeUpperBulgeRadius = dLXeUpperBulgeRadius;
	const G4double dVetoGXeUpperBulgeHeight = dLXeUpperBulgeHeight;
	const G4double dVetoGXeUpperBulgeCutOffsetZ = dVetoGXeHalfZ-dVetoGXeUpperBulgeRadius;

	G4Tubs *pVetoGXeTubs = new G4Tubs("VetoGXeTube", 0.*cm, dVetoGXeRadius, dVetoGXeHalfZ, 0.*deg, 360.*deg);

	G4Sphere *pVetoGXeUpperBulgeCut = new G4Sphere("VetoGXeUpperBulgeCut", dVetoGXeUpperBulgeRadius,
		dVetoGXeUpperBulgeRadius+2*dVetoGXeUpperBulgeHeight, 0.*deg, 360.*deg, 0.*deg, 90.*deg);

	G4SubtractionSolid *pVetoGXeWithCut1 = new G4SubtractionSolid("VetoGXeWithCut1",
		pVetoGXeTubs, pVetoGXeUpperBulgeCut, 0, G4ThreeVector(0., 0., dVetoGXeUpperBulgeCutOffsetZ));

	m_pVetoGXeLogicalVolume = new G4LogicalVolume(pVetoGXeWithCut1, GXe, "VetoGXeLogicalVolume", 0, 0, 0);

	m_pVetoGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dVetoGXeOffsetZ),
		m_pVetoGXeLogicalVolume, "VetoGXe", m_pLXeLogicalVolume, false, 0);

	//--------------------------------- inside bell ---------------------------------
	const G4double dGXeRadius = dBellRadius-dBellSideThickness;
	const G4double dGXeHalfZ = 0.5*dGXeHeight;
	const G4double dGXeOffsetZ = dLXeHalfZ-dLidInnerWallDomeTopToLiquidLevel+dGXeHalfZ;

	G4Tubs *pGXeTubs = new G4Tubs("GXeTube", 0.*cm, dGXeRadius, dGXeHalfZ, 0.*deg, 360.*deg);

	m_pGXeLogicalVolume = new G4LogicalVolume(pGXeTubs, GXe, "GXeLogicalVolume", 0, 0, 0);

	m_pGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dGXeOffsetZ),
		m_pGXeLogicalVolume, "GXe", m_pLXeLogicalVolume, false, 0);

	//------------------------------ inside lid in tube -----------------------------
	const G4double dLidInTubeGXeRadius = dLidInTubeRadius-dLidInTubeThickness;
	const G4double dLidInTubeGXeHalfZ = 0.5*(dLidInTubeHeightInside-(dLidInnerWallDomeTopToVetoLiquidLevel-(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))+dBellTopThickness);
	const G4double dLidInTubeGXeOffsetZ = dLXeHalfZ-dLidInnerWallDomeTopToLiquidLevel+dGXeHeight+dLidInTubeGXeHalfZ;

	G4Tubs *pLidInTubeGXeTubs = new G4Tubs("LidInTubeGXeTube", 0.*cm, dLidInTubeGXeRadius,
		dLidInTubeGXeHalfZ, 0.*deg, 360.*deg);

	m_pLidInTubeGXeLogicalVolume = new G4LogicalVolume(pLidInTubeGXeTubs, GXe, "LidInTubeGXeLogicalVolume", 0, 0, 0);

	m_pLidInTubeGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeGXeOffsetZ),
		m_pLidInTubeGXeLogicalVolume, "LidInTubeGXe", m_pLXeLogicalVolume, false, 0);

	//------------------------------ xenon sensitivity ------------------------------
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	Xenon100LXeSensitiveDetector *pLXeSD = new Xenon100LXeSensitiveDetector("Xenon100/LXeSD");
	pSDManager->AddNewDetector(pLXeSD);
	m_pLXeLogicalVolume->SetSensitiveDetector(pLXeSD);

	//================================== attributes =================================
	//G4Colour hLXeColour(0.5, 0.0, 0.0);			// xlred
	//G4Colour hGXeColour(0.0, 0.75, 0.0);		// lgreen

	G4Colour hLXeColour(0.0, 0.0, 0.0);	
	G4Colour hGXeColour(0.0, 0.0, 0.0);	
	
	G4VisAttributes *pLXeVisAtt = new G4VisAttributes(hLXeColour);
	//pLXeVisAtt->SetVisibility(true);
	pLXeVisAtt->SetVisibility(false);
	m_pLXeLogicalVolume->SetVisAttributes(pLXeVisAtt);

	G4VisAttributes *pVetoGXeVisAtt = new G4VisAttributes(hGXeColour);
	//pVetoGXeVisAtt->SetVisibility(true);
	pVetoGXeVisAtt->SetVisibility(false);
	m_pVetoGXeLogicalVolume->SetVisAttributes(pVetoGXeVisAtt);

	G4VisAttributes *pGXeVisAtt = new G4VisAttributes(hGXeColour);
	//pGXeVisAtt->SetVisibility(true);
	pGXeVisAtt->SetVisibility(false);
	m_pGXeLogicalVolume->SetVisAttributes(pGXeVisAtt);
	m_pLidInTubeGXeLogicalVolume->SetVisAttributes(pGXeVisAtt);
}

void
Xenon100DetectorConstruction::ConstructBell()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bell >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
	const G4double dGXeHeight = GetGeometryParameter("GXeHeight");

	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");

	const G4double dLidInTubeRadius = GetGeometryParameter("CryostatLidInTubeRadius");
	const G4double dLidInTubeThickness = GetGeometryParameter("CryostatLidInTubeThickness");
	const G4double dLidInTubeOffsetX = GetGeometryParameter("CryostatLidInTubeOffsetX");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");
	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");

	const G4double dBellHeight = GetGeometryParameter("BellHeight");
	const G4double dBellRadius = GetGeometryParameter("BellRadius");
	const G4double dBellTopThickness = GetGeometryParameter("BellTopThickness");
	const G4double dBellSideThickness = GetGeometryParameter("BellSideThickness");

	const G4double dTopPmtTeflonHolderRadius = GetGeometryParameter("TopPmtTeflonHolderRadius");
	const G4double dTopPmtTeflonHolderThickness = GetGeometryParameter("TopPmtTeflonHolderThickness");
	const G4double dTopPmtTeflonHolderCutDepth = GetGeometryParameter("TopPmtTeflonHolderCutDepth");
	const G4double dTopPmtTeflonHolderCutWidth = GetGeometryParameter("TopPmtTeflonHolderCutWidth");
	const G4double dTopPmtTeflonHolderBottomCutWidth = GetGeometryParameter("TopPmtTeflonHolderBottomCutWidth");
	const G4double dTopPmtTeflonHolderBottomCutAngle = GetGeometryParameter("TopPmtTeflonHolderBottomCutAngle");
	const G4double dTopPmtTeflonHolderToBell = GetGeometryParameter("TopPmtTeflonHolderToBell");

	const G4double dBellPmtSupportRingThickness = GetGeometryParameter("BellPmtSupportRingThickness");
	const G4double dBellPmtSupportRingRadius = GetGeometryParameter("BellPmtSupportRingRadius");
	const G4double dBellPmtSupportRingHeight = GetGeometryParameter("BellPmtSupportRingHeight");

	const G4double dBellSupportCylinderRadius = GetGeometryParameter("BellSupportCylinderRadius");
	const G4double dBellSupportCylinderHeight = GetGeometryParameter("BellSupportCylinderHeight");

	const G4double dTopVetoAngleThickness = GetGeometryParameter("TopVetoAngleThickness");
	const G4double dTopVetoAngleLowerHeight = GetGeometryParameter("TopVetoAngleLowerHeight");
	const G4double dTopVetoAngleUpperHeight = GetGeometryParameter("TopVetoAngleUpperHeight");
	const G4double dTopVetoAngleLength = GetGeometryParameter("TopVetoAngleLength");
	const G4double dTopVetoAngleWidth = GetGeometryParameter("TopVetoAngleWidth");

	const G4double dUpperSideVetoAngleThickness = GetGeometryParameter("UpperSideVetoAngleThickness");
	const G4double dUpperSideVetoAngleHeight = GetGeometryParameter("UpperSideVetoAngleHeight");
	const G4double dUpperSideVetoAngleLength = GetGeometryParameter("UpperSideVetoAngleLength");
	const G4double dUpperSideVetoAngleWidth = GetGeometryParameter("UpperSideVetoAngleWidth");

	const G4double dLidInnerWallDomeTopToVetoLiquidLevel = GetGeometryParameter("CryostatLidInnerWallDomeTopToVetoLiquidLevel");
	const G4double dLidInnerWallDomeRadius = GetGeometryParameter("CryostatLidInnerWallDomeRadius");

	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
	G4Material *Copper = G4Material::GetMaterial("Copper");	
	G4Material *Teflon = G4Material::GetMaterial("Teflon");

	//===================================== bell ====================================
	const G4double dLXeHalfZ = 0.5*dLXeHeight;

	const G4double dBellHalfZ = 0.5*dBellHeight;
	const G4double dBellCutRadius = dBellRadius-dBellSideThickness;
	const G4double dBellCutHalfZ = dBellHeight-dBellTopThickness;
	const G4double dBellCutOffsetZ = -dBellHalfZ;
	const G4double dBellOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHalfZ;
	const G4double dBellInTubeCutRadius = dLidInTubeRadius-dLidInTubeThickness;
	
	G4Tubs *pBellTubs = new G4Tubs("BellTube", 0.*cm, dBellRadius, dBellHalfZ, 0.*deg, 360.*deg);

	G4Tubs *pBellCutTubs = new G4Tubs("BellCutTube", 0.*cm, dBellCutRadius, dBellCutHalfZ, 0.*deg, 360.*deg);

	G4SubtractionSolid *pBellWithCut1 = new G4SubtractionSolid("BellWithCut1",
		pBellTubs, pBellCutTubs, 0, G4ThreeVector(0., 0., dBellCutOffsetZ));

	G4Tubs *pBellInTubeCutTubs = new G4Tubs("BellInTubeCutTube", 0.*cm, dBellInTubeCutRadius,
		2*dBellHalfZ, 0.*deg, 360.*deg);

	G4SubtractionSolid *pBellWithCut2 = new G4SubtractionSolid("BellWithCut2",
		pBellWithCut1, pBellInTubeCutTubs, 0, G4ThreeVector(-dLidInTubeOffsetX, 0., 0.));

	m_pBellLogicalVolume = new G4LogicalVolume(pBellWithCut2, SS304LSteel, "BellLogicalVolume", 0, 0, 0);

	m_pBellPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBellOffsetZ),
		m_pBellLogicalVolume, "Bell", m_pLXeLogicalVolume, false, 0);

	//============================= top pmt teflon holder ===========================
	const G4double dGXeHalfZ = 0.5*dGXeHeight;

	const G4double dTopPmtTeflonHolderHalfZ = 0.5*dTopPmtTeflonHolderThickness;
	const G4double dTopPmtTeflonHolderOffsetZ = dGXeHalfZ-dTopPmtTeflonHolderHalfZ-dTopPmtTeflonHolderToBell;

	const G4double dTopPmtTeflonHolderCutHalfX = 0.5*dTopPmtTeflonHolderCutWidth;
	const G4double dTopPmtTeflonHolderCutHalfZ = dTopPmtTeflonHolderCutDepth;

	const G4double dTopPmtTeflonEdgeThickness = dTopPmtTeflonHolderThickness-dTopPmtTeflonHolderCutDepth;
	const G4double dTopPmtTeflonEdgeDelta = 2.*dTopPmtTeflonEdgeThickness/tan(M_PI/180.*dTopPmtTeflonHolderBottomCutAngle);
	const G4double dTopPmtTeflonHolderWindowCutTopHalfX = 0.5*dTopPmtTeflonHolderBottomCutWidth-dTopPmtTeflonEdgeDelta;
	const G4double dTopPmtTeflonHolderWindowCutBottomHalfX = 0.5*dTopPmtTeflonHolderBottomCutWidth+dTopPmtTeflonEdgeDelta;
	const G4double dTopPmtTeflonHolderWindowCutHalfZ = 2*dTopPmtTeflonEdgeThickness;
	
	G4Tubs *pTopPmtTeflonHolderTubs = new G4Tubs("TopPmtTeflonHolderTube", 0.,
		dTopPmtTeflonHolderRadius, dTopPmtTeflonHolderHalfZ, 0.*deg, 360.*deg);

	G4Box *pTopPmtTeflonHolderPmtCut = new G4Box("TopPmtTeflonHolderPmtCut",
		dTopPmtTeflonHolderCutHalfX, dTopPmtTeflonHolderCutHalfX, dTopPmtTeflonHolderCutHalfZ);

	G4Trd *pTopPmtTeflonHolderWindowCut = new G4Trd("TopPmtTeflonHolderWindowCut",
		dTopPmtTeflonHolderWindowCutBottomHalfX, dTopPmtTeflonHolderWindowCutTopHalfX,
		dTopPmtTeflonHolderWindowCutBottomHalfX, dTopPmtTeflonHolderWindowCutTopHalfX,
		dTopPmtTeflonHolderWindowCutHalfZ);

	G4SubtractionSolid *pTopPmtTeflonHolderWithCuts, *pTopPmtTeflonHolderWithPreviousCuts = (G4SubtractionSolid *) pTopPmtTeflonHolderTubs;
	for(G4int iPmtNb=0; iPmtNb<iNbTopPmts; iPmtNb++)
	{
		stringstream hSolidName;
		hSolidName << "TopPmtTeflonHolderWithCut" << iPmtNb << "a";

		G4ThreeVector hPos = GetPmtPosition(iPmtNb, PMT_CASING);
		hPos.setZ(dTopPmtTeflonHolderHalfZ);

		pTopPmtTeflonHolderWithCuts = new G4SubtractionSolid(hSolidName.str(),
			pTopPmtTeflonHolderWithPreviousCuts, pTopPmtTeflonHolderPmtCut,
			GetPmtRotation(iPmtNb, PMT_CASING), hPos);

		hSolidName.str("");
		hSolidName << "TopPmtTeflonHolderWithCut" << iPmtNb << "b";

		hPos.setZ(-dTopPmtTeflonHolderHalfZ+dTopPmtTeflonEdgeThickness);

		pTopPmtTeflonHolderWithPreviousCuts = new G4SubtractionSolid(hSolidName.str(),
			pTopPmtTeflonHolderWithCuts, pTopPmtTeflonHolderWindowCut,
			GetPmtRotation(iPmtNb, PMT_CASING), hPos);
	}

	m_pTopPmtTeflonHolderLogicalVolume = new G4LogicalVolume(pTopPmtTeflonHolderWithPreviousCuts, Teflon,
		"TopPmtTeflonHolderLogicalVolume", 0, 0, 0);

	m_pTopPmtTeflonHolderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopPmtTeflonHolderOffsetZ),
		m_pTopPmtTeflonHolderLogicalVolume, "TopPmtTeflonHolder", m_pGXeLogicalVolume, false, 0);

	//============================= bell pmt support ring ===========================
	const G4double dBellPmtSupportRingInnerRadius = dBellPmtSupportRingRadius-dBellPmtSupportRingThickness;
	const G4double dBellPmtSupportRingHalfZ = 0.5*dBellPmtSupportRingHeight;
	const G4double dBellPmtSupportRingOffsetZ = dBellOffsetZ+dBellHalfZ+dBellPmtSupportRingHalfZ;

	G4Tubs *pBellPmtSupportRingTubs = new G4Tubs("BellPmtSupportRingTube", dBellPmtSupportRingInnerRadius,
		dBellPmtSupportRingRadius, dBellPmtSupportRingHalfZ, 0.*deg, 360.*deg);

	m_pBellPmtSupportRingLogicalVolume = new G4LogicalVolume(pBellPmtSupportRingTubs, SS304LSteel,
		"BellPmtSupportRingLogicalVolume", 0, 0, 0);

	m_pBellPmtSupportRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBellPmtSupportRingOffsetZ),
		m_pBellPmtSupportRingLogicalVolume, "BellPmtSupportRing", m_pLXeLogicalVolume, false, 0);
	

	//> added 24.11.2010, Alex from Elizabeth's code
	//============================= bell support cylinders in LXe ===========================
	const G4double dBellSupportCylinderOffsetX = GetGeometryParameter("BellSupportCylinderOffsetX");
	const G4double dBellSupportCylinderOffsetY = GetGeometryParameter("BellSupportCylinderOffsetY");

	const G4double dVetoGXeHalfZ = 0.5*dLidInnerWallDomeTopToVetoLiquidLevel;
	const G4double dVetoGXeOffsetZ = dLXeHalfZ-dVetoGXeHalfZ;
	const G4double dBellSupportCylinderInLXeHeight = (dVetoGXeOffsetZ-dVetoGXeHalfZ)-(dBellOffsetZ+dBellHalfZ);
    const G4double dBellSupportCylinderInLXeHalfZ = 0.5*dBellSupportCylinderInLXeHeight;
	const G4double dBellSupportCylinderInLXeOffsetZ = dBellOffsetZ+dBellHalfZ+dBellSupportCylinderInLXeHalfZ;

	G4Tubs *pBellSupportCylinderInLXeTubs = new G4Tubs("BellSupportCylinderInLXeTube", 0, dBellSupportCylinderRadius, dBellSupportCylinderInLXeHalfZ, 0.*deg, 360.*deg);

	m_pBellSupportCylinderInLXeLogicalVolume = new G4LogicalVolume(pBellSupportCylinderInLXeTubs, SS304LSteel,"BellSupportCylinderInLXeLogicalVolume", 0, 0, 0);

    //two copies of the logical volume
    m_hBellSupportCylinderInLXePhysicalVolumes.push_back(new G4PVPlacement(0,G4ThreeVector(dBellSupportCylinderOffsetX,dBellSupportCylinderOffsetY,dBellSupportCylinderInLXeOffsetZ), m_pBellSupportCylinderInLXeLogicalVolume, "BellSupportCylinderInLXe", m_pLXeLogicalVolume, false, 0));
    m_hBellSupportCylinderInLXePhysicalVolumes.push_back(new G4PVPlacement(0,G4ThreeVector(dBellSupportCylinderOffsetX,-dBellSupportCylinderOffsetY,dBellSupportCylinderInLXeOffsetZ), m_pBellSupportCylinderInLXeLogicalVolume, "BellSupportCylinderInLXe", m_pLXeLogicalVolume, false, 1));

	//============================= bell support cylinders in (Veto) GXe ===========================
	const G4double dBellSupportCylinderInVetoGXeHeight 	= dLidInnerWallDomeTopToVetoLiquidLevel;
	const G4double dBellSupportCylinderInVetoGXeHalfZ 	= 0.5*dBellSupportCylinderInVetoGXeHeight;
	const G4double dBellSupportCylinderInVetoGXeOffsetZ = dBellSupportCylinderInLXeOffsetZ+dBellSupportCylinderInLXeHalfZ+dBellSupportCylinderInVetoGXeHalfZ;

	const G4double dVetoGXeUpperBulgeRadius 	= dLidInnerWallDomeRadius;
	const G4double dVetoGXeUpperBulgeHeight 	= dLidInnerWallDomeHeight;
	const G4double dVetoGXeUpperBulgeCutOffsetZ = dVetoGXeHalfZ-dVetoGXeUpperBulgeRadius;

	G4Tubs *pBellSupportCylinderInVetoGXeTubs 					= new G4Tubs("BellSupportCylinderInVetoGXeTube", 0, dBellSupportCylinderRadius, dBellSupportCylinderInVetoGXeHalfZ, 0.*deg, 360.*deg);
	G4Sphere *pVetoGXeUpperBulgeCut 							= new G4Sphere("VetoGXeUpperBulgeCut", dVetoGXeUpperBulgeRadius, dVetoGXeUpperBulgeRadius+2*dVetoGXeUpperBulgeHeight, 0.*deg, 360.*deg, 0.*deg, 90.*deg);
	G4SubtractionSolid *pBellSupportCylinderInVetoGXeWithCut1 	= new G4SubtractionSolid("BellSupportCylinderInVetoGXeWithCut1", pBellSupportCylinderInVetoGXeTubs, pVetoGXeUpperBulgeCut, 0, G4ThreeVector(-dBellSupportCylinderOffsetX, -dBellSupportCylinderOffsetY, dVetoGXeUpperBulgeCutOffsetZ));
	G4SubtractionSolid *pBellSupportCylinderInVetoGXeWithCut2  	= new G4SubtractionSolid("BellSupportCylinderInVetoGXeWithCut2", pBellSupportCylinderInVetoGXeTubs, pVetoGXeUpperBulgeCut, 0, G4ThreeVector(-dBellSupportCylinderOffsetX, dBellSupportCylinderOffsetY, dVetoGXeUpperBulgeCutOffsetZ));

	m_pBellSupportCylinderInVetoGXe1LogicalVolume  	= new G4LogicalVolume(pBellSupportCylinderInVetoGXeWithCut1, SS304LSteel, "BellSupportCylinderInVetoGXe1LogicalVolume", 0, 0, 0);
	m_pBellSupportCylinderInVetoGXe2LogicalVolume 	= new G4LogicalVolume(pBellSupportCylinderInVetoGXeWithCut2, SS304LSteel, "BellSupportCylinderInVetoGXe2LogicalVolume", 0, 0, 0);
	m_pBellSupportCylinderInVetoGXe1PhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(dBellSupportCylinderOffsetX,dBellSupportCylinderOffsetY, dBellSupportCylinderInVetoGXeOffsetZ-dVetoGXeOffsetZ), m_pBellSupportCylinderInVetoGXe1LogicalVolume, "BellSupportCylinderInVetoGXe1", m_pVetoGXeLogicalVolume, false, 0);
	m_pBellSupportCylinderInVetoGXe2PhysicalVolume 	= new G4PVPlacement(0,G4ThreeVector(dBellSupportCylinderOffsetX,-dBellSupportCylinderOffsetY,dBellSupportCylinderInVetoGXeOffsetZ-dVetoGXeOffsetZ), m_pBellSupportCylinderInVetoGXe2LogicalVolume, "BellSupportCylinderInVetoGXe2", m_pVetoGXeLogicalVolume, false, 0);	
	//< added 24.11.2010, Alex from Elizabeth's code


	//=============================== top veto angles ===============================
	const G4double dTopVetoAngleHalfX = 0.5*dTopVetoAngleLength;
	const G4double dTopVetoAngleHalfY = 0.5*dTopVetoAngleWidth;
	const G4double dTopVetoAngleHalfZ = 0.5*(dTopVetoAngleLowerHeight+dTopVetoAngleUpperHeight-dTopVetoAngleThickness);
	const G4double dTopVetoAngleOffsetX = dBellPmtSupportRingRadius+dTopVetoAngleHalfX;
	const G4double dTopVetoAngleOffsetZ = dBellPmtSupportRingOffsetZ-dBellPmtSupportRingHalfZ+dTopVetoAngleHalfZ;

	const G4double dTopVetoAngleBottomCutHalfX = dTopVetoAngleLength-dTopVetoAngleThickness;
	const G4double dTopVetoAngleBottomCutHalfZ = dTopVetoAngleLowerHeight-dTopVetoAngleThickness;

	const G4double dTopVetoAngleTopCutHalfX = dTopVetoAngleBottomCutHalfX;
	const G4double dTopVetoAngleTopCutHalfZ = dTopVetoAngleUpperHeight-dTopVetoAngleThickness;

	G4Box *pTopVetoAngleBox 						= new G4Box("TopVetoAngleBox", dTopVetoAngleHalfX, dTopVetoAngleHalfY, dTopVetoAngleHalfZ);
	G4Box *pTopVetoAngleBottomCut 					= new G4Box("TopVetoAngleBottomCut", dTopVetoAngleBottomCutHalfX, dTopVetoAngleWidth, dTopVetoAngleBottomCutHalfZ);
	G4SubtractionSolid *pTopVetoAngleBoxWithCut1 	= new G4SubtractionSolid("TopVetoAngleBoxWithCut1", pTopVetoAngleBox, pTopVetoAngleBottomCut, 0, G4ThreeVector(-dTopVetoAngleHalfX, 0, -dTopVetoAngleHalfZ));
	G4Box *pTopVetoAngleTopCut 						= new G4Box("TopVetoAngleTopCut", dTopVetoAngleTopCutHalfX, dTopVetoAngleWidth, dTopVetoAngleTopCutHalfZ);
	G4SubtractionSolid *pTopVetoAngleBoxWithCut2 	= new G4SubtractionSolid("TopVetoAngleBoxWithCut2", pTopVetoAngleBoxWithCut1, pTopVetoAngleTopCut, 0, G4ThreeVector(dTopVetoAngleHalfX, 0, dTopVetoAngleHalfZ));

	m_pTopVetoAngleLogicalVolume = new G4LogicalVolume(pTopVetoAngleBoxWithCut2, Copper, "TopVetoAngleLogicalVolume", 0, 0, 0);

	for(int iTopVetoAngleId = 0; iTopVetoAngleId < 16; iTopVetoAngleId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(2*iTopVetoAngleId+1, 32, dTopVetoAngleOffsetX);
		hPos.setZ(dTopVetoAngleOffsetZ);

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ((2*iTopVetoAngleId+1)*360.*deg/32);

		m_hTopVetoAnglePhysicalVolumes.push_back(new G4PVPlacement(pRotationMatrix, hPos, m_pTopVetoAngleLogicalVolume, "TopVetoAngle", m_pLXeLogicalVolume, false, 0));
	}

	//============================ upper side veto angles ===========================
	const G4double dUpperSideVetoAngleHalfX = 0.5*dUpperSideVetoAngleLength;
	const G4double dUpperSideVetoAngleHalfY = 0.5*dUpperSideVetoAngleWidth;
	const G4double dUpperSideVetoAngleHalfZ = 0.5*dUpperSideVetoAngleHeight;
	const G4double dUpperSideVetoAngleOffsetX = dBellPmtSupportRingRadius+dUpperSideVetoAngleHalfX;
	const G4double dUpperSideVetoAngleOffsetZ = dBellPmtSupportRingOffsetZ-dBellPmtSupportRingHalfZ+dUpperSideVetoAngleHalfZ;

	const G4double dUpperSideVetoAngleBottomCutHalfX = dUpperSideVetoAngleLength-dUpperSideVetoAngleThickness;
	const G4double dUpperSideVetoAngleBottomCutHalfZ = dUpperSideVetoAngleHeight-dUpperSideVetoAngleThickness;

	G4Box *pUpperSideVetoAngleBox = new G4Box("UpperSideVetoAngleBox",
		dUpperSideVetoAngleHalfX, dUpperSideVetoAngleHalfY, dUpperSideVetoAngleHalfZ);

	G4Box *pUpperSideVetoAngleBottomCut = new G4Box("UpperSideVetoAngleBottomCut",
		dUpperSideVetoAngleBottomCutHalfX, dUpperSideVetoAngleWidth, dUpperSideVetoAngleBottomCutHalfZ);

	G4SubtractionSolid *pUpperSideVetoAngleBoxWithCut1 = new G4SubtractionSolid("UpperSideVetoAngleBoxWithCut1",
		pUpperSideVetoAngleBox, pUpperSideVetoAngleBottomCut,
		0, G4ThreeVector(-dUpperSideVetoAngleHalfX, 0, -dUpperSideVetoAngleHalfZ));

	m_pUpperSideVetoAngleLogicalVolume = new G4LogicalVolume(pUpperSideVetoAngleBoxWithCut1,
		Copper, "UpperSideVetoAngleLogicalVolume", 0, 0, 0);

	for(int iUpperSideVetoAngleId = 0; iUpperSideVetoAngleId < 16; iUpperSideVetoAngleId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(2*iUpperSideVetoAngleId, 32, dUpperSideVetoAngleOffsetX);
		hPos.setZ(dUpperSideVetoAngleOffsetZ);

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ((2*iUpperSideVetoAngleId)*360.*deg/32);

		m_hUpperSideVetoAnglePhysicalVolumes.push_back(new G4PVPlacement(pRotationMatrix, hPos,
			m_pUpperSideVetoAngleLogicalVolume, "UpperSideVetoAngle", m_pLXeLogicalVolume, false, 0));
	}

	//================================== attributes =================================
	//G4Colour hSteelColour(0.0, 0.0, 1.0);		// xlblue
	//G4Colour hCopperColour(1.0, 1.0, 0.0);		// yellow
	G4Colour hSteelColour(0.0, 0.0, 0.0);		// xlblue
	G4Colour hCopperColour(0.0, 0.0, 0.0);		// yellow
	
	G4VisAttributes *pBellVisAtt = new G4VisAttributes(hSteelColour);
	pBellVisAtt->SetVisibility(true);
	m_pBellLogicalVolume->SetVisAttributes(pBellVisAtt);
	m_pBellPmtSupportRingLogicalVolume->SetVisAttributes(pBellVisAtt);
	m_pBellSupportCylinderInLXeLogicalVolume->SetVisAttributes(pBellVisAtt);
	m_pBellSupportCylinderInVetoGXe1LogicalVolume->SetVisAttributes(pBellVisAtt);
	m_pBellSupportCylinderInVetoGXe2LogicalVolume->SetVisAttributes(pBellVisAtt);

	G4VisAttributes *pTopVetoAngleVisAtt = new G4VisAttributes(hCopperColour);
	pTopVetoAngleVisAtt->SetVisibility(true);
	m_pTopVetoAngleLogicalVolume->SetVisAttributes(pTopVetoAngleVisAtt);

	G4VisAttributes *pUpperSideVetoAngleVisAtt = new G4VisAttributes(hCopperColour);
	pUpperSideVetoAngleVisAtt->SetVisibility(true);
	m_pUpperSideVetoAngleLogicalVolume->SetVisAttributes(pUpperSideVetoAngleVisAtt);
}

void
Xenon100DetectorConstruction::ConstructFieldCage()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< field cage >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dGXeHeight 	= GetGeometryParameter("GXeHeight");
	const G4double dLXeHeight 	= GetGeometryParameter("LXeHeight");

	const G4double dLidInnerWallDomeHeight 		= GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInTubeHeightInside 		= GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference 	= GetGeometryParameter("CryostatLidInTubeOffsetDifference");

	const G4double dBellHeight 				= GetGeometryParameter("BellHeight");
	const G4double dTopPlateThickness 		= GetGeometryParameter("TopPlateThickness");
	const G4double dTeflonPanelHeight 		= GetGeometryParameter("TeflonPanelHeight");
	const G4double dBottomPlateThickness 	= GetGeometryParameter("BottomPlateThickness");

	const G4double dGridToAnode 		= GetGeometryParameter("GridToAnode");
	const G4double dAnodeToLiquidLevel 	= GetGeometryParameter("AnodeToLiquidLevel");
	const G4double dLiquidLevelToGrid 	= GetGeometryParameter("LiquidLevelToGrid");

	const G4double dTopGridRingRadius 		= GetGeometryParameter("TopGridRingRadius");
	const G4double dAnodeGridRingRadius 	= GetGeometryParameter("AnodeGridRingRadius");
	const G4double dBottomGridRingRadius 	= GetGeometryParameter("BottomGridRingRadius");
	const G4double dCathodeGridRingRadius 	= GetGeometryParameter("CathodeGridRingRadius");
	const G4double dGridRingWidth 			= GetGeometryParameter("GridRingWidth");
	const G4double dGridRingThickness 		= GetGeometryParameter("GridRingThickness");
	const G4double dGridMeshThickness 		= GetGeometryParameter("GridMeshThickness");
    const G4double dScreenMeshRadius 		= GetGeometryParameter("LowerTeflonPanelRadius")-GetGeometryParameter("LowerTeflonPanelThickness");
    const G4double dCathodeToScreenMesh 	= GetGeometryParameter("CathodeToScreenMesh");

	const G4double dGridHolderRadius 	= GetGeometryParameter("GridHolderRadius");
	const G4double dGridHolderWidth 	= GetGeometryParameter("GridHolderWidth");
	const G4double dGridHolderCutHeight = GetGeometryParameter("GridHolderCutHeight");

	const G4double dCopperScreenThickness 	= GetGeometryParameter("CopperScreenThickness");
	const G4double dCopperScreenOuterRadius = GetGeometryParameter("CopperScreenOuterRadius");
	const G4double dCopperScreenInnerRadius = GetGeometryParameter("TopGridRingRadius");

	G4Material *GridMeshAluminium 	= G4Material::GetMaterial("GridMeshAluminium");
	G4Material *SS304LSteel 		= G4Material::GetMaterial("SS304LSteel");
	G4Material *Teflon 				= G4Material::GetMaterial("Teflon");
	G4Material *Copper 				= G4Material::GetMaterial("Copper");

	//================================== top grid ===================================

	//------------------------------- top grid mesh ---------------------------------
	const G4double dTopGridMeshOffsetZ = -0.5*dGXeHeight+dAnodeToLiquidLevel+dGridToAnode;
	
	G4Tubs *pTopGridMeshTubs = new G4Tubs("TopGridMeshTube", 0., dTopGridRingRadius,
		0.5*dGridMeshThickness, 0.*deg, 360.*deg);

	m_pTopGridMeshLogicalVolume = new G4LogicalVolume(pTopGridMeshTubs, GridMeshAluminium, "TopGridMeshLogicalVolume");

	m_pTopGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopGridMeshOffsetZ),
		m_pTopGridMeshLogicalVolume, "TopGridMesh", m_pGXeLogicalVolume, false, 0);

	//------------------------------- top grid ring ---------------------------------
	const G4double dTopGridRingOffsetZ = dTopGridMeshOffsetZ+0.5*dGridRingThickness;

	G4Tubs *pTopGridRingTubs = new G4Tubs("TopGridRingTube", dTopGridRingRadius,
		dTopGridRingRadius+dGridRingWidth, 0.5*dGridRingThickness, 0.*deg, 360.*deg);

	m_pTopGridRingLogicalVolume = new G4LogicalVolume(pTopGridRingTubs, SS304LSteel, "TopGridRingLogicalVolume");

	m_pTopGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopGridRingOffsetZ),
		m_pTopGridRingLogicalVolume, "TopGridRing", m_pGXeLogicalVolume, false, 0);

	//==================================== anode ====================================
	
	//--------------------------------- anode mesh ----------------------------------
	const G4double dAnodeGridMeshOffsetZ = -0.5*dGXeHeight+dAnodeToLiquidLevel;
	
	G4Tubs *pAnodeGridMeshTubs = new G4Tubs("AnodeGridMeshTube", 0., dAnodeGridRingRadius,
		0.5*dGridMeshThickness, 0.*deg, 360.*deg);

	m_pAnodeGridMeshLogicalVolume = new G4LogicalVolume(pAnodeGridMeshTubs, GridMeshAluminium, "AnodeGridMeshLogicalVolume");

	m_pAnodeGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dAnodeGridMeshOffsetZ),
		m_pAnodeGridMeshLogicalVolume, "AnodeGridMesh", m_pGXeLogicalVolume, false, 0);

	//--------------------------------- anode ring ----------------------------------
	const G4double dAnodeGridRingOffsetZ = dAnodeGridMeshOffsetZ+0.5*dGridRingThickness;

	G4Tubs *pAnodeGridRingTubs = new G4Tubs("AnodeGridRingTube", dAnodeGridRingRadius,
		dAnodeGridRingRadius+dGridRingWidth, 0.5*dGridRingThickness, 0.*deg, 360.*deg);

	m_pAnodeGridRingLogicalVolume = new G4LogicalVolume(pAnodeGridRingTubs, SS304LSteel, "AnodeGridRingLogicalVolume");

	m_pAnodeGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dAnodeGridRingOffsetZ),
		m_pAnodeGridRingLogicalVolume, "AnodeGridRing", m_pGXeLogicalVolume, false, 0);

	//================================= bottom grid =================================

	//------------------------------ bottom grid mesh -------------------------------
	const G4double dLXeHalfZ = 0.5*dLXeHeight;	
	const G4double dBottomGridMeshOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight+dGridRingThickness;
	
	G4Tubs *pBottomGridMeshTubs = new G4Tubs("BottomGridMeshTube", 0., dBottomGridRingRadius,
		0.5*dGridMeshThickness, 0.*deg, 360.*deg);

	m_pBottomGridMeshLogicalVolume = new G4LogicalVolume(pBottomGridMeshTubs, GridMeshAluminium, "BottomGridMeshLogicalVolume");

	m_pBottomGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomGridMeshOffsetZ),
		m_pBottomGridMeshLogicalVolume, "BottomGridMesh", m_pLXeLogicalVolume, false, 0);

	//------------------------------ bottom grid ring -------------------------------
	const G4double dBottomGridRingOffsetZ = dBottomGridMeshOffsetZ-0.5*dGridRingThickness;

	G4Tubs *pBottomGridRingTubs = new G4Tubs("BottomGridRingTube", dBottomGridRingRadius,
		dBottomGridRingRadius+dGridRingWidth, 0.5*dGridRingThickness, 0.*deg, 360.*deg);

	m_pBottomGridRingLogicalVolume = new G4LogicalVolume(pBottomGridRingTubs, SS304LSteel, "BottomGridRingLogicalVolume");

	m_pBottomGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomGridRingOffsetZ),
		m_pBottomGridRingLogicalVolume, "BottomGridRing", m_pLXeLogicalVolume, false, 0);

	//=================================== cathode ===================================

	//-------------------------------- cathode mesh ---------------------------------
	const G4double dCathodeGridMeshOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight-dTopPlateThickness-dTeflonPanelHeight-dBottomPlateThickness+dGridRingThickness;
	
	G4Tubs *pCathodeGridMeshTubs = new G4Tubs("CathodeGridMeshTube", 0., dCathodeGridRingRadius,
		0.5*dGridMeshThickness, 0.*deg, 360.*deg);

	m_pCathodeGridMeshLogicalVolume = new G4LogicalVolume(pCathodeGridMeshTubs, GridMeshAluminium, "CathodeGridMeshLogicalVolume");

	m_pCathodeGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dCathodeGridMeshOffsetZ),
		m_pCathodeGridMeshLogicalVolume, "CathodeGridMesh", m_pLXeLogicalVolume, false, 0);

	//-------------------------------- cathode ring ---------------------------------
	const G4double dCathodeGridRingOffsetZ = dCathodeGridMeshOffsetZ-0.5*dGridRingThickness;

	G4Tubs *pCathodeGridRingTubs = new G4Tubs("CathodeGridRingTube", dCathodeGridRingRadius,
		dCathodeGridRingRadius+dGridRingWidth, 0.5*dGridRingThickness, 0.*deg, 360.*deg);

	m_pCathodeGridRingLogicalVolume = new G4LogicalVolume(pCathodeGridRingTubs, SS304LSteel, "CathodeGridRingLogicalVolume");

	m_pCathodeGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dCathodeGridRingOffsetZ),
		m_pCathodeGridRingLogicalVolume, "CathodeGridRing", m_pLXeLogicalVolume, false, 0);


	//=================================== screening mesh ===================================
	//-------------------------------- screen mesh ---------------------------------
	const G4double dScreenMeshOffsetZ = dCathodeGridMeshOffsetZ - dCathodeToScreenMesh;	
        G4cout<<"dScreenMeshOffsetZ: "<<dScreenMeshOffsetZ<<G4endl;
	G4Tubs *pScreenMeshTubs = new G4Tubs("ScreenMeshTube", 0., dScreenMeshRadius, 0.5*dGridMeshThickness, 0.*deg, 360.*deg);

	m_pScreenMeshLogicalVolume = new G4LogicalVolume(pScreenMeshTubs, GridMeshAluminium, "ScreenMeshLogicalVolume");
	m_pScreenMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dScreenMeshOffsetZ),
                                                        m_pScreenMeshLogicalVolume, "ScreenMesh", m_pLXeLogicalVolume, false, 0);

	//================================= grid holder =================================

	//----------------------------- grid holder in lxe ------------------------------
	const G4double dGridHolderInLXeHeight = dLiquidLevelToGrid+dGridRingThickness;
	const G4double dGridHolderInLXeOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight+0.5*dGridHolderInLXeHeight;

	G4Tubs *pGridHolderInLXeTubs 	= new G4Tubs("GridHolderInLXeTube", dGridHolderRadius, dGridHolderRadius+dGridHolderWidth, 0.5*dGridHolderInLXeHeight, 0.*deg, 360.*deg);
	G4Tubs *pGridHolderInLXeCut 	= new G4Tubs("GridHolderInLXeCut", 0., dBottomGridRingRadius, dGridHolderCutHeight, 0.*deg, 360.*deg);
	G4SubtractionSolid *pGridHolderInLXeWithCut = new G4SubtractionSolid("GridHolderInLXeWithCut", pGridHolderInLXeTubs, pGridHolderInLXeCut, 0, G4ThreeVector(0., 0., -0.5*dGridHolderInLXeHeight));

	m_pGridHolderInLXeLogicalVolume = new G4LogicalVolume(pGridHolderInLXeWithCut, Teflon, "GridHolderInLXeLogicalVolume");
	m_pGridHolderInLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dGridHolderInLXeOffsetZ), m_pGridHolderInLXeLogicalVolume, "GridHolderInLXe", m_pLXeLogicalVolume, false, 0);

	//----------------------------- grid holder in gxe ------------------------------
	const G4double dGXeHalfZ = 0.5*dGXeHeight;	
	const G4double dGridHolderInGXeHeight = dAnodeToLiquidLevel+dGridToAnode;
	const G4double dGridHolderInGXeOffsetZ = -dGXeHalfZ+0.5*dGridHolderInGXeHeight;
	const G4double dGridHolderInGXeCut1OffsetZ = 0.5*dGridHolderInGXeHeight-dGridToAnode-dGridRingThickness+0.5*dGridHolderCutHeight;
	const G4double dGridHolderInGXeCut2OffsetZ = 0.5*dGridHolderInGXeHeight;

	G4Tubs *pGridHolderInGXeTubs = new G4Tubs("GridHolderInGXeTube", dGridHolderRadius,
		dGridHolderRadius+dGridHolderWidth, 0.5*dGridHolderInGXeHeight, 0.*deg, 360.*deg);

	G4Tubs *pGridHolderInGXeCut1 = new G4Tubs("GridHolderInGXeCut1", 0., dBottomGridRingRadius,
		0.5*dGridHolderCutHeight, 0.*deg, 360.*deg);

	G4SubtractionSolid *pGridHolderInGXeWithCut1 = new G4SubtractionSolid("GridHolderInGXeWithCut1",
		pGridHolderInGXeTubs, pGridHolderInGXeCut1, 0, G4ThreeVector(0., 0., dGridHolderInGXeCut1OffsetZ));

	G4Tubs *pGridHolderInGXeCut2 = new G4Tubs("GridHolderInGXeCut2", 0., dBottomGridRingRadius,
		dGridHolderCutHeight-dGridRingThickness, 0.*deg, 360.*deg);

	G4SubtractionSolid *pGridHolderInGXeWithCut2 = new G4SubtractionSolid("GridHolderInGXeWithCut2",
		pGridHolderInGXeWithCut1, pGridHolderInGXeCut2, 0, G4ThreeVector(0., 0., dGridHolderInGXeCut2OffsetZ));

	m_pGridHolderInGXeLogicalVolume = new G4LogicalVolume(pGridHolderInGXeWithCut2, Teflon, "GridHolderInGXeLogicalVolume");

	m_pGridHolderInGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dGridHolderInGXeOffsetZ),
		m_pGridHolderInGXeLogicalVolume, "GridHolderInGXe", m_pGXeLogicalVolume, false, 0);

        //----------------------------- copper screen on top of the grids ---------------
    const G4double dCopperScreenOffsetZ = dTopGridRingOffsetZ + 0.5*dGridRingThickness + 0.5*dCopperScreenThickness;
    G4Tubs *pCopperScreenTube 		= new G4Tubs("CopperScreenTube", dCopperScreenInnerRadius, dCopperScreenOuterRadius, 0.5*dCopperScreenThickness, 0.*deg, 360.*deg);
    m_pCopperScreenLogicalVolume 	= new G4LogicalVolume(pCopperScreenTube, Copper, "CopperScreenLogicalVolume");
    m_pCopperScreenPhysicalVolume 	= new G4PVPlacement(0, G4ThreeVector(0., 0., dCopperScreenOffsetZ), m_pCopperScreenLogicalVolume, "CopperScreen", m_pGXeLogicalVolume, false, 0);

	//=============================== optical surfaces ==============================
	G4double dSigmaAlpha = 0.1;
	G4OpticalSurface *pTeflonOpticalSurface = new G4OpticalSurface("TeflonOpticalSurface", unified, ground, dielectric_metal, dSigmaAlpha);
    G4OpticalSurface *pCopperOpticalSurface = new G4OpticalSurface("CopperOpticalSurface", unified, ground, dielectric_metal, 0);

	pTeflonOpticalSurface->SetMaterialPropertiesTable(Teflon->GetMaterialPropertiesTable());
	pCopperOpticalSurface->SetMaterialPropertiesTable(Copper->GetMaterialPropertiesTable());

	new G4LogicalBorderSurface("TeflonGridHolderInLXeLogicalBorderSurface", m_pGXePhysicalVolume, m_pGridHolderInLXePhysicalVolume, pTeflonOpticalSurface);
	new G4LogicalBorderSurface("TeflonGridHolderInGXeLogicalBorderSurface", m_pGXePhysicalVolume, m_pGridHolderInGXePhysicalVolume, pTeflonOpticalSurface);
    new G4LogicalBorderSurface("CopperScreenLogicalBorderSurface", m_pGXePhysicalVolume, m_pCopperScreenPhysicalVolume, pCopperOpticalSurface);

	//================================== attributes =================================
	//G4Colour hGridMeshColor(0.42, 0.72, 0.0, 1.);
	//G4Colour hGridRingColor(0.04, 0.72, 0.0, 1.);
	//G4Colour hGridColour(0.5, 0.5, 0.5);		// grey
	//G4Colour hCopperColour(1.0, 1.0, 0.0);		// yellow
	G4Colour hGridColour(0.0, 0.0, 0.0);		// grey
	G4Colour hCopperColour(0.0, 0.0, 0.0);		// yellow

	G4VisAttributes *pGridMeshVisAtt = new G4VisAttributes(hGridColour);
	pGridMeshVisAtt->SetVisibility(true);
	m_pTopGridMeshLogicalVolume		->SetVisAttributes(pGridMeshVisAtt);
	m_pAnodeGridMeshLogicalVolume	->SetVisAttributes(pGridMeshVisAtt);
	m_pBottomGridMeshLogicalVolume	->SetVisAttributes(pGridMeshVisAtt);
	m_pCathodeGridMeshLogicalVolume	->SetVisAttributes(pGridMeshVisAtt);
	m_pScreenMeshLogicalVolume		->SetVisAttributes(pGridMeshVisAtt);

	G4VisAttributes *pGridRingVisAtt = new G4VisAttributes(hGridColour);
	pGridRingVisAtt->SetVisibility(true);
	m_pTopGridRingLogicalVolume		->SetVisAttributes(pGridRingVisAtt);
	m_pAnodeGridRingLogicalVolume	->SetVisAttributes(pGridRingVisAtt);
	m_pBottomGridRingLogicalVolume	->SetVisAttributes(pGridRingVisAtt);
	m_pCathodeGridRingLogicalVolume	->SetVisAttributes(pGridRingVisAtt);

	G4VisAttributes *pCopperScreenVisAtt = new G4VisAttributes(hCopperColour);
	pCopperScreenVisAtt->SetVisibility(true);
	m_pCopperScreenLogicalVolume	->SetVisAttributes(pCopperScreenVisAtt);

	//OffsetZ
	G4cout<<"OffsetZ: "<< G4endl;
	G4cout<<"---------------------------------------------------"<< G4endl;
	G4cout<<"TopGridMesh     :     "<<dTopGridMeshOffsetZ<<G4endl;
	G4cout<<"AnodeGridMesh   :   "<<dAnodeGridMeshOffsetZ<<G4endl;
	G4cout<<"BottomGridMesh  :  "<<dBottomGridMeshOffsetZ<<G4endl;
	G4cout<<"CathodeGridMesh : "<<dCathodeGridMeshOffsetZ<<G4endl;
	G4cout<<"ScreenMesh      :      "<<dScreenMeshOffsetZ<<G4endl;
	G4cout<<"---------------------------------------------------: "<< G4endl;
}

void
Xenon100DetectorConstruction::ConstructTPC()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< tpc >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

	const G4double dPmtWidth = GetGeometryParameter("PmtWidth");
	const G4double dPmtSpacing = GetGeometryParameter("PmtSpacing");

	const G4double dCathodeGridRingRadius = GetGeometryParameter("CathodeGridRingRadius");
	const G4double dGridRingWidth = GetGeometryParameter("GridRingWidth");
	const G4double dGridRingThickness = GetGeometryParameter("GridRingThickness");

	const G4double dVesselInnerWallDomeHeight = GetGeometryParameter("CryostatVesselInnerWallDomeHeight");
	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInnerWallDomeTopToVetoLiquidLevel = GetGeometryParameter("CryostatLidInnerWallDomeTopToVetoLiquidLevel");

	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");

	const G4double dBellHeight = GetGeometryParameter("BellHeight");

	const G4double dTopPlateThickness = GetGeometryParameter("TopPlateThickness");
	const G4double dTopPlateRadius = GetGeometryParameter("TopPlateRadius");
	const G4double dTopPlateWidth = GetGeometryParameter("TopPlateWidth");

	const G4double dTeflonPanelRadius = GetGeometryParameter("TeflonPanelRadius");
	const G4double dTeflonPanelThickness = GetGeometryParameter("TeflonPanelThickness");
	const G4double dTeflonPanelHeight = GetGeometryParameter("TeflonPanelHeight");

	const G4double dTeflonRodsRadius = GetGeometryParameter("TeflonRodsRadius");
	const G4double dTeflonRodsHeight = GetGeometryParameter("TeflonRodsHeight");
	const G4double dTeflonRodsOffsetX = GetGeometryParameter("TeflonRodsOffsetX");

	const G4double dBottomPlateThickness = GetGeometryParameter("BottomPlateThickness");
	const G4double dBottomPlateRadius = GetGeometryParameter("BottomPlateRadius");
	const G4double dBottomPlateWidth = GetGeometryParameter("BottomPlateWidth");
	const G4double dBottomPlateAnchorLength = GetGeometryParameter("BottomPlateAnchorLength");
	const G4double dBottomPlateAnchorWidth = GetGeometryParameter("BottomPlateAnchorWidth");

	const G4double dUpperBasePlateThickness = GetGeometryParameter("UpperBasePlateThickness");
	const G4double dUpperBasePlateRadius = GetGeometryParameter("UpperBasePlateRadius");
	const G4double dUpperBasePlateWidth = GetGeometryParameter("UpperBasePlateWidth");
	const G4double dUpperBasePlateAnchorLength = GetGeometryParameter("UpperBasePlateAnchorLength");
	const G4double dUpperBasePlateAnchorWidth = GetGeometryParameter("UpperBasePlateAnchorWidth");

	const G4double dLowerTeflonPanelRadius = GetGeometryParameter("LowerTeflonPanelRadius");
	const G4double dLowerTeflonPanelThickness = GetGeometryParameter("LowerTeflonPanelThickness");
	const G4double dLowerTeflonPanelHeight = GetGeometryParameter("LowerTeflonPanelHeight");

	const G4double dLowerBasePlateRadius = GetGeometryParameter("LowerBasePlateRadius");
	const G4double dLowerBasePlateThickness = GetGeometryParameter("LowerBasePlateThickness");
	const G4double dLowerBasePlateSupportRingRadius = GetGeometryParameter("LowerBasePlateSupportRingRadius");

	const G4double dBottomPmtPlateRadius = GetGeometryParameter("BottomPmtPlateRadius");
	const G4double dBottomPmtPlateThickness = GetGeometryParameter("BottomPmtPlateThickness");

	const G4double dBottomVetoAngleThickness = GetGeometryParameter("BottomVetoAngleThickness");
	const G4double dBottomVetoAngleHeight = GetGeometryParameter("BottomVetoAngleHeight");
	const G4double dBottomVetoAngleLength = GetGeometryParameter("BottomVetoAngleLength");
	const G4double dBottomVetoAngleWidth = GetGeometryParameter("BottomVetoAngleWidth");

	const G4double dLowerSideVetoAngleThickness = GetGeometryParameter("LowerSideVetoAngleThickness");
	const G4double dLowerSideVetoAngleHeight = GetGeometryParameter("LowerSideVetoAngleHeight");
	const G4double dLowerSideVetoAngleShortLength = GetGeometryParameter("LowerSideVetoAngleShortLength");
	const G4double dLowerSideVetoAngleLongLength = GetGeometryParameter("LowerSideVetoAngleLongLength");
	const G4double dLowerSideVetoAngleWidth = GetGeometryParameter("LowerSideVetoAngleWidth");

	const G4double dTeflonSideVetoLiningRadius = GetGeometryParameter("TeflonSideVetoLiningRadius");
	const G4double dTeflonSideVetoLiningThickness = GetGeometryParameter("TeflonSideVetoLiningThickness");

 	const G4double dTeflonLedBlockWidth = GetGeometryParameter("TeflonLedBlockWidth");
 	const G4double dTeflonLedBlockHeight = GetGeometryParameter("TeflonLedBlockHeight");
 	const G4double dTeflonLedBlockDepth = GetGeometryParameter("TeflonLedBlockDepth");
 	const G4double dTeflonLedBlockOffsetX = GetGeometryParameter("TeflonLedBlockOffsetX");
 	const G4double dTeflonLedBlockOffsetY = GetGeometryParameter("TeflonLedBlockOffsetY");
 	const G4double dTeflonLedBlockRotationZ = GetGeometryParameter("TeflonLedBlockRotationZ");

 	const G4double dTeflonBottomReflectorThickness = GetGeometryParameter("TeflonBottomReflectorThickness");
 	const G4double dTeflonBottomReflectorTopRadius = GetGeometryParameter("TeflonBottomReflectorTopRadius");
 	const G4double dTeflonBottomReflectorBottomRadius = GetGeometryParameter("TeflonBottomReflectorBottomRadius");
 	const G4double dBottomReflectorHolderRadius = GetGeometryParameter("BottomReflectorHolderRadius");
 	const G4double dBottomReflectorHolderHeight = GetGeometryParameter("BottomReflectorHolderHeight");
 	const G4double dBottomReflectorHolderOffsetX = GetGeometryParameter("BottomReflectorHolderOffsetX");

	G4Material *Copper = G4Material::GetMaterial("Copper");	
	G4Material *Teflon = G4Material::GetMaterial("Teflon");

	//================================== top plate ==================================
	const G4double dLXeHalfZ = 0.5*dLXeHeight;
	const G4double dTopPlateHalfZ = 0.5*dTopPlateThickness;
	const G4double dTopPlateOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight-dTopPlateHalfZ;

	G4Tubs *pTopPlateTubs = new G4Tubs("TopPlateTube", dTopPlateRadius-dTopPlateWidth,
		dTopPlateRadius, dTopPlateHalfZ, 0.*deg, 360.*deg);

	m_pTopPlateLogicalVolume = new G4LogicalVolume(pTopPlateTubs,
		Copper, "TopPlateLogicalVolume", 0, 0, 0);

	m_pTopPlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopPlateOffsetZ),
		m_pTopPlateLogicalVolume, "TopPlate", m_pLXeLogicalVolume, false, 0);

	//================================= teflon panels ===============================
	const G4double dTeflonPanelHalfZ = 0.5*dTeflonPanelHeight;
	const G4double dTeflonPanelOffsetZ = dTopPlateOffsetZ-dTopPlateHalfZ-dTeflonPanelHalfZ;

	const G4double pdTeflonPanelZPlanes[] = {-dTeflonPanelHalfZ, dTeflonPanelHalfZ};
	G4double pdTeflonPanelInnerRadii[24], pdTeflonPanelOuterRadii[24];
	
	for(int i = 0; i < 24; i++)
	{
		pdTeflonPanelInnerRadii[i] = dTeflonPanelRadius;
		pdTeflonPanelOuterRadii[i] = dTeflonPanelRadius+dTeflonPanelThickness;
	}

	G4Polyhedra *pTeflonPanelPolyhedra = new G4Polyhedra("TeflonPanelPolyhedra", 0.*deg, 360.*deg,
		24, 2, pdTeflonPanelZPlanes, pdTeflonPanelInnerRadii, pdTeflonPanelOuterRadii);

	m_pTeflonPanelLogicalVolume = new G4LogicalVolume(pTeflonPanelPolyhedra, Teflon,
			"TeflonPanelLogicalVolume", 0, 0, 0);

	m_pTeflonPanelPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTeflonPanelOffsetZ),
		m_pTeflonPanelLogicalVolume, "TeflonPanel", m_pLXeLogicalVolume, false, 0);

	//====================================== rods ===================================
	const G4double dTeflonRodsHalfZ = 0.5*dTeflonRodsHeight;
	G4Tubs *pTeflonRodsTubs = new G4Tubs("TeflonRodsTube", 0, dTeflonRodsRadius, dTeflonRodsHalfZ, 0.*deg, 360.*deg);

	m_pTeflonRodsLogicalVolume = new G4LogicalVolume(pTeflonRodsTubs, Teflon,
			"TeflonRodsLogicalVolume", 0, 0, 0);

	for(int iRodId = 0; iRodId < 16; iRodId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(iRodId, 16, dTeflonRodsOffsetX);
		hPos.setZ(dTeflonPanelOffsetZ);
		hPos.rotateZ(360.*deg/32);

		m_hTeflonRodsPhysicalVolumes.push_back(new G4PVPlacement(0, hPos,
			m_pTeflonRodsLogicalVolume, "TeflonRods", m_pLXeLogicalVolume, false, 0));
	}
	
	//==============================================================================
	//================================= resistors ==================================
	G4Material *Cirlex = G4Material::GetMaterial("Ceramics");
	const G4double dResistorChainHalfZ		= 0.5*dTeflonRodsHeight;
	const G4double dResistorChainRadius		= 1*mm;
	const G4double dResistorChainOffsetZ	= dTeflonPanelOffsetZ;
	const G4double dResistorChainOffsetX	= dTeflonPanelRadius + dTeflonPanelThickness + 2*mm;
	G4Tubs *pResistorChain = new G4Tubs("ResistorChainTube", 0, dResistorChainRadius, dResistorChainHalfZ, 0.*deg, 360.*deg);
	m_pResistorChainLogicalVolume = new G4LogicalVolume(pResistorChain, Cirlex, "ResistorChainLogicalVolume", 0, 0 ,0);
	m_pResistorChainPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dResistorChainOffsetX, 0., dResistorChainOffsetZ), 
	m_pResistorChainLogicalVolume, "ResistorChain", m_pLXeLogicalVolume, false, 0); 
	
	G4Colour hResistorChainColour (1.0, 0.0, 0.0);

	G4VisAttributes *pResistorChainVisAtt = new G4VisAttributes(hResistorChainColour);
	pResistorChainVisAtt->SetVisibility(true);
	m_pResistorChainLogicalVolume->SetVisAttributes(pResistorChainVisAtt);
	
/*	//==============================================================================
	//================================== cables (bottom veto bundle) ===============
	//G4Material *Copper = G4Material::GetMaterial("Copper");
	const G4double dCableHalfZ		= 0.5*dTeflonRodsHeight;
	const G4double dCableRadius		= 8.0 *mm;
	const G4double dCableOffsetZ	= dTeflonPanelOffsetZ;
	const G4double dCableOffsetX	= 190.0 *mm;
	G4Tubs *pCable = new G4Tubs("CableTube", 0, dCableRadius, dCableHalfZ, 0.*deg, 360.*deg);
	m_pCableLogicalVolume = new G4LogicalVolume(pCable, Copper, "CableLogicalVolume", 0, 0, 0);
	m_pCablePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dCableOffsetX, 0., dCableOffsetZ),
	m_pCableLogicalVolume, "Cable", m_pLXeLogicalVolume, false, 0);
	
	G4Colour hCableColour (1.0, 1.0, 0.0);
	
	G4VisAttributes *pCableVisAtt = new G4VisAttributes(hCableColour);
	pCableVisAtt->SetVisibility(true);
	m_pCableLogicalVolume->SetVisAttributes(pCableVisAtt);
*/	
	//================================= bottom plate ===============================
	const G4double dBottomPlateHalfZ = 0.5*dBottomPlateThickness;
	const G4double dBottomPlateInnerRadius = dBottomPlateRadius-dBottomPlateWidth;
	const G4double dBottomPlateOuterRadius = dBottomPlateRadius+dBottomPlateAnchorLength;
	const G4double dBottomPlateOffsetZ = dTeflonPanelOffsetZ-dTeflonPanelHalfZ-dBottomPlateHalfZ;

	const G4double dBottomPlateCutInnerRadius = dBottomPlateRadius;
	const G4double dBottomPlateCutOuterRadius = dBottomPlateRadius+2*dBottomPlateAnchorLength;
	const G4double dBottomPlateCutAngularSize = (-dBottomPlateAnchorWidth/dBottomPlateRadius*180./M_PI+360./16)*deg;
	const G4double dBottomPlateCutHalfZ = dBottomPlateThickness;

	G4Tubs *pBottomPlateTubs = new G4Tubs("BottomPlateTube",
		dBottomPlateInnerRadius, dBottomPlateOuterRadius, dBottomPlateHalfZ, 0.*deg, 360.*deg);

	G4Tubs *pBottomPlateSideCut = new G4Tubs("BottomPlateSideCut",
		dBottomPlateCutInnerRadius, dBottomPlateCutOuterRadius, dBottomPlateCutHalfZ,
		0, dBottomPlateCutAngularSize);

	G4SubtractionSolid *pBottomPlateWithCuts, *pBottomPlateWithPreviousCuts = (G4SubtractionSolid *) pBottomPlateTubs;
	for(int iCutId = 0; iCutId < 16; iCutId++)
	{
		stringstream hStr;
		hStr << "BottomPlateTubeWithCut" << iCutId+1;

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ(iCutId*360.*deg/16+0.5*dBottomPlateCutAngularSize);

		pBottomPlateWithCuts = new G4SubtractionSolid(hStr.str(),
			pBottomPlateWithPreviousCuts, pBottomPlateSideCut,
			pRotationMatrix, G4ThreeVector(0., 0., 0.));

		pBottomPlateWithPreviousCuts = pBottomPlateWithCuts;
	}

	G4Tubs *pBottomPlateCathodeCut = new G4Tubs("BottomPlateCathodeCut",
		0, dCathodeGridRingRadius+dGridRingWidth, dGridRingThickness, 0.*deg, 360.*deg);

	G4SubtractionSolid *pBottomPlateWithAllCuts = new G4SubtractionSolid("BottomPlateTubeWithAllCuts",
		pBottomPlateWithCuts, pBottomPlateCathodeCut, 0, G4ThreeVector(0., 0., -dBottomPlateHalfZ));
	
	m_pBottomPlateLogicalVolume = new G4LogicalVolume(pBottomPlateWithAllCuts,
		Copper, "BottomPlateLogicalVolume", 0, 0, 0);

	m_pBottomPlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomPlateOffsetZ),
		m_pBottomPlateLogicalVolume, "BottomPlate", m_pLXeLogicalVolume, false, 0);

	//=============================== upper base plate ==============================
	const G4double dUpperBasePlateHalfZ = 0.5*dUpperBasePlateThickness;
	const G4double dUpperBasePlateInnerRadius = dUpperBasePlateRadius-dUpperBasePlateWidth;
	const G4double dUpperBasePlateOuterRadius = dUpperBasePlateRadius+dUpperBasePlateAnchorLength;
	const G4double dUpperBasePlateOffsetZ = dBottomPlateOffsetZ-dBottomPlateHalfZ-dUpperBasePlateHalfZ;

	const G4double dUpperBasePlateCutInnerRadius = dUpperBasePlateRadius;
	const G4double dUpperBasePlateCutOuterRadius = dUpperBasePlateRadius+2*dUpperBasePlateAnchorLength;
	const G4double dUpperBasePlateCutAngularSize = (-dUpperBasePlateAnchorWidth/dUpperBasePlateRadius*180./M_PI+360./16)*deg;
	const G4double dUpperBasePlateCutHalfZ = dUpperBasePlateThickness;

	G4Tubs *pUpperBasePlateTubs = new G4Tubs("UpperBasePlateTube",
		dUpperBasePlateInnerRadius, dUpperBasePlateOuterRadius, dUpperBasePlateHalfZ, 0.*deg, 360.*deg);

	G4Tubs *pUpperBasePlateSideCut = new G4Tubs("UpperBasePlateSideCut",
		dUpperBasePlateCutInnerRadius, dUpperBasePlateCutOuterRadius, dUpperBasePlateCutHalfZ,
		0, dUpperBasePlateCutAngularSize);

	G4SubtractionSolid *pUpperBasePlateWithCuts, *pUpperBasePlateWithPreviousCuts = (G4SubtractionSolid *) pUpperBasePlateTubs;
	for(int iCutId = 0; iCutId < 16; iCutId++)
	{
		stringstream hStr;
		hStr << "UpperBasePlateTubeWithCut" << iCutId+1;

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ(iCutId*360.*deg/16+0.5*dUpperBasePlateCutAngularSize);

		pUpperBasePlateWithCuts = new G4SubtractionSolid(hStr.str(),
			pUpperBasePlateWithPreviousCuts, pUpperBasePlateSideCut,
			pRotationMatrix, G4ThreeVector(0., 0., 0.));

		pUpperBasePlateWithPreviousCuts = pUpperBasePlateWithCuts;
	}
	
	m_pUpperBasePlateLogicalVolume = new G4LogicalVolume(pUpperBasePlateWithCuts,
		Copper, "UpperBasePlateLogicalVolume", 0, 0, 0);

	m_pUpperBasePlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dUpperBasePlateOffsetZ),
		m_pUpperBasePlateLogicalVolume, "UpperBasePlate", m_pLXeLogicalVolume, false, 0);

	//============================== lower teflon panels ============================
	const G4double dLowerTeflonPanelHalfZ = 0.5*dLowerTeflonPanelHeight;
	const G4double dLowerTeflonPanelOffsetZ = dUpperBasePlateOffsetZ-dUpperBasePlateHalfZ-dLowerTeflonPanelHalfZ;

	const G4double pdLowerTeflonPanelZPlanes[] = {-dLowerTeflonPanelHalfZ, dLowerTeflonPanelHalfZ};
	G4double pdLowerTeflonPanelInnerRadii[32], pdLowerTeflonPanelOuterRadii[32];
	
	for(int i = 0; i < 32; i++)
	{
		pdLowerTeflonPanelInnerRadii[i] = dLowerTeflonPanelRadius;
		pdLowerTeflonPanelOuterRadii[i] = dLowerTeflonPanelRadius+dLowerTeflonPanelThickness;
	}

	G4Polyhedra *pLowerTeflonPanelPolyhedra = new G4Polyhedra("LowerTeflonPanelPolyhedra", 0.*deg, 360.*deg,
		32, 2, pdLowerTeflonPanelZPlanes, pdLowerTeflonPanelInnerRadii, pdLowerTeflonPanelOuterRadii);

	m_pLowerTeflonPanelLogicalVolume = new G4LogicalVolume(pLowerTeflonPanelPolyhedra, Teflon,
			"LowerTeflonPanelLogicalVolume", 0, 0, 0);

	G4RotationMatrix *pLowerTelfonPanelRotationMatrix = new G4RotationMatrix();
	pLowerTelfonPanelRotationMatrix->rotateZ(360.*deg/64);

	m_pLowerTeflonPanelPhysicalVolume = new G4PVPlacement(pLowerTelfonPanelRotationMatrix,
		G4ThreeVector(0., 0., dLowerTeflonPanelOffsetZ), m_pLowerTeflonPanelLogicalVolume,
		"LowerTeflonPanel", m_pLXeLogicalVolume, false, 0);

	//=============================== lower base plate ==============================
	const G4double dLowerBasePlateHalfZ = 0.5*dLowerBasePlateThickness;
	const G4double dLowerBasePlateOffsetZ = dLowerTeflonPanelOffsetZ-dLowerTeflonPanelHalfZ-dLowerBasePlateHalfZ;

	const G4double dPmtDistance = dPmtWidth+dPmtSpacing;

	G4Tubs *pLowerBasePlateTubs = new G4Tubs("LowerBasePlateTube", 0., dLowerBasePlateRadius,
		dLowerBasePlateHalfZ, 0.*deg, 360.*deg);

	G4Box *pLowerBasePlateFirstRowCut = new G4Box("LowerBasePlateFirstRowCut",
		2*dPmtDistance, 5*dPmtDistance, dLowerBasePlateThickness);

	G4SubtractionSolid *pLowerBasePlateWithCut1 = new G4SubtractionSolid("LowerBasePlateWithCut1",
		pLowerBasePlateTubs, pLowerBasePlateFirstRowCut, 0, G4ThreeVector(0., 0., 0.));

	G4Box *pLowerBasePlateSecondRowCut = new G4Box("LowerBasePlateSecondRowCut",
		3.5*dPmtDistance, 4*dPmtDistance, dLowerBasePlateThickness);

	G4SubtractionSolid *pLowerBasePlateWithCut2 = new G4SubtractionSolid("LowerBasePlateWithCut2",
		pLowerBasePlateWithCut1, pLowerBasePlateSecondRowCut, 0, G4ThreeVector(0., 0., 0.));

	G4Box *pLowerBasePlateThirdRowCut = new G4Box("LowerBasePlateThirdRowCut",
		4.5*dPmtDistance, 3*dPmtDistance, dLowerBasePlateThickness);

	G4SubtractionSolid *pLowerBasePlateWithCut3 = new G4SubtractionSolid("LowerBasePlateWithCut3",
		pLowerBasePlateWithCut2, pLowerBasePlateThirdRowCut, 0, G4ThreeVector(0., 0., 0.));

	G4Box *pLowerBasePlateFourthRowCut = new G4Box("LowerBasePlateFourthRowCut",
		5*dPmtDistance, 2*dPmtDistance, dLowerBasePlateThickness);

	G4SubtractionSolid *pLowerBasePlateWithCut4 = new G4SubtractionSolid("LowerBasePlateWithCut4",
		pLowerBasePlateWithCut3, pLowerBasePlateFourthRowCut, 0, G4ThreeVector(0., 0., 0.));

	m_pLowerBasePlateLogicalVolume = new G4LogicalVolume(pLowerBasePlateWithCut4, Copper,
		"LowerBasePlateLogicalVolume", 0, 0, 0);

	m_pLowerBasePlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLowerBasePlateOffsetZ),
		m_pLowerBasePlateLogicalVolume, "LowerBasePlate", m_pLXeLogicalVolume, false, 0);

	//================================ pmt base plate ===============================
	const G4double dBottomPmtPlateHalfZ = 0.5*dBottomPmtPlateThickness;
	const G4double dBottomPmtPlateOffsetZ = dLowerBasePlateOffsetZ-dLowerBasePlateHalfZ-dBottomPmtPlateHalfZ;

	G4Tubs *pBottomPmtPlateTubs = new G4Tubs("BottomPmtPlateTube", 0., dBottomPmtPlateRadius,
		dBottomPmtPlateHalfZ, 0.*deg, 360.*deg);

	m_pBottomPmtPlateLogicalVolume = new G4LogicalVolume(pBottomPmtPlateTubs, Copper,
		"BottomPmtPlateLogicalVolume", 0, 0, 0);

	m_pBottomPmtPlatePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomPmtPlateOffsetZ),
		m_pBottomPmtPlateLogicalVolume, "BottomPmtPlate", m_pLXeLogicalVolume, false, 0);

	//============================== bottom veto angles =============================
	const G4double dBottomVetoAngleHalfX = 0.5*dBottomVetoAngleLength;
	const G4double dBottomVetoAngleHalfY = 0.5*dBottomVetoAngleWidth;
	const G4double dBottomVetoAngleHalfZ = 0.5*dBottomVetoAngleHeight;
	const G4double dBottomVetoAngleOffsetX = dLowerBasePlateSupportRingRadius+dBottomVetoAngleHalfX;
	const G4double dBottomVetoAngleOffsetZ = dLowerBasePlateOffsetZ-dLowerBasePlateHalfZ-dBottomVetoAngleHalfZ;

	const G4double dBottomVetoAngleBottomCutHalfX = dBottomVetoAngleLength-dBottomVetoAngleThickness;
	const G4double dBottomVetoAngleBottomCutHalfZ = dBottomVetoAngleHeight-dBottomVetoAngleThickness;

	G4Box *pBottomVetoAngleBox = new G4Box("BottomVetoAngleBox",
		dBottomVetoAngleHalfX, dBottomVetoAngleHalfY, dBottomVetoAngleHalfZ);

	G4Box *pBottomVetoAngleBottomCut = new G4Box("BottomVetoAngleBottomCut",
		dBottomVetoAngleBottomCutHalfX, dBottomVetoAngleWidth, dBottomVetoAngleBottomCutHalfZ);

	G4SubtractionSolid *pBottomVetoAngleBoxWithCut1 = new G4SubtractionSolid("BottomVetoAngleBoxWithCut1",
		pBottomVetoAngleBox, pBottomVetoAngleBottomCut,
		0, G4ThreeVector(dBottomVetoAngleHalfX, 0, -dBottomVetoAngleHalfZ));

	m_pBottomVetoAngleLogicalVolume = new G4LogicalVolume(pBottomVetoAngleBoxWithCut1,
		Copper, "BottomVetoAngleLogicalVolume", 0, 0, 0);

	for(int iBottomVetoAngleId = 0; iBottomVetoAngleId < 16; iBottomVetoAngleId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(2*iBottomVetoAngleId+1, 32, dBottomVetoAngleOffsetX);
		hPos.setZ(dBottomVetoAngleOffsetZ);

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ((2*iBottomVetoAngleId+1)*360.*deg/32);

		m_hBottomVetoAnglePhysicalVolumes.push_back(new G4PVPlacement(pRotationMatrix, hPos,
			m_pBottomVetoAngleLogicalVolume, "BottomVetoAngle", m_pLXeLogicalVolume, false, 0));
	}

	//============================= lower side veto angles ==========================
	const G4double dLowerSideVetoAngleHalfX = 0.5*(dLowerSideVetoAngleShortLength+dLowerSideVetoAngleLongLength-dLowerSideVetoAngleThickness);
	const G4double dLowerSideVetoAngleHalfY = 0.5*dLowerSideVetoAngleWidth;
	const G4double dLowerSideVetoAngleHalfZ = 0.5*dLowerSideVetoAngleHeight;
	const G4double dLowerSideVetoAngleOffsetX = dLowerBasePlateSupportRingRadius+dLowerSideVetoAngleHalfX;
	const G4double dLowerSideVetoAngleOffsetZ = dLowerBasePlateOffsetZ-dLowerBasePlateHalfZ-dLowerSideVetoAngleHalfZ;

	const G4double dLowerSideVetoAngleShortCutHalfX = dLowerSideVetoAngleShortLength-dLowerSideVetoAngleThickness;
	const G4double dLowerSideVetoAngleShortCutHalfZ = dLowerSideVetoAngleHeight-dLowerSideVetoAngleThickness;

	const G4double dLowerSideVetoAngleLongCutHalfX = dLowerSideVetoAngleLongLength-dLowerSideVetoAngleThickness;
	const G4double dLowerSideVetoAngleLongCutHalfZ = dLowerSideVetoAngleShortCutHalfZ;

	G4Box *pLowerSideVetoAngleBox = new G4Box("LowerSideVetoAngleBox",
		dLowerSideVetoAngleHalfX, dLowerSideVetoAngleHalfY, dLowerSideVetoAngleHalfZ);

	G4Box *pLowerSideVetoAngleShortCut = new G4Box("LowerSideVetoAngleShortCut",
		dLowerSideVetoAngleShortCutHalfX, dLowerSideVetoAngleWidth, dLowerSideVetoAngleShortCutHalfZ);

	G4SubtractionSolid *pLowerSideVetoAngleBoxWithCut1 = new G4SubtractionSolid("LowerSideVetoAngleBoxWithCut1",
		pLowerSideVetoAngleBox, pLowerSideVetoAngleShortCut,
		0, G4ThreeVector(dLowerSideVetoAngleHalfX, 0, -dLowerSideVetoAngleHalfZ));

	G4Box *pLowerSideVetoAngleLongCut = new G4Box("LowerSideVetoAngleLongCut",
		dLowerSideVetoAngleLongCutHalfX, dLowerSideVetoAngleWidth, dLowerSideVetoAngleLongCutHalfZ);

	G4SubtractionSolid *pLowerSideVetoAngleBoxWithCut2 = new G4SubtractionSolid("LowerSideVetoAngleBoxWithCut2",
		pLowerSideVetoAngleBoxWithCut1, pLowerSideVetoAngleLongCut,
		0, G4ThreeVector(-dLowerSideVetoAngleHalfX, 0, dLowerSideVetoAngleHalfZ));

	m_pLowerSideVetoAngleLogicalVolume = new G4LogicalVolume(pLowerSideVetoAngleBoxWithCut2,
		Copper, "LowerSideVetoAngleLogicalVolume", 0, 0, 0);

	for(int iLowerSideVetoAngleId = 0; iLowerSideVetoAngleId < 16; iLowerSideVetoAngleId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(2*iLowerSideVetoAngleId, 32, dLowerSideVetoAngleOffsetX);
		hPos.setZ(dLowerSideVetoAngleOffsetZ);

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ((2*iLowerSideVetoAngleId)*360.*deg/32);

		m_hLowerSideVetoAnglePhysicalVolumes.push_back(new G4PVPlacement(pRotationMatrix, hPos,
			m_pLowerSideVetoAngleLogicalVolume, "LowerSideVetoAngle", m_pLXeLogicalVolume, false, 0));
	}

	//================================== veto lining ================================
	const G4double dTeflonSideVetoLiningInsideGXeHalfZ = max(0., dLidInnerWallDomeTopToVetoLiquidLevel-dLidInnerWallDomeHeight);
	const G4double dTeflonSideVetoLiningHalfZ = 0.5*(dLXeHeight-dVesselInnerWallDomeHeight-dLidInnerWallDomeHeight-dTeflonSideVetoLiningInsideGXeHalfZ);
	const G4double dTeflonSideVetoLiningOffsetZ = -dLXeHalfZ+dVesselInnerWallDomeHeight+dTeflonSideVetoLiningHalfZ;

	G4Tubs *pTeflonSideVetoLiningTubs = new G4Tubs("TeflonSideVetoLiningTube",
		dTeflonSideVetoLiningRadius-dTeflonSideVetoLiningThickness,
		dTeflonSideVetoLiningRadius, dTeflonSideVetoLiningHalfZ, 0.*deg, 360.*deg);

	m_pTeflonSideVetoLiningLogicalVolume = new G4LogicalVolume(pTeflonSideVetoLiningTubs, Teflon,
			"TeflonSideVetoLiningLogicalVolume", 0, 0, 0);

	m_pTeflonSideVetoLiningPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTeflonSideVetoLiningOffsetZ),
		m_pTeflonSideVetoLiningLogicalVolume, "TeflonSideVetoLining", m_pLXeLogicalVolume, false, 0);

	//=============================== teflon led block ==============================
	const G4double dTeflonLedBlockHalfX = 0.5*dTeflonLedBlockWidth;
	const G4double dTeflonLedBlockHalfY = 0.5*dTeflonLedBlockHeight;
	const G4double dTeflonLedBlockHalfZ = 0.5*dTeflonLedBlockDepth;

	const G4double dTeflonLedBlockOffsetZ = dLowerBasePlateOffsetZ+dLowerBasePlateHalfZ+dTeflonLedBlockHalfZ;

	G4Box *pTeflonLedBlock = new G4Box("TeflonLedBlock", dTeflonLedBlockHalfX, dTeflonLedBlockHalfY, dTeflonLedBlockHalfZ);

	m_pTeflonLedBlockLogicalVolume = new G4LogicalVolume(pTeflonLedBlock, Teflon, "LowerBasePlateLogicalVolume", 0, 0, 0);

        // two copies of Logical Volume 
        m_pRotationZTeflonLedBlock = new G4RotationMatrix();                                                                    
        m_pRotationZTeflonLedBlock->rotateZ(dTeflonLedBlockRotationZ);

        m_hTeflonLedBlockPhysicalVolumes.push_back
                (new G4PVPlacement(m_pRotationZTeflonLedBlock,G4ThreeVector(dTeflonLedBlockOffsetX,dTeflonLedBlockOffsetY,dTeflonLedBlockOffsetZ),
                                   m_pTeflonLedBlockLogicalVolume, "TeflonLedBlock",m_pLXeLogicalVolume,false,0));
        m_hTeflonLedBlockPhysicalVolumes.push_back
                (new G4PVPlacement(m_pRotationZTeflonLedBlock,G4ThreeVector(-dTeflonLedBlockOffsetX,-dTeflonLedBlockOffsetY,dTeflonLedBlockOffsetZ),
                                   m_pTeflonLedBlockLogicalVolume, "TeflonLedBlock",m_pLXeLogicalVolume,false,1));
        
                                    
       //====================== bottom PTFE reflector and Holders  ===================
        const G4double dBottomReflectorHolderHalfZ = 0.5*dBottomReflectorHolderHeight;
        const G4double dTeflonBottomReflectorTopOffsetZ = dBottomPmtPlateOffsetZ-0.5*dBottomPlateThickness-0.5*dTeflonBottomReflectorThickness;
        const G4double dBottomReflectorHolderOffsetZ = dTeflonBottomReflectorTopOffsetZ-0.5*dTeflonBottomReflectorThickness-dBottomReflectorHolderHalfZ;
	const G4double dTeflonBottomReflectorBottomOffsetZ = dBottomReflectorHolderOffsetZ-dBottomReflectorHolderHalfZ-0.5*dTeflonBottomReflectorThickness;

        G4cout<<"dTeflonBottomReflectorTopOffsetZ: "<<dTeflonBottomReflectorTopOffsetZ<<G4endl;
        G4cout<<"dTeflonBottomReflectorBottomOffsetZ: "<<dTeflonBottomReflectorBottomOffsetZ<<G4endl;

	G4Tubs *pTeflonBottomReflectorTop = new G4Tubs("TeflonBottomReflectorTop", 0, dTeflonBottomReflectorTopRadius, 
                                                       0.5*dTeflonBottomReflectorThickness, 0.*deg, 360.*deg);
	G4Tubs *pTeflonBottomReflectorBottom = new G4Tubs("TeflonBottomReflectorBottom", 0, dTeflonBottomReflectorBottomRadius, 
                                                          0.5*dTeflonBottomReflectorThickness, 0.*deg, 360.*deg);
	G4Tubs *pBottomReflectorHolder = new G4Tubs("BottomReflectorHolder",0., dBottomReflectorHolderRadius, 
                                                          dBottomReflectorHolderHalfZ, 0.*deg, 360.*deg);
        
	m_pTeflonBottomReflectorTopLogicalVolume = new G4LogicalVolume(pTeflonBottomReflectorTop, Teflon, "TeflonBottomReflectorTopLogicalVolume", 0, 0, 0);
	m_pTeflonBottomReflectorBottomLogicalVolume = new G4LogicalVolume(pTeflonBottomReflectorBottom, Teflon, "TeflonBottomReflectorBottomLogicalVolume", 0, 0, 0);
	m_pBottomReflectorHolderLogicalVolume = new G4LogicalVolume(pBottomReflectorHolder, Copper, "BottomReflectorHolderLogicalVolume", 0, 0, 0);

	
	m_pTeflonBottomReflectorTopPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0.,0.,dTeflonBottomReflectorTopOffsetZ), 
                                                                      m_pTeflonBottomReflectorTopLogicalVolume, "TeflonBottomReflectorTop", m_pLXeLogicalVolume, false, 0);
	m_pTeflonBottomReflectorBottomPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0.,0.,dTeflonBottomReflectorBottomOffsetZ), 
                                                                         m_pTeflonBottomReflectorBottomLogicalVolume, "TeflonBottomReflectorBottom", m_pLXeLogicalVolume, false, 0);

	for(int iBottomReflectorHolderId = 0; iBottomReflectorHolderId < 4; iBottomReflectorHolderId++)
	{
		G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(iBottomReflectorHolderId, 4, dBottomReflectorHolderOffsetX);
		hPos.setZ(dBottomReflectorHolderOffsetZ);

		G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();
		pRotationMatrix->rotateZ((iBottomReflectorHolderId)*360.*deg/8);

		m_hBottomReflectorHolderPhysicalVolumes.push_back(new G4PVPlacement(pRotationMatrix, hPos,
			m_pBottomReflectorHolderLogicalVolume, "BottomReflectorHolder", m_pLXeLogicalVolume, false, 0));
	}

	//=============================== optical surfaces ==============================
	G4double dSigmaAlpha = 0.1;
	G4OpticalSurface *pTeflonOpticalSurface = new G4OpticalSurface("TeflonOpticalSurface",
		unified, ground, dielectric_metal, dSigmaAlpha);
	pTeflonOpticalSurface->SetMaterialPropertiesTable(Teflon->GetMaterialPropertiesTable());	

	new G4LogicalBorderSurface("TeflonPanelLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pTeflonPanelPhysicalVolume, pTeflonOpticalSurface);

	for(int iRodId = 0; iRodId < 16; iRodId++)
	{
		new G4LogicalBorderSurface("TeflonRodsLogicalBorderSurface",
			m_pLXePhysicalVolume, m_hTeflonRodsPhysicalVolumes[iRodId], pTeflonOpticalSurface);
	}
	new G4LogicalBorderSurface("TeflonSideVetoLiningLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pTeflonSideVetoLiningPhysicalVolume, pTeflonOpticalSurface);	

	new G4LogicalBorderSurface("LowerTeflonPanelLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pLowerTeflonPanelPhysicalVolume, pTeflonOpticalSurface);	

	//================================== attributes =================================
	//G4Colour hCopperColor(0.835, 0.424, 0.059, 0.50);
	//G4Colour hTeflonColor(1., 1., 1., 0.5);
	//G4Colour hCopperColour(1.0, 1.0, 0.0);		// yellow
	//G4Colour hTeflonColour(1.0, 0.0, 1.0);		// magenta
	G4Colour hCopperColour(0.0, 0.0, 0.0);
	G4Colour hTeflonColour(0.0, 0.0, 0.0);


	G4VisAttributes *pTopPlateVisAtt = new G4VisAttributes(hCopperColour);
	pTopPlateVisAtt->SetVisibility(true);
	m_pTopPlateLogicalVolume->SetVisAttributes(pTopPlateVisAtt);

	G4VisAttributes *pTeflonPanelVisAtt = new G4VisAttributes(hTeflonColour);
	pTeflonPanelVisAtt->SetVisibility(true);
	m_pTeflonPanelLogicalVolume->SetVisAttributes(pTeflonPanelVisAtt);

	G4VisAttributes *pTeflonRodsVisAtt = new G4VisAttributes(hTeflonColour);
	pTeflonRodsVisAtt->SetVisibility(true);
	m_pTeflonRodsLogicalVolume->SetVisAttributes(pTeflonRodsVisAtt);

	G4VisAttributes *pTeflonSideVetoLiningVisAtt = new G4VisAttributes(hTeflonColour);
	pTeflonSideVetoLiningVisAtt->SetVisibility(false);
	m_pTeflonSideVetoLiningLogicalVolume->SetVisAttributes(pTeflonSideVetoLiningVisAtt);

	G4VisAttributes *pBottomPlateVisAtt = new G4VisAttributes(hCopperColour);
	pBottomPlateVisAtt->SetVisibility(true);
	m_pBottomPlateLogicalVolume->SetVisAttributes(pBottomPlateVisAtt);

	G4VisAttributes *pUpperBasePlateVisAtt = new G4VisAttributes(hCopperColour);
	pUpperBasePlateVisAtt->SetVisibility(true);
	m_pUpperBasePlateLogicalVolume->SetVisAttributes(pUpperBasePlateVisAtt);

	G4VisAttributes *pLowerTeflonPanelVisAtt = new G4VisAttributes(hTeflonColour);
	pLowerTeflonPanelVisAtt->SetVisibility(true);
	m_pLowerTeflonPanelLogicalVolume->SetVisAttributes(pLowerTeflonPanelVisAtt);

	G4VisAttributes *pLowerBasePlateVisAtt = new G4VisAttributes(hCopperColour);
	pLowerBasePlateVisAtt->SetVisibility(true);
	m_pLowerBasePlateLogicalVolume->SetVisAttributes(pLowerBasePlateVisAtt);

	G4VisAttributes *pBottomPmtPlateVisAtt = new G4VisAttributes(hCopperColour);
	pBottomPmtPlateVisAtt->SetVisibility(true);
	m_pBottomPmtPlateLogicalVolume->SetVisAttributes(pBottomPmtPlateVisAtt);

	G4VisAttributes *pBottomVetoAngleVisAtt = new G4VisAttributes(hCopperColour);
	pBottomVetoAngleVisAtt->SetVisibility(true);
	m_pBottomVetoAngleLogicalVolume->SetVisAttributes(pBottomVetoAngleVisAtt);

	G4VisAttributes *pLowerSideVetoAngleVisAtt = new G4VisAttributes(hCopperColour);
	pLowerSideVetoAngleVisAtt->SetVisibility(true);
	m_pLowerSideVetoAngleLogicalVolume->SetVisAttributes(pLowerSideVetoAngleVisAtt);

	G4VisAttributes *pTeflonLedBlockVisAtt = new G4VisAttributes(hTeflonColour);
	pTeflonLedBlockVisAtt->SetVisibility(true);
	m_pTeflonLedBlockLogicalVolume->SetVisAttributes(pTeflonLedBlockVisAtt);

	G4VisAttributes *pTeflonBottomReflectorVisAtt = new G4VisAttributes(hTeflonColour);
	pTeflonBottomReflectorVisAtt->SetVisibility(true);
	m_pTeflonBottomReflectorTopLogicalVolume->SetVisAttributes(pTeflonBottomReflectorVisAtt);
	m_pTeflonBottomReflectorBottomLogicalVolume->SetVisAttributes(pTeflonBottomReflectorVisAtt);

	G4VisAttributes *pBottomReflectorHolderVisAtt = new G4VisAttributes(hCopperColour);
	pBottomReflectorHolderVisAtt->SetVisibility(true);
	m_pBottomReflectorHolderLogicalVolume->SetVisAttributes(pBottomReflectorHolderVisAtt);
}

void
Xenon100DetectorConstruction::ConstructPmtArrays()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< pmt arrays >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dPmtWindowWidth = GetGeometryParameter("PmtWindowWidth");
	const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");

	const G4double dPmtCasingWidth = GetGeometryParameter("PmtCasingWidth");
	const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
	const G4double dPmtCasingThickness = GetGeometryParameter("PmtCasingThickness");

	const G4double dPmtPhotoCathodeWidth = GetGeometryParameter("PmtPhotoCathodeWidth");
	const G4double dPmtPhotoCathodeThickness = GetGeometryParameter("PmtPhotoCathodeThickness");

	const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

	G4Material *Quartz = G4Material::GetMaterial("Quartz");
	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
	G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
	G4Material *Aluminium = G4Material::GetMaterial("PhotoCathodeAluminium");
	G4Material *Cirlex = G4Material::GetMaterial("Cirlex");

	//==================================== pmts =====================================

	//--------------------------------- pmt window ----------------------------------
	const G4double dPmtWindowHalfX = 0.5*dPmtWindowWidth;
	const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;

	G4Box *pPmtWindowBox = new G4Box("PmtWindowBox", dPmtWindowHalfX, dPmtWindowHalfX, dPmtWindowHalfZ);

	m_pPmtWindowLogicalVolume = new G4LogicalVolume(pPmtWindowBox, Quartz, "PmtWindowLogicalVolume", 0, 0, 0);

	//--------------------------------- pmt casing ----------------------------------
	const G4double dPmtCasingHalfX = 0.5*dPmtCasingWidth;
	const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;

	G4Box *pPmtCasingBox = new G4Box("PmtCasingBox", dPmtCasingHalfX, dPmtCasingHalfX, dPmtCasingHalfZ);

	m_pPmtCasingLogicalVolume = new G4LogicalVolume(pPmtCasingBox, SS304LSteel, "PmtCasingLogicalVolume", 0, 0, 0);

	//-------------------------------- pmt interior ---------------------------------
	const G4double dPmtInteriorHalfX = dPmtCasingHalfX-dPmtCasingThickness;
	const G4double dPmtInteriorHalfZ = dPmtCasingHalfZ-dPmtCasingThickness;

	G4Box *pPmtInteriorBox = new G4Box("PmtInteriorBox", dPmtInteriorHalfX, dPmtInteriorHalfX, dPmtInteriorHalfZ);

	m_pPmtInteriorLogicalVolume = new G4LogicalVolume(pPmtInteriorBox, Vacuum, "PmtInteriorLogicalVolume", 0, 0, 0);

	m_pPmtInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPmtInteriorLogicalVolume, "PmtInterior", m_pPmtCasingLogicalVolume, false, 0);

	//------------------------------ pmt photocathode -------------------------------
	const G4double dPmtPhotoCathodeHalfX = 0.5*dPmtPhotoCathodeWidth;
	const G4double dPmtPhotoCathodeHalfZ = 0.5*dPmtPhotoCathodeThickness;
	const G4double dPmtPhotoCathodeOffsetZ = dPmtCasingHalfZ-dPmtPhotoCathodeHalfZ;

	G4Box *pPmtPhotoCathodeBox = new G4Box("PmtPhotoCathodeBox", dPmtPhotoCathodeHalfX, dPmtPhotoCathodeHalfX, dPmtPhotoCathodeHalfZ);

	m_pPmtPhotoCathodeLogicalVolume = new G4LogicalVolume(pPmtPhotoCathodeBox, Aluminium, "PmtPhotoCathodeLogicalVolume", 0, 0, 0);

	m_pPmtPhotoCathodePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., -dPmtPhotoCathodeOffsetZ),
		m_pPmtPhotoCathodeLogicalVolume, "PmtPhotoCathode", m_pPmtCasingLogicalVolume, false, 0);

	//---------------------------------- pmt base -----------------------------------
	const G4double dPmtBaseHalfX = dPmtCasingHalfX;
	const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

	G4Box *pPmtBaseBox = new G4Box("PmtBaseBox", dPmtBaseHalfX, dPmtBaseHalfX, dPmtBaseHalfZ);

	m_pPmtBaseLogicalVolume = new G4LogicalVolume(pPmtBaseBox, Cirlex, "PmtBaseLogicalVolume", 0, 0, 0);

	//================================== top array ==================================
	G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");

	stringstream hVolumeName;
	for(G4int iPmtNb=0; iPmtNb<iNbTopPmts; iPmtNb++)
	{
		hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;

		m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
			GetPmtPosition(iPmtNb, PMT_WINDOW), m_pPmtWindowLogicalVolume,
			hVolumeName.str(), m_pGXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

		m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
			GetPmtPosition(iPmtNb, PMT_CASING), m_pPmtCasingLogicalVolume,
			hVolumeName.str(), m_pGXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

		m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
			GetPmtPosition(iPmtNb, PMT_BASE), m_pPmtBaseLogicalVolume, hVolumeName.str(),
			m_pGXeLogicalVolume, false, iPmtNb));
	}
	
	//================================ bottom array =================================
	G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");

	for(G4int iPmtNb=iNbTopPmts; iPmtNb<iNbTopPmts+iNbBottomPmts; iPmtNb++)
	{
		hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;

		m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
			GetPmtPosition(iPmtNb, PMT_WINDOW), m_pPmtWindowLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

		m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
			GetPmtPosition(iPmtNb, PMT_CASING), m_pPmtCasingLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

		m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
			GetPmtPosition(iPmtNb, PMT_BASE), m_pPmtBaseLogicalVolume, hVolumeName.str(),
			m_pLXeLogicalVolume, false, iPmtNb));
	}

	//=============================== top veto array ================================
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");

	for(G4int iPmtNb=iNbTopPmts+iNbBottomPmts; iPmtNb<iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts; iPmtNb++)
	{
		hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;

		m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
			GetPmtPosition(iPmtNb, PMT_WINDOW), m_pPmtWindowLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

		m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
			GetPmtPosition(iPmtNb, PMT_CASING), m_pPmtCasingLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

		m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
			GetPmtPosition(iPmtNb, PMT_BASE), m_pPmtBaseLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));
	}

	//============================== bottom veto array ==============================
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");

	for(G4int iPmtNb=iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts; iPmtNb<iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts; iPmtNb++)
	{
		hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;

		m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
			GetPmtPosition(iPmtNb, PMT_WINDOW), m_pPmtWindowLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

		m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
			GetPmtPosition(iPmtNb, PMT_CASING), m_pPmtCasingLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));

		hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

		m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
			GetPmtPosition(iPmtNb, PMT_BASE), m_pPmtBaseLogicalVolume,
			hVolumeName.str(), m_pLXeLogicalVolume, false, iPmtNb));
	}

	//------------------------------- pmt sensitivity -------------------------------
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	Xenon100PmtSensitiveDetector *pPmtSD = new Xenon100PmtSensitiveDetector("Xenon100/PmtSD");
	pSDManager->AddNewDetector(pPmtSD);
	m_pPmtPhotoCathodeLogicalVolume->SetSensitiveDetector(pPmtSD);

	//---------------------------------- attributes ---------------------------------
	//G4Colour hPMTColour(1.0, 0.0, 0.0);			// red
	//G4Colour hQuartzColour(0.0, 0.5, 0.0);		// xlgreen
	//G4Colour hGridColour(0.5, 0.5, 0.5);		// grey
	//G4Colour hPhotoCathodeColour(1.0, 0.5, 0.0);// orange
	//G4Colour hCirlexColour(0.7, 0.4, 0.1);		// brown

	G4Colour hPMTColour(0.0, 0.0, 0.0);			// red
	G4Colour hQuartzColour(0.0, 0.0, 0.0);		// xlgreen
	G4Colour hGridColour(0.0, 0.0, 0.0);		// grey
	G4Colour hPhotoCathodeColour(0.0, 0.0, 0.0);// orange
	G4Colour hCirlexColour(0.0, 0.0, 0.0);		// brown
	
	//G4Colour hPmtWindowColor(1., 0.757, 0.024);
	G4VisAttributes *pPmtWindowVisAtt = new G4VisAttributes(hQuartzColour);
	pPmtWindowVisAtt->SetVisibility(true);
	m_pPmtWindowLogicalVolume->SetVisAttributes(pPmtWindowVisAtt);

	//G4Colour hPmtCasingColor(1., 0.486, 0.027);
	G4VisAttributes *pPmtCasingVisAtt = new G4VisAttributes(hPMTColour);
	pPmtCasingVisAtt->SetVisibility(true);
	m_pPmtCasingLogicalVolume->SetVisAttributes(pPmtCasingVisAtt);

	G4VisAttributes *pPmtInteriorVisAtt = new G4VisAttributes();
	pPmtInteriorVisAtt->SetVisibility(false);
	m_pPmtInteriorLogicalVolume->SetVisAttributes(pPmtInteriorVisAtt);	

	//G4Colour hPmtPhotoCathodeColor(1., 0.082, 0.011);
	G4VisAttributes *pPmtPhotoCathodeVisAtt = new G4VisAttributes(hPhotoCathodeColour);
	pPmtPhotoCathodeVisAtt->SetVisibility(true);
	m_pPmtPhotoCathodeLogicalVolume->SetVisAttributes(pPmtPhotoCathodeVisAtt);

	//G4Colour hPmtBaseColor(0.788, 0.188, 0.780);
	G4VisAttributes *pPmtBaseVisAtt = new G4VisAttributes(hCirlexColour);
	pPmtBaseVisAtt->SetVisibility(true);
	m_pPmtBaseLogicalVolume->SetVisAttributes(pPmtBaseVisAtt);
}

void
Xenon100DetectorConstruction::ConstructCryostat()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cryostat >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

	const G4double dVesselFlangeToZero = GetGeometryParameter("CryostatVesselFlangeToZero");

	const G4double dVesselFlangeRadius = GetGeometryParameter("CryostatVesselFlangeRadius");
	const G4double dVesselFlangeThickness = GetGeometryParameter("CryostatVesselFlangeThickness");

	const G4double dVesselInnerWallRadius = GetGeometryParameter("CryostatVesselInnerWallRadius");
	const G4double dVesselInnerWallHeight = GetGeometryParameter("CryostatVesselInnerWallHeight");
	const G4double dVesselInnerWallThickness = GetGeometryParameter("CryostatVesselInnerWallThickness");
	const G4double dVesselInnerWallDomeHeight = GetGeometryParameter("CryostatVesselInnerWallDomeHeight");
	const G4double dVesselInnerWallDomeRadius = GetGeometryParameter("CryostatVesselInnerWallDomeRadius");

	const G4double dVesselOuterWallRadius = GetGeometryParameter("CryostatVesselOuterWallRadius");
	const G4double dVesselOuterWallThickness = GetGeometryParameter("CryostatVesselOuterWallThickness");

	const G4double dLidFlangeRadius = GetGeometryParameter("CryostatLidFlangeRadius");
	const G4double dLidFlangeThickness = GetGeometryParameter("CryostatLidFlangeThickness");

	const G4double dLidInnerWallRadius = GetGeometryParameter("CryostatLidInnerWallRadius");
	const G4double dLidInnerWallHeight = GetGeometryParameter("CryostatLidInnerWallHeight");
	const G4double dLidInnerWallThickness = GetGeometryParameter("CryostatLidInnerWallThickness");
	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInnerWallDomeRadius = GetGeometryParameter("CryostatLidInnerWallDomeRadius");

	const G4double dLidInnerWallDomeTopToVetoLiquidLevel = GetGeometryParameter("CryostatLidInnerWallDomeTopToVetoLiquidLevel");

	const G4double dLidInTubeRadius = GetGeometryParameter("CryostatLidInTubeRadius");
	const G4double dLidInTubeThickness = GetGeometryParameter("CryostatLidInTubeThickness");
	const G4double dLidInTubeHeight = GetGeometryParameter("CryostatLidInTubeHeight");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");
	const G4double dLidInTubeElbowLength = GetGeometryParameter("CryostatLidInTubeElbowLength");
	const G4double dLidInTubeTopLength = GetGeometryParameter("CryostatLidInTubeTopLength");
	const G4double dLidInTubeFlangeRadius = GetGeometryParameter("CryostatLidInTubeFlangeRadius");
	const G4double dLidInTubeFlangeThickness = GetGeometryParameter("CryostatLidInTubeFlangeThickness");
	const G4double dLidInTubeTopFlangeRadius = GetGeometryParameter("CryostatLidInTubeTopFlangeRadius");
	const G4double dLidInTubeTopFlangeThickness = GetGeometryParameter("CryostatLidInTubeTopFlangeThickness");
	const G4double dLidInTubeJacketRadius = GetGeometryParameter("CryostatLidInTubeJacketRadius");
	const G4double dLidInTubeJacketHeight = GetGeometryParameter("CryostatLidInTubeJacketHeight");
	const G4double dLidInTubeOuterJacketRadius = GetGeometryParameter("CryostatLidInTubeOuterJacketRadius");
	const G4double dLidInTubeOuterJacketHeight = GetGeometryParameter("CryostatLidInTubeOuterJacketHeight");
	const G4double dLidInTubeOffsetX = GetGeometryParameter("CryostatLidInTubeOffsetX");

	const G4double dLidOutTubeRadius = GetGeometryParameter("CryostatLidOutTubeRadius");
	const G4double dLidOutTubeThickness = GetGeometryParameter("CryostatLidOutTubeThickness");
	const G4double dLidOutTubeHeight = GetGeometryParameter("CryostatLidOutTubeHeight");
	const G4double dLidOutTubeElbowLength = GetGeometryParameter("CryostatLidOutTubeElbowLength");
	const G4double dLidOutTubeTopLength = GetGeometryParameter("CryostatLidOutTubeTopLength");
	const G4double dLidOutTubeJacketRadius = GetGeometryParameter("CryostatLidOutTubeJacketRadius");
	const G4double dLidOutTubeOffsetX = GetGeometryParameter("CryostatLidOutTubeOffsetX");

	const G4double dLidOutTubeTopFlangeRadius = GetGeometryParameter("CryostatLidOutTubeTopFlangeRadius");
	const G4double dLidOutTubeTopFlangeThickness = GetGeometryParameter("CryostatLidOutTubeTopFlangeThickness");

	const G4double dLidOutTubeJacketHeight = GetGeometryParameter("CryostatLidOutTubeJacketHeight");

	const G4double dLidOuterWallRadius = GetGeometryParameter("CryostatLidOuterWallRadius");
	const G4double dLidOuterWallThickness = GetGeometryParameter("CryostatLidOuterWallThickness");

	const G4double dZeroOffset = (m_pMotherLogicalVolume == m_pShieldCavityLogicalVolume)?(GetGeometryParameter("ShieldZeroOffset")):(0.);

	const G4double dShieldHalfZ = 0.5*GetGeometryParameter("ShieldHeight");
	const G4double dShieldHalfX = 0.5*GetGeometryParameter("ShieldWidth");
	const G4double dShieldHalfY = 0.5*GetGeometryParameter("ShieldDepth");

	const G4double dPolishLeadThickness = GetGeometryParameter("ShieldPolishLeadThickness");
	const G4double dFrenchLeadThickness = GetGeometryParameter("ShieldFrenchLeadThickness");
	const G4double dPolyethyleneThickness = GetGeometryParameter("ShieldPolyethyleneThickness");
	const G4double dCopperTopThickness = GetGeometryParameter("ShieldCopperTopThickness");
	const G4double dCopperSideThickness = GetGeometryParameter("ShieldCopperThicknessThick");
	const G4double dCopperBottomThickness = GetGeometryParameter("ShieldCopperBottomThickness");

	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
	G4Material *Vacuum = G4Material::GetMaterial("Vacuum");

	//=============================== cryostat vessel ===============================

	//-------------------------------- vessel flange --------------------------------
	const G4double dVesselFlangeHalfZ = 0.5*dVesselFlangeThickness;
	const G4double dVesselFlangeOffsetZ = dVesselFlangeToZero-dVesselFlangeHalfZ;

	G4Tubs *pCryostatVesselFlangeTubs = new G4Tubs("CryostatVesselFlangeTube", dVesselInnerWallRadius,
		dVesselFlangeRadius, dVesselFlangeHalfZ, 0.*deg, 360.*deg);

	m_pCryostatVesselFlangeLogicalVolume = new G4LogicalVolume(pCryostatVesselFlangeTubs,
		SS304LSteel, "CryostatVesselFlangeLogicalVolume", 0, 0, 0);

	m_pCryostatVesselFlangePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dVesselFlangeOffsetZ),
		m_pCryostatVesselFlangeLogicalVolume, "CryostatVesselFlange", m_pMotherLogicalVolume, false, 0);

	//------------------------------ vessel inner wall ------------------------------
	const G4double dVesselInnerWallHalfZ = 0.5*dVesselInnerWallHeight;
	const G4double dVesselInnerWallOffsetZ = dVesselFlangeOffsetZ-dVesselFlangeHalfZ-dVesselInnerWallHalfZ;

	G4Tubs *pCryostatVesselInnerWallTubs = new G4Tubs("CryostatVesselInnerWallTube", dVesselInnerWallRadius,
		dVesselInnerWallRadius+dVesselInnerWallThickness, dVesselInnerWallHalfZ, 0.*deg, 360.*deg);

	m_pCryostatVesselInnerWallLogicalVolume = new G4LogicalVolume(pCryostatVesselInnerWallTubs,
		SS304LSteel, "CryostatVesselInnerWallLogicalVolume", 0, 0, 0);

	m_pCryostatVesselInnerWallPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dVesselInnerWallOffsetZ),
		m_pCryostatVesselInnerWallLogicalVolume, "CryostatVesselInnerWall", m_pMotherLogicalVolume, false, 0);

	//--------------------------- vessel inner wall dome ----------------------------
	const G4double dVesselInnerWallDomeThickness = sqrt(pow(dVesselInnerWallDomeRadius-dVesselInnerWallDomeHeight,2)+pow(dVesselInnerWallRadius+dVesselInnerWallThickness,2))-dVesselInnerWallDomeRadius;
	const G4double dVesselInnerWallDomeOffsetZ = dVesselFlangeOffsetZ-dVesselFlangeHalfZ-dVesselInnerWallHeight+dVesselInnerWallDomeRadius-dVesselInnerWallDomeHeight;

	G4Sphere *pVesselInnerWallDomeSphere = new G4Sphere("VesselInnerWallDomeSphere",
		dVesselInnerWallDomeRadius, dVesselInnerWallDomeRadius+dVesselInnerWallDomeThickness,
		0.*deg, 360.*deg, 90.*deg, 270.*deg);

	G4Tubs *pVesselInnerWallDomeSphereCut = new G4Tubs("VesselInnerWallDomeSphereCut",
		0., dVesselInnerWallRadius+dVesselInnerWallThickness, dVesselInnerWallDomeHeight+dVesselInnerWallDomeThickness, 0.*deg, 360.*deg);

	G4IntersectionSolid *pCryostatVesselInnerWallDome = new G4IntersectionSolid("CryostatVesselInnerWallDome",
		pVesselInnerWallDomeSphere, pVesselInnerWallDomeSphereCut, 0,
		G4ThreeVector(0., 0., -dVesselInnerWallDomeRadius-dVesselInnerWallDomeThickness));

	m_pCryostatVesselInnerWallDomeLogicalVolume = new G4LogicalVolume(pCryostatVesselInnerWallDome,
		SS304LSteel, "CryostatVesselInnerWallDomeLogicalVolume", 0, 0, 0);

	m_pCryostatVesselInnerWallDomePhysicalVolume = new G4PVPlacement(0,
		G4ThreeVector(0., 0., dVesselInnerWallDomeOffsetZ), m_pCryostatVesselInnerWallDomeLogicalVolume,
		"CryostatVesselInnerWallDome", m_pMotherLogicalVolume, false, 0);

	//--------------------------- vessel outer wall dome ----------------------------
	const G4double dVesselOuterWallDomeRadius = dVesselInnerWallDomeRadius+dVesselInnerWallDomeThickness+25.4*mm;
	const G4double dVesselOuterWallHeightDifference = sqrt(pow(dVesselOuterWallDomeRadius,2)-pow(dVesselOuterWallRadius,2))-(dVesselInnerWallDomeRadius-dVesselInnerWallDomeHeight);
	const G4double dVesselOuterWallHeight = dVesselInnerWallHeight+dVesselOuterWallHeightDifference;
	const G4double dVesselOuterWallDomeHeight = dVesselInnerWallDomeHeight+dVesselInnerWallDomeThickness+25.4*mm-dVesselOuterWallHeightDifference;
	const G4double dVesselOuterWallDomeThickness = sqrt(pow(dVesselOuterWallDomeRadius-dVesselOuterWallDomeHeight,2)+pow(dVesselOuterWallRadius+dVesselOuterWallThickness,2))-dVesselOuterWallDomeRadius;
	const G4double dVesselOuterWallDomeOffsetZ = dVesselInnerWallDomeOffsetZ;

	G4Sphere *pOuterWallDomeSphere = new G4Sphere("OuterWallDomeSphere",
		dVesselOuterWallDomeRadius, dVesselOuterWallDomeRadius+dVesselOuterWallDomeThickness,
		0.*deg, 360.*deg, 90.*deg, 270.*deg);

	G4Tubs *pOuterWallDomeSphereCut = new G4Tubs("OuterWallDomeSphereCut",
		0., dVesselOuterWallRadius+dVesselOuterWallThickness, dVesselOuterWallDomeHeight+dVesselOuterWallDomeThickness, 0.*deg, 360.*deg);

	G4IntersectionSolid *pCryostatVesselOuterWallDome = new G4IntersectionSolid("CryostatVesselOuterWallDome",
		pOuterWallDomeSphere, pOuterWallDomeSphereCut, 0,
		G4ThreeVector(0., 0., -dVesselOuterWallDomeRadius-dVesselOuterWallDomeThickness));

	m_pCryostatVesselOuterWallDomeLogicalVolume = new G4LogicalVolume(pCryostatVesselOuterWallDome,
		SS304LSteel, "CryostatVesselOuterWallDomeLogicalVolume", 0, 0, 0);

	m_pCryostatVesselOuterWallDomePhysicalVolume = new G4PVPlacement(0,
		G4ThreeVector(0., 0., dVesselOuterWallDomeOffsetZ), m_pCryostatVesselOuterWallDomeLogicalVolume,
		"CryostatVesselOuterWallDome", m_pMotherLogicalVolume, false, 0);

	//------------------------------ vessel outer wall ------------------------------
	const G4double dVesselOuterWallHalfZ = 0.5*dVesselOuterWallHeight;
	const G4double dVesselOuterWallOffsetZ = dVesselFlangeOffsetZ-dVesselFlangeHalfZ-dVesselOuterWallHalfZ;

	G4Tubs *pCryostatVesselOuterWallTubs = new G4Tubs("CryostatVesselOuterWallTube", dVesselOuterWallRadius,
		dVesselOuterWallRadius+dVesselOuterWallThickness, dVesselOuterWallHalfZ, 0.*deg, 360.*deg);

	m_pCryostatVesselOuterWallLogicalVolume = new G4LogicalVolume(pCryostatVesselOuterWallTubs,
		SS304LSteel, "CryostatVesselOuterWallLogicalVolume", 0, 0, 0);

	m_pCryostatVesselOuterWallPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dVesselOuterWallOffsetZ),
		m_pCryostatVesselOuterWallLogicalVolume, "CryostatVesselOuterWall", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- vessel vacuum --------------------------------
	const G4double dVacuumVesselInnerRadius = dVesselInnerWallRadius+dVesselInnerWallThickness;
	const G4double dVacuumVesselOuterRadius = dVesselOuterWallRadius;
	const G4double dVacuumVesselHalfZ = dVesselOuterWallHalfZ;

	const G4double dVacuumVesselDomeRadius = dVesselInnerWallDomeRadius+dVesselInnerWallDomeThickness;
	const G4double dVacuumVesselDomeTheta = M_PI-asin(dVesselOuterWallRadius/dVesselOuterWallDomeRadius);
	const G4double dVacuumVesselDomeOffsetZ = -dVesselOuterWallHalfZ+dVesselOuterWallDomeRadius-dVesselOuterWallDomeHeight;

	const G4double dVacuumVesselOffsetZ = dVesselOuterWallOffsetZ;
	
	G4Tubs *pVacuumCryostatVesselTubs = new G4Tubs("VacuumCryostatVesselTubs",
		dVacuumVesselInnerRadius, dVacuumVesselOuterRadius, dVacuumVesselHalfZ, 0.*deg, 360.*deg);

	G4Sphere *pVacuumCryostatVesselDomeSphere = new G4Sphere("VacuumCryostatVesselDomeSphere",
		dVacuumVesselDomeRadius, dVacuumVesselDomeRadius+25.4*mm,
		0.*deg, 360.*deg, dVacuumVesselDomeTheta*rad, M_PI*rad);

	G4UnionSolid *pVacuumCryostatVesselSolid = new G4UnionSolid("VacuumCryostatVesselSolid",
		pVacuumCryostatVesselTubs, pVacuumCryostatVesselDomeSphere, 0,
		G4ThreeVector(0., 0., dVacuumVesselDomeOffsetZ));

	m_pVacuumCryostatVesselLogicalVolume = new G4LogicalVolume(pVacuumCryostatVesselSolid,
		Vacuum, "VacuumCryostatVesselLogicalVolume", 0, 0, 0);

	m_pVacuumCryostatVesselPhysicalVolume = new G4PVPlacement(0,
		G4ThreeVector(0., 0., dVacuumVesselOffsetZ),
		m_pVacuumCryostatVesselLogicalVolume, "VacuumCryostatVessel", m_pMotherLogicalVolume, false, 0);

	//================================= cryostat lid ================================

	//---------------------------------- lid flange ---------------------------------
	const G4double dLidFlangeInnerRadius = dLidOuterWallRadius-dLidOuterWallThickness;
	const G4double dLidFlangeHalfZ = 0.5*dLidFlangeThickness;
	const G4double dLidFlangeOffsetZ = dVesselFlangeToZero+dLidFlangeHalfZ;

	G4Tubs *pCryostatLidFlangeTubs = new G4Tubs("CryostatLidFlangeTube", dLidFlangeInnerRadius,
		dLidFlangeRadius, dLidFlangeHalfZ, 0.*deg, 360.*deg);

	m_pCryostatLidFlangeLogicalVolume = new G4LogicalVolume(pCryostatLidFlangeTubs,
		SS304LSteel, "CryostatLidFlangeLogicalVolume", 0, 0, 0);

	m_pCryostatLidFlangePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLidFlangeOffsetZ),
		m_pCryostatLidFlangeLogicalVolume, "CryostatLidFlange", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- lid inner wall -------------------------------
	const G4double dLidInnerWallInnerRadius = dLidInnerWallRadius-dLidInnerWallThickness;
	const G4double dLidInnerWallHalfZ = 0.5*dLidInnerWallHeight;
	const G4double dLidInnerWallOffsetZ = dLidFlangeOffsetZ-dLidFlangeHalfZ-dLidInnerWallHalfZ;

	const G4double dLidInnerWallDomeThickness = dLidInnerWallThickness;

	G4Tubs *pCryostatLidInnerWallTubs = new G4Tubs("CryostatLidInnerWallTube", dLidInnerWallInnerRadius,
		dLidInnerWallRadius, dLidInnerWallHalfZ, 0.*deg, 360.*deg);

	G4Orb *pCryostatLidInnerWallBottomCut = new G4Orb("CryostatLidInnerWallBottomCut", dLidInnerWallDomeRadius+dLidInnerWallDomeThickness);

	G4SubtractionSolid *pCryostatLidInnerWallWithCut1 = new G4SubtractionSolid("CryostatLidInnerWallWithCut1",
		pCryostatLidInnerWallTubs, pCryostatLidInnerWallBottomCut, 0,
		G4ThreeVector(0., 0., -dLidInnerWallHalfZ-dLidInnerWallDomeRadius+dLidInnerWallDomeHeight));

	m_pCryostatLidInnerWallLogicalVolume = new G4LogicalVolume(pCryostatLidInnerWallWithCut1,
		SS304LSteel, "CryostatLidInnerWallLogicalVolume", 0, 0, 0);

	m_pCryostatLidInnerWallPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLidInnerWallOffsetZ),
		m_pCryostatLidInnerWallLogicalVolume, "CryostatLidInnerWall", m_pMotherLogicalVolume, false, 0);

	//----------------------------- lid inner wall dome -----------------------------
	const G4double dLidInnerWallDomeOffsetZ = dLidInnerWallOffsetZ-dLidInnerWallHalfZ-dLidInnerWallDomeRadius+dLidInnerWallDomeHeight;

	const G4double dLidInnerWallDomeInTubeCutRadius = dLidInTubeJacketRadius;
	const G4double dLidInnerWallDomeInTubeCutOffsetX = -dLidInTubeOffsetX;
	const G4double dLidInnerWallDomeOutTubeCutRadius = dLidOutTubeRadius;
	const G4double dLidInnerWallDomeOutTubeCutOffsetX = dLidOutTubeOffsetX;

	G4Sphere *pCryostatLidInnerWallDomeSphere = new G4Sphere("CryostatLidInnerWallDomeSphere",
		dLidInnerWallDomeRadius, dLidInnerWallDomeRadius+dLidInnerWallDomeThickness,
		0.*deg, 360.*deg, 0.*deg, 90.*deg);

	G4Tubs *pCryostatLidInnerWallDomeSphereCut = new G4Tubs("CryostatLidInnerWallDomeSphereCut",
		0., dLidInnerWallRadius, 2*dLidInnerWallDomeHeight, 0.*deg, 360.*deg);

	G4IntersectionSolid *pCryostatLidInnerWallDome = new G4IntersectionSolid("CryostatLidInnerWallDome",
		pCryostatLidInnerWallDomeSphere, pCryostatLidInnerWallDomeSphereCut, 0,
		G4ThreeVector(0., 0., dLidInnerWallDomeRadius+dLidInnerWallDomeThickness));

	G4Tubs *pCryostatLidInnerWallInTubeCut = new G4Tubs("CryostatLidInnerWallInTubeCut",
		0., dLidInnerWallDomeInTubeCutRadius, dLidInnerWallDomeRadius, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidInnerWallDomeWithCut1 = new G4SubtractionSolid(
		"CryostatLidInnerWallWithCut1", pCryostatLidInnerWallDome, pCryostatLidInnerWallInTubeCut,
		0, G4ThreeVector(dLidInnerWallDomeInTubeCutOffsetX, 0., 0.));

	G4Tubs *pCryostatLidInnerWallOutTubeCut = new G4Tubs("CryostatLidInnerWallOutTubeCut",
		0., dLidInnerWallDomeOutTubeCutRadius, dLidInnerWallDomeRadius, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidInnerWallDomeWithCut2 = new G4SubtractionSolid(
		"CryostatLidInnerWallWithCut2", pCryostatLidInnerWallDomeWithCut1, pCryostatLidInnerWallOutTubeCut,
		0, G4ThreeVector(dLidInnerWallDomeOutTubeCutOffsetX, 0., 0.));

	m_pCryostatLidInnerWallDomeLogicalVolume = new G4LogicalVolume(pCryostatLidInnerWallDomeWithCut2,
		SS304LSteel, "CryostatLidInnerWallDomeLogicalVolume", 0, 0, 0);

	m_pCryostatLidInnerWallDomePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLidInnerWallDomeOffsetZ),
		m_pCryostatLidInnerWallDomeLogicalVolume, "CryostatLidInnerWallDome", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- lid outer wall -------------------------------
	const G4double dLidOuterWallDomeRadius = dLidInnerWallDomeRadius+dLidInnerWallDomeThickness+25.4*mm;
	const G4double dLidOuterWallDomeHeight = dLidOuterWallDomeRadius-sqrt(pow(dLidOuterWallDomeRadius,2)-pow(dLidOuterWallRadius,2));
	const G4double dLidOuterWallHeightDifference = dLidOuterWallDomeRadius-dLidOuterWallDomeHeight-(dLidInnerWallDomeRadius-dLidInnerWallDomeHeight);
	const G4double dLidOuterWallHeight = dLidInnerWallHeight-dLidOuterWallHeightDifference;	
	const G4double dLidOuterWallDomeThickness = dLidOuterWallThickness;

	const G4double dLidOuterWallInnerRadius = dLidOuterWallRadius-dLidOuterWallThickness;
	const G4double dLidOuterWallHalfZ = 0.5*dLidOuterWallHeight;
	const G4double dLidOuterWallOffsetZ = dLidFlangeOffsetZ-dLidFlangeHalfZ-dLidOuterWallHalfZ;

	G4Tubs *pCryostatLidOuterWallTubs = new G4Tubs("CryostatLidOuterWallTube", dLidOuterWallInnerRadius,
		dLidOuterWallRadius, dLidOuterWallHalfZ, 0.*deg, 360.*deg);

	G4Orb *pCryostatLidOuterWallBottomCut = new G4Orb("CryostatLidOuterWallBottomCut", dLidOuterWallDomeRadius+dLidOuterWallDomeThickness);

	G4SubtractionSolid *pCryostatLidOuterWallWithCut1 = new G4SubtractionSolid("CryostatLidOuterWallWithCut1",
		pCryostatLidOuterWallTubs, pCryostatLidOuterWallBottomCut, 0,
		G4ThreeVector(0., 0., -0.5*dLidOuterWallHeight-dLidOuterWallDomeRadius+dLidOuterWallDomeHeight));

	m_pCryostatLidOuterWallLogicalVolume = new G4LogicalVolume(pCryostatLidOuterWallWithCut1,
		SS304LSteel, "CryostatLidOuterWallLogicalVolume", 0, 0, 0);

	m_pCryostatLidOuterWallPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLidOuterWallOffsetZ),
		m_pCryostatLidOuterWallLogicalVolume, "CryostatLidOuterWall", m_pMotherLogicalVolume, false, 0);

	//----------------------------- lid outer wall dome -----------------------------
	const G4double dLidOuterWallDomeOffsetZ = dLidOuterWallOffsetZ-0.5*dLidOuterWallHeight-dLidOuterWallDomeRadius+dLidOuterWallDomeHeight;

	const G4double dLidOuterWallDomeInTubeCutRadius = dLidInTubeOuterJacketRadius;
	const G4double dLidOuterWallDomeInTubeCutOffsetX = -dLidInTubeOffsetX;
	const G4double dLidOuterWallDomeOutTubeCutRadius = dLidOutTubeJacketRadius;
	const G4double dLidOuterWallDomeOutTubeCutOffsetX = dLidOutTubeOffsetX;

	G4Sphere *pCryostatLidOuterWallDomeSphere = new G4Sphere("CryostatLidOuterWallDomeSphere",
		dLidOuterWallDomeRadius, dLidOuterWallDomeRadius+dLidOuterWallDomeThickness,
		0.*deg, 360.*deg, 0.*deg, 90.*deg);

	G4Tubs *pCryostatLidOuterWallDomeSphereCut = new G4Tubs("CryostatLidOuterWallDomeSphereCut",
		0., dLidOuterWallRadius, 2*dLidOuterWallDomeHeight, 0.*deg, 360.*deg);

	G4IntersectionSolid *pCryostatLidOuterWallDome = new G4IntersectionSolid("CryostatLidOuterWallDome",
		pCryostatLidOuterWallDomeSphere, pCryostatLidOuterWallDomeSphereCut, 0,
		G4ThreeVector(0., 0., dLidOuterWallDomeRadius+dLidOuterWallDomeThickness));

	G4Tubs *pCryostatLidOuterWallInTubeCut = new G4Tubs("CryostatLidOuterWallInTubeCut",
		0., dLidOuterWallDomeInTubeCutRadius, dLidOuterWallDomeRadius, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidOuterWallDomeWithCut1 = new G4SubtractionSolid(
		"CryostatLidOuterWallDomeWithCut1", pCryostatLidOuterWallDome, pCryostatLidOuterWallInTubeCut,
		0, G4ThreeVector(dLidOuterWallDomeInTubeCutOffsetX, 0., 0.));

	G4Tubs *pCryostatLidOuterWallOutTubeCut = new G4Tubs("CryostatLidOuterWallOutTubeCut",
		0., dLidOuterWallDomeOutTubeCutRadius, dLidOuterWallDomeRadius, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidOuterWallDomeWithCut2 = new G4SubtractionSolid(
		"CryostatLidOuterWallWithCut2", pCryostatLidOuterWallDomeWithCut1, pCryostatLidOuterWallOutTubeCut,
		0, G4ThreeVector(dLidOuterWallDomeOutTubeCutOffsetX, 0., 0.));

	m_pCryostatLidOuterWallDomeLogicalVolume = new G4LogicalVolume(pCryostatLidOuterWallDomeWithCut2,
		SS304LSteel, "CryostatLidOuterWallDomeLogicalVolume", 0, 0, 0);

	m_pCryostatLidOuterWallDomePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLidOuterWallDomeOffsetZ),
		m_pCryostatLidOuterWallDomeLogicalVolume, "CryostatLidOuterWallDome", m_pMotherLogicalVolume, false, 0);

	//--------------------------------- lid in tube ---------------------------------
	const G4double dLidInTubeHalfZ = dLidInTubeHeight;
	const G4double dLidInTubeInnerRadius = dLidInTubeRadius-dLidInTubeThickness;
	const G4double dLidInTubeOffsetZ = dLidInnerWallOffsetZ-0.5*dLidInnerWallHeight+dLidInTubeOffsetDifference;

	const G4double dLidInTubeAngleCutBoxHalfX = 2*dLidInTubeRadius;
	const G4double dLidInTubeAngleCutBoxHalfZ = dLidInTubeRadius*sin(M_PI/8.);

	G4Tubs *pCryostatLidInTubeTubs = new G4Tubs("CryostatLidInTubeTube", dLidInTubeInnerRadius,
		dLidInTubeRadius, dLidInTubeHalfZ, 0.*deg, 360.*deg);

	G4Box *pCryostatLidInTubeAngleCut = new G4Box("CryostatLidInTubeAngleCut", dLidInTubeAngleCutBoxHalfX,
		dLidInTubeAngleCutBoxHalfX, dLidInTubeAngleCutBoxHalfZ);

	G4SubtractionSolid *pCryostatLidInTubeWithCut1 = new G4SubtractionSolid("CryostatLidInTubeWithCut1",
		pCryostatLidInTubeTubs, pCryostatLidInTubeAngleCut, m_pRotationXMinus225,
		G4ThreeVector(0., 0., dLidInTubeHalfZ));

	G4Orb *pCryostatLidInTubeBottomCut = new G4Orb("CryostatLidInTubeBottomCut", dLidInnerWallDomeRadius);

	G4SubtractionSolid *pCryostatLidInTubeWithCut2 = new G4SubtractionSolid("CryostatLidInTubeWithCut2",
		pCryostatLidInTubeWithCut1, pCryostatLidInTubeBottomCut, 0,
		G4ThreeVector(dLidInTubeOffsetX, 0., -dLidInnerWallDomeRadius-dLidInTubeOffsetDifference+dLidInnerWallDomeHeight));

	m_pCryostatLidInTubeLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeWithCut2,
		SS304LSteel, "CryostatLidInTubeLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeOffsetZ),
		m_pCryostatLidInTubeLogicalVolume, "CryostatLidInTube", m_pMotherLogicalVolume, false, 0);

	//-------------------------- lid in tube inside veto GXe ------------------------
	const G4double dLidInTubeInsideVetoGXeHalfZ = dLidInnerWallDomeTopToVetoLiquidLevel-(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference);
	const G4double dLidInTubeInsideVetoGXeOffsetZ = 0.5*dLidInnerWallDomeTopToVetoLiquidLevel-(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference);

	G4Tubs *pCryostatLidInTubeInsideVetoGXeTubs = new G4Tubs("CryostatLidInTubeInsideVetoGXeTube", dLidInTubeInnerRadius,
		dLidInTubeRadius, dLidInTubeInsideVetoGXeHalfZ, 0.*deg, 360.*deg);

	G4Orb *pCryostatLidInTubeInsideVetoGXeTopCut = new G4Orb("CryostatLidInTubeInsideVetoGXeBottomCut", dLidInnerWallDomeRadius);

	G4IntersectionSolid *pCryostatLidInTubeInsideVetoGXeWithCut1 = new G4IntersectionSolid("CryostatLidInTubeInsideVetoGXeWithCut2",
		pCryostatLidInTubeInsideVetoGXeTubs, pCryostatLidInTubeInsideVetoGXeTopCut, 0,
		G4ThreeVector(dLidInTubeOffsetX, 0., -dLidInnerWallDomeRadius-dLidInTubeOffsetDifference+dLidInnerWallDomeHeight));

	m_pCryostatLidInTubeInsideVetoGXeLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeInsideVetoGXeWithCut1,
		SS304LSteel, "CryostatLidInTubeInsideVetoGXeLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeInsideVetoGXePhysicalVolume = new G4PVPlacement(0,
		G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeInsideVetoGXeOffsetZ),
		m_pCryostatLidInTubeInsideVetoGXeLogicalVolume, "CryostatLidInTubeInsideVetoGXe",
		m_pVetoGXeLogicalVolume, false, 0);

	//------------------------- lid in tube inside veto LXe -------------------------
	const G4double dLXeHalfZ = 0.5*dLXeHeight;
	const G4double dLidInTubeInsideVetoLXeHalfZ = 0.5*(dLidInTubeHeightInside-dLidInTubeInsideVetoGXeHalfZ);
	const G4double dLidInTubeInsideVetoLXeOffsetZ = dLXeHalfZ-dLidInnerWallDomeTopToVetoLiquidLevel-dLidInTubeInsideVetoLXeHalfZ;

	G4Tubs *pCryostatLidInTubeInsideVetoLXeTubs = new G4Tubs("CryostatLidInTubeInsideVetoLXeTube", dLidInTubeInnerRadius,
		dLidInTubeRadius, dLidInTubeInsideVetoLXeHalfZ, 0.*deg, 360.*deg);

	m_pCryostatLidInTubeInsideVetoLXeLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeInsideVetoLXeTubs,
		SS304LSteel, "CryostatLidInTubeInsideVetoLXeLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeInsideVetoLXePhysicalVolume = new G4PVPlacement(0,
		G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeInsideVetoLXeOffsetZ),
		m_pCryostatLidInTubeInsideVetoLXeLogicalVolume, "CryostatLidInTubeInsideVetoLXe",
		m_pLXeLogicalVolume, false, 0);

	//----------------------------- lid in tube flange ------------------------------
	const G4double dLidInTubeFlangeHalfZ = 0.5*dLidInTubeFlangeThickness;
	const G4double dLidInTubeFlangeOffsetZ = dLidInTubeInsideVetoLXeOffsetZ-dLidInTubeInsideVetoLXeHalfZ+dLidInTubeFlangeHalfZ;

	G4Tubs *pCryostatLidInTubeFlangeTubs = new G4Tubs("CryostatLidInTubeFlangeTube", dLidInTubeInnerRadius,
		dLidInTubeFlangeRadius, dLidInTubeFlangeHalfZ, 0.*deg, 360.*deg);

	m_pCryostatLidInTubeFlangeLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeFlangeTubs,
		SS304LSteel, "CryostatLidInTubeFlangeLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeFlangePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeFlangeOffsetZ),
		m_pCryostatLidInTubeFlangeLogicalVolume, "CryostatLidInTubeFlange", m_pLXeLogicalVolume, false, 0);

	//----------------------------- lid in tube jacket ------------------------------
	const G4double dLidInTubeJacketHalfZ = dLidInTubeJacketHeight;
	const G4double dLidInTubeJacketThickness = dLidInnerWallThickness;
	const G4double dLidInTubeJacketInnerRadius = dLidInTubeJacketRadius-dLidInTubeJacketThickness;
	const G4double dLidInTubeJacketOffsetDifference = sqrt(pow(dLidInnerWallDomeRadius,2)-pow(-dLidInTubeOffsetX,2))-(dLidInnerWallDomeRadius-dLidInnerWallDomeHeight);
	const G4double dLidInTubeJacketOffsetZ = dLidInnerWallOffsetZ-0.5*dLidInnerWallHeight+dLidInTubeJacketOffsetDifference;
	
	G4Tubs *pCryostatLidInTubeJacketTubs = new G4Tubs("CryostatLidInTubeJacketTube", dLidInTubeJacketInnerRadius,
		dLidInTubeJacketRadius, dLidInTubeJacketHalfZ, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidInTubeJacketWithCut1 = new G4SubtractionSolid("CryostatLidInTubeJacketWithCut1",
		pCryostatLidInTubeJacketTubs, pCryostatLidInTubeBottomCut, 0,
		G4ThreeVector(dLidInTubeOffsetX, 0., -dLidInnerWallDomeRadius-dLidInTubeJacketOffsetDifference+dLidInnerWallDomeHeight));

	m_pCryostatLidInTubeJacketLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeJacketWithCut1,
		SS304LSteel, "CryostatLidInTubeJacketLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeJacketPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeJacketOffsetZ),
		m_pCryostatLidInTubeJacketLogicalVolume, "CryostatLidInTubeJacket", m_pMotherLogicalVolume, false, 0);

	//-------------------------- lid in tube jacket cover ---------------------------
	const G4double dLidInTubeJacketCoverOffsetZ = dLidInTubeJacketOffsetZ+dLidInTubeJacketHalfZ-0.5*dLidInTubeJacketThickness;

	G4Tubs *pCryostatLidInTubeJacketCoverTubs = new G4Tubs("CryostatLidInTubeJacketCoverTubs",
		dLidInTubeRadius, dLidInTubeJacketInnerRadius, 0.5*dLidInTubeJacketThickness, 0.*deg, 360.*deg);

	m_pCryostatLidInTubeJacketCoverLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeJacketCoverTubs,
		SS304LSteel, "CryostatLidInTubeJacketCoverLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeJacketCoverPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeJacketCoverOffsetZ),
		m_pCryostatLidInTubeJacketCoverLogicalVolume, "CryostatLidInTubeJacketCover", m_pMotherLogicalVolume, false, 0);

	//-------------------------- lid in tube outer jacket ---------------------------
	const G4double dLidInTubeOuterJacketHalfZ = dLidInTubeOuterJacketHeight;
	const G4double dLidInTubeOuterJacketThickness = dLidOuterWallThickness;
	const G4double dLidInTubeOuterJacketInnerRadius = dLidInTubeOuterJacketRadius-dLidInTubeOuterJacketThickness;
	const G4double dLidInTubeOuterJacketOffsetDifference = sqrt(pow(dLidOuterWallDomeRadius,2)-pow(dLidInTubeOffsetX,2))-(dLidOuterWallDomeRadius-dLidOuterWallDomeHeight);
	const G4double dLidInTubeOuterJacketOffsetZ = dLidOuterWallOffsetZ-0.5*dLidOuterWallHeight+dLidInTubeOuterJacketOffsetDifference;

	G4Tubs *pCryostatLidInTubeOuterJacketTubs = new G4Tubs("CryostatLidInTubeOuterJacketTube", dLidInTubeOuterJacketInnerRadius,
		dLidInTubeOuterJacketRadius, dLidInTubeOuterJacketHalfZ, 0.*deg, 360.*deg);

	G4Orb *pCryostatLidInTubeOuterJacketBottomCut = new G4Orb("CryostatLidInTubeOuterJacketBottomCut", dLidOuterWallDomeRadius);

	G4SubtractionSolid *pCryostatLidInTubeOuterJacketWithCut1 = new G4SubtractionSolid("CryostatLidInTubeOuterJacketWithCut1",
		pCryostatLidInTubeOuterJacketTubs, pCryostatLidInTubeOuterJacketBottomCut, 0,
		G4ThreeVector(dLidInTubeOffsetX, 0., -dLidOuterWallDomeRadius-dLidInTubeOuterJacketOffsetDifference+dLidOuterWallDomeHeight));

	m_pCryostatLidInTubeOuterJacketLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeOuterJacketWithCut1,
		SS304LSteel, "CryostatLidInTubeOuterJacketLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeOuterJacketPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeOuterJacketOffsetZ),
		m_pCryostatLidInTubeOuterJacketLogicalVolume, "CryostatLidInTubeOuterJacket", m_pMotherLogicalVolume, false, 0);

	//----------------------- lid in tube outer jacket cover ------------------------
	const G4double dLidInTubeOuterJacketCoverOffsetZ = dLidInTubeOuterJacketOffsetZ+dLidInTubeOuterJacketHalfZ-0.5*dLidInTubeOuterJacketThickness;

	G4Tubs *pCryostatLidInTubeOuterJacketCoverTubs = new G4Tubs("CryostatLidInTubeOuterJacketCoverTubs",
		dLidInTubeJacketRadius, dLidInTubeOuterJacketInnerRadius, 0.5*dLidInTubeOuterJacketThickness, 0.*deg, 360.*deg);

	m_pCryostatLidInTubeOuterJacketCoverLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeOuterJacketCoverTubs,
		SS304LSteel, "CryostatLidInTubeOuterJacketCoverLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeOuterJacketCoverPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(-dLidInTubeOffsetX, 0., dLidInTubeOuterJacketCoverOffsetZ),
		m_pCryostatLidInTubeOuterJacketCoverLogicalVolume, "CryostatLidInTubeOuterJacketCover", m_pMotherLogicalVolume, false, 0);

	//------------------------------- lid in tube elbow -----------------------------
	const G4double dLidInTubeElbowOffsetY = -(M_SQRT1_2*dLidInTubeRadius+M_SQRT1_2*0.5*dLidInTubeElbowLength-dLidInTubeRadius);
	const G4double dLidInTubeElbowOffsetZ = dLidInTubeHalfZ-M_SQRT1_2*dLidInTubeRadius+M_SQRT1_2*0.5*dLidInTubeElbowLength;

	G4Tubs *pCryostatLidInTubeElbowTubs = new G4Tubs("CryostatLidInTubeElbowTube", dLidInTubeInnerRadius,
		dLidInTubeRadius, 0.5*dLidInTubeElbowLength, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidInTubeElbowWithCut1 = new G4SubtractionSolid("CryostatLidInTubeElbowWithCut1",
		pCryostatLidInTubeElbowTubs, pCryostatLidInTubeAngleCut, m_pRotationXMinus225,
		G4ThreeVector(0., 0., 0.5*dLidInTubeElbowLength));

	G4SubtractionSolid *pCryostatLidInTubeElbowWithCut2 = new G4SubtractionSolid("CryostatLidInTubeElbowWithCut2",
		pCryostatLidInTubeElbowWithCut1, pCryostatLidInTubeAngleCut, m_pRotationXPlus225,
		G4ThreeVector(0., 0., -0.5*dLidInTubeElbowLength));

	m_pCryostatLidInTubeElbowLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeElbowWithCut2,
		SS304LSteel, "CryostatLidInTubeElbowLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeElbowPhysicalVolume = new G4PVPlacement(m_pRotationXMinus45,
		G4ThreeVector(-dLidInTubeOffsetX, dLidInTubeElbowOffsetY, dLidInTubeOffsetZ+dLidInTubeElbowOffsetZ),
		m_pCryostatLidInTubeElbowLogicalVolume, "CryostatLidInTubeElbow", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- lid in tube top ------------------------------
	const G4double dLidInTubeTopOffsetY = -(M_SQRT1_2*dLidInTubeElbowLength+0.5*dLidInTubeTopLength-dLidInTubeRadius);
	const G4double dLidInTubeTopOffsetZ = dLidInTubeHalfZ+M_SQRT1_2*dLidInTubeElbowLength-dLidInTubeRadius;

	G4Tubs *pCryostatLidInTubeTopTubs = new G4Tubs("CryostatLidInTubeTopTube", dLidInTubeInnerRadius,
		dLidInTubeRadius, 0.5*dLidInTubeTopLength, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidInTubeTopWithCut1 = new G4SubtractionSolid("CryostatLidInTubeTopWithCut1",
		pCryostatLidInTubeTopTubs, pCryostatLidInTubeAngleCut, m_pRotationXPlus225,
		G4ThreeVector(0., 0., -0.5*dLidInTubeTopLength));

	m_pCryostatLidInTubeTopLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeTopWithCut1,
		SS304LSteel, "CryostatLidInTubeTopLogicalVolume", 0, 0, 0);

	m_pCryostatLidInTubeTopPhysicalVolume = new G4PVPlacement(m_pRotationXMinus90,
		G4ThreeVector(-dLidInTubeOffsetX, dLidInTubeTopOffsetY, dLidInTubeOffsetZ+dLidInTubeTopOffsetZ),
		m_pCryostatLidInTubeTopLogicalVolume, "CryostatLidInTubeTop", m_pMotherLogicalVolume, false, 0);

	//------------------------------- lid in tube top flange ------------------------------
	const G4double dLidInTubeTopFlangeOffsetY = dLidInTubeTopOffsetY-0.5*dLidInTubeTopLength-0.5*dLidInTubeTopFlangeThickness;
	const G4double dLidInTubeTopFlangeOffsetZ = dLidInTubeTopOffsetZ;

	G4Tubs *pCryostatLidInTubeTopFlangeTubs = new G4Tubs("CryostatLidInTubeTopFlangeTube", dLidInTubeInnerRadius, dLidInTubeTopFlangeRadius, 0.5*dLidInTubeTopFlangeThickness, 0.*deg, 360.*deg);
	m_pCryostatLidInTubeTopFlangeLogicalVolume = new G4LogicalVolume(pCryostatLidInTubeTopFlangeTubs, SS304LSteel, "CryostatLidInTubeTopFlangeLogicalVolume", 0, 0, 0);
	m_pCryostatLidInTubeTopFlangePhysicalVolume = new G4PVPlacement(m_pRotationXMinus90,G4ThreeVector(-dLidInTubeOffsetX, dLidInTubeTopFlangeOffsetY, dLidInTubeOffsetZ+dLidInTubeTopFlangeOffsetZ), m_pCryostatLidInTubeTopFlangeLogicalVolume, "CryostatLidInTubeTopFlange", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- lid in tube top to shield ------------------------------
	const G4double dPolishLeadHalfY 	= dShieldHalfY;
	const G4double dFrenchLeadHalfY 	= dPolishLeadHalfY - dPolishLeadThickness;
	const G4double dPolyethyleneHalfY 	= dFrenchLeadHalfY - dFrenchLeadThickness;
	const G4double dCopperHalfY 		= dPolyethyleneHalfY - dPolyethyleneThickness;
	const G4double dCavityHalfY 		= dCopperHalfY - GetGeometryParameter("dCopperThicknessThick");

	const G4double dLidInTubeTopToShieldLength = dCavityHalfY + (dLidInTubeTopFlangeOffsetY-0.5*dLidInTubeTopFlangeThickness) - dCopperSideThickness;
	const G4double dLidInTubeTopToShieldOffsetX = -dLidInTubeOffsetX;
	const G4double dLidInTubeTopToShieldOffsetY = dLidInTubeTopFlangeOffsetY-0.5*dLidInTubeTopFlangeThickness-0.5*dLidInTubeTopToShieldLength;
	const G4double dLidInTubeTopToShieldOffsetZ = dLidInTubeOffsetZ+dLidInTubeTopOffsetZ;

	G4Tubs *pCryostatLidInTubeTopToShieldTube 		= new G4Tubs("CryostatLidInTubeTopToShieldTube", dLidInTubeInnerRadius, dLidInTubeRadius, 0.5*dLidInTubeTopToShieldLength, 0.*deg, 360.*deg);
	m_pCryostatLidInTubeTopToShieldLogicalVolume 	= new G4LogicalVolume(pCryostatLidInTubeTopToShieldTube, SS304LSteel, "CryostatLidInTubeTopToShieldLogicalVolume", 0, 0, 0);
	m_pCryostatLidInTubeTopToShieldPhysicalVolume 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(dLidInTubeTopToShieldOffsetX, dLidInTubeTopToShieldOffsetY, dLidInTubeTopToShieldOffsetZ), m_pCryostatLidInTubeTopToShieldLogicalVolume, "CryostatLidInTubeTopToShield", m_pMotherLogicalVolume, false, 0);

	const G4double dLidInTubeTopToShieldOutsideLength 	= 10 *cm;
	const G4double dLidInTubeTopToShieldOutsideOffsetX 	= dLidInTubeTopToShieldOffsetX;
	const G4double dLidInTubeTopToShieldOutsideOffsetY 	= -dShieldHalfY - dLidInTubeTopToShieldOutsideLength/2;
	const G4double dLidInTubeTopToShieldOutsideOffsetZ 	= dLidInTubeTopToShieldOffsetZ-21.5*mm;

	G4Tubs *pCryostatLidInTubeTopToShieldTubeOutside 		= new G4Tubs("CryostatLidInTubeTopToShieldTubeOutside", dLidInTubeInnerRadius, dLidInTubeRadius, 0.5*dLidInTubeTopToShieldOutsideLength, 0.*deg, 360.*deg);
	m_pCryostatLidInTubeTopToShieldOutsideLogicalVolume 	= new G4LogicalVolume(pCryostatLidInTubeTopToShieldTubeOutside, SS304LSteel, "CryostatLidInTubeTopToShieldLogicalVolume", 0, 0, 0);
	m_pCryostatLidInTubeTopToShieldOutsidePhysicalVolume 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(dLidInTubeTopToShieldOutsideOffsetX, dLidInTubeTopToShieldOutsideOffsetY, dLidInTubeTopToShieldOutsideOffsetZ), m_pCryostatLidInTubeTopToShieldOutsideLogicalVolume, "CryostatLidInTubeTopToShieldOutside", m_pLabLogicalVolume, false, 0);

	//--------------------------------- lid out tube --------------------------------
	const G4double dLidOutTubeHalfZ = dLidOutTubeHeight;
	const G4double dLidOutTubeInnerRadius = dLidOutTubeRadius-dLidOutTubeThickness;
	const G4double dLidOutTubeOffsetDifference = sqrt(pow(dLidInnerWallDomeRadius,2)-pow(dLidOutTubeOffsetX,2))-(dLidInnerWallDomeRadius-dLidInnerWallDomeHeight);
	const G4double dLidOutTubeOffsetZ = dLidInnerWallOffsetZ-0.5*dLidInnerWallHeight+dLidOutTubeOffsetDifference;

	const G4double dLidOutTubeAngleCutBoxHalfX = 2*dLidOutTubeRadius;
	const G4double dLidOutTubeAngleCutBoxHalfZ = dLidOutTubeRadius*sin(M_PI/8.);

	G4Tubs *pCryostatLidOutTubeTubs = new G4Tubs("CryostatLidOutTubeTube", dLidOutTubeInnerRadius,
		dLidOutTubeRadius, dLidOutTubeHalfZ, 0.*deg, 360.*deg);

	G4Box *pCryostatLidOutTubeAngleCut = new G4Box("CryostatLidOutTubeAngleCut", dLidOutTubeAngleCutBoxHalfX,
		dLidOutTubeAngleCutBoxHalfX, dLidOutTubeAngleCutBoxHalfZ);

	G4SubtractionSolid *pCryostatLidOutTubeWithCut1 = new G4SubtractionSolid("CryostatLidOutTubeWithCut1",
		pCryostatLidOutTubeTubs, pCryostatLidOutTubeAngleCut, m_pRotationXMinus225,
		G4ThreeVector(0., 0., dLidOutTubeHalfZ));

	G4Orb *pCryostatLidOutTubeBottomCut = new G4Orb("CryostatLidOutTubeBottomCut", dLidInnerWallDomeRadius);

	G4SubtractionSolid *pCryostatLidOutTubeWithCut2 = new G4SubtractionSolid("CryostatLidOutTubeWithCut2", pCryostatLidOutTubeWithCut1, pCryostatLidOutTubeBottomCut, 0, G4ThreeVector(-dLidOutTubeOffsetX, 0., -dLidInnerWallDomeRadius-dLidOutTubeOffsetDifference+dLidInnerWallDomeHeight));

	m_pCryostatLidOutTubeLogicalVolume = new G4LogicalVolume(pCryostatLidOutTubeWithCut2, SS304LSteel, "CryostatLidOutTubeLogicalVolume", 0, 0, 0);
	m_pCryostatLidOutTubePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dLidOutTubeOffsetX, 0., dLidOutTubeOffsetZ), m_pCryostatLidOutTubeLogicalVolume, "CryostatLidOutTube", m_pMotherLogicalVolume, false, 0);

	//----------------------------- lid out tube jacket -----------------------------
	const G4double dLidOutTubeJacketHalfZ = dLidOutTubeJacketHeight;
	const G4double dLidOutTubeJacketThickness = dLidOuterWallThickness;
	const G4double dLidOutTubeJacketInnerRadius = dLidOutTubeJacketRadius-dLidOutTubeJacketThickness;
	const G4double dLidOutTubeJacketOffsetDifference = sqrt(pow(dLidOuterWallDomeRadius,2)-pow(dLidOutTubeOffsetX,2))-(dLidOuterWallDomeRadius-dLidOuterWallDomeHeight);
	const G4double dLidOutTubeJacketOffsetZ = dLidOuterWallOffsetZ-0.5*dLidOuterWallHeight+dLidOutTubeJacketOffsetDifference;

	G4Tubs *pCryostatLidOutTubeJacketTubs = new G4Tubs("CryostatLidOutTubeJacketTube", dLidOutTubeJacketInnerRadius, dLidOutTubeJacketRadius, dLidOutTubeJacketHalfZ, 0.*deg, 360.*deg);

	G4Orb *pCryostatLidOutTubeJacketBottomCut = new G4Orb("CryostatLidOutTubeJacketBottomCut", dLidOuterWallDomeRadius);

	G4SubtractionSolid *pCryostatLidOutTubeJacketWithCut1 = new G4SubtractionSolid("CryostatLidOutTubeJacketWithCut1", pCryostatLidOutTubeJacketTubs, pCryostatLidOutTubeJacketBottomCut, 0, G4ThreeVector(-dLidOutTubeOffsetX, 0., -dLidOuterWallDomeRadius-dLidOutTubeJacketOffsetDifference+dLidOuterWallDomeHeight));

	m_pCryostatLidOutTubeJacketLogicalVolume = new G4LogicalVolume(pCryostatLidOutTubeJacketWithCut1, SS304LSteel, "CryostatLidOutTubeJacketLogicalVolume", 0, 0, 0);
	m_pCryostatLidOutTubeJacketPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dLidOutTubeOffsetX, 0., dLidOutTubeJacketOffsetZ), m_pCryostatLidOutTubeJacketLogicalVolume, "CryostatLidOutTubeJacket", m_pMotherLogicalVolume, false, 0);

	//-------------------------- lid out tube jacket cover --------------------------
	const G4double dLidOutTubeJacketCoverOffsetZ = dLidOutTubeJacketOffsetZ+dLidOutTubeJacketHalfZ-0.5*dLidOutTubeJacketThickness;

	G4Tubs *pCryostatLidOutTubeJacketCoverTubs = new G4Tubs("CryostatLidOutTubeJacketCoverTubs",
		dLidOutTubeRadius, dLidOutTubeJacketInnerRadius, 0.5*dLidOutTubeJacketThickness, 0.*deg, 360.*deg);

	m_pCryostatLidOutTubeJacketCoverLogicalVolume = new G4LogicalVolume(pCryostatLidOutTubeJacketCoverTubs,
		SS304LSteel, "CryostatLidOutTubeJacketCoverLogicalVolume", 0, 0, 0);

	m_pCryostatLidOutTubeJacketCoverPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dLidOutTubeOffsetX, 0., dLidOutTubeJacketCoverOffsetZ),
		m_pCryostatLidOutTubeJacketCoverLogicalVolume, "CryostatLidOutTubeJacketCover", m_pMotherLogicalVolume, false, 0);

//    G4Tubs *pTempTubs = new G4Tubs("TempTube", 0.*cm, 10.*cm, 1.*mm, 0.*deg, 360.*deg);

//    G4LogicalVolume *pTempLogicalVolume = new G4LogicalVolume(pTempTubs,
//        SS304LSteel, "TempLogicalVolume", 0, 0, 0);
//    pTempLogicalVolume->SetVisAttributes(G4VisAttributes(G4Colour(1,0,0)));

//    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
//        pTempLogicalVolume, "Temp", m_pMotherLogicalVolume, false, 0);

	//------------------------------ lid out tube elbow -----------------------------
	const G4double dLidOutTubeElbowOffsetY = -(M_SQRT1_2*dLidOutTubeRadius+M_SQRT1_2*0.5*dLidOutTubeElbowLength-dLidOutTubeRadius);
	const G4double dLidOutTubeElbowOffsetZ = dLidOutTubeHalfZ-M_SQRT1_2*dLidOutTubeRadius+M_SQRT1_2*0.5*dLidOutTubeElbowLength;

	G4Tubs *pCryostatLidOutTubeElbowTubs = new G4Tubs("CryostatLidOutTubeElbowTube", dLidOutTubeInnerRadius,
		dLidOutTubeRadius, 0.5*dLidOutTubeElbowLength, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidOutTubeElbowWithCut1 = new G4SubtractionSolid("CryostatLidOutTubeElbowWithCut1",
		pCryostatLidOutTubeElbowTubs, pCryostatLidOutTubeAngleCut, m_pRotationXMinus225,
		G4ThreeVector(0., 0., 0.5*dLidOutTubeElbowLength));

	G4SubtractionSolid *pCryostatLidOutTubeElbowWithCut2 = new G4SubtractionSolid("CryostatLidOutTubeElbowWithCut2",
		pCryostatLidOutTubeElbowWithCut1, pCryostatLidOutTubeAngleCut, m_pRotationXPlus225,
		G4ThreeVector(0., 0., -0.5*dLidOutTubeElbowLength));

	m_pCryostatLidOutTubeElbowLogicalVolume = new G4LogicalVolume(pCryostatLidOutTubeElbowWithCut2,
		SS304LSteel, "CryostatLidOutTubeElbowLogicalVolume", 0, 0, 0);

	m_pCryostatLidOutTubeElbowPhysicalVolume = new G4PVPlacement(m_pRotationXMinus45,
		G4ThreeVector(dLidOutTubeOffsetX, dLidOutTubeElbowOffsetY, dLidOutTubeOffsetZ+dLidOutTubeElbowOffsetZ),
		m_pCryostatLidOutTubeElbowLogicalVolume, "CryostatLidOutTubeElbow", m_pMotherLogicalVolume, false, 0);

	//------------------------------- lid out tube top ------------------------------
	const G4double dLidOutTubeTopOffsetY = -(M_SQRT1_2*dLidOutTubeElbowLength+0.5*dLidOutTubeTopLength-dLidOutTubeRadius);
	const G4double dLidOutTubeTopOffsetZ = dLidOutTubeHalfZ+M_SQRT1_2*dLidOutTubeElbowLength-dLidOutTubeRadius;

	G4Tubs *pCryostatLidOutTubeTopTubs = new G4Tubs("CryostatLidOutTubeTopTube", dLidOutTubeInnerRadius,
		dLidOutTubeRadius, 0.5*dLidOutTubeTopLength, 0.*deg, 360.*deg);

	G4SubtractionSolid *pCryostatLidOutTubeTopWithCut1 = new G4SubtractionSolid("CryostatLidOutTubeTopWithCut1",
		pCryostatLidOutTubeTopTubs, pCryostatLidOutTubeAngleCut, m_pRotationXPlus225,
		G4ThreeVector(0., 0., -0.5*dLidOutTubeTopLength));

	m_pCryostatLidOutTubeTopLogicalVolume = new G4LogicalVolume(pCryostatLidOutTubeTopWithCut1,
		SS304LSteel, "CryostatLidOutTubeTopLogicalVolume", 0, 0, 0);

	m_pCryostatLidOutTubeTopPhysicalVolume = new G4PVPlacement(m_pRotationXMinus90,
		G4ThreeVector(dLidOutTubeOffsetX, dLidOutTubeTopOffsetY, dLidOutTubeOffsetZ+dLidOutTubeTopOffsetZ),
		m_pCryostatLidOutTubeTopLogicalVolume, "CryostatLidOutTubeTop", m_pMotherLogicalVolume, false, 0);

	//------------------------------- lid out tube top flange ------------------------------
	const G4double dLidOutTubeTopFlangeOffsetY = dLidOutTubeTopOffsetY-0.5*dLidOutTubeTopLength-0.5*dLidOutTubeTopFlangeThickness;
	const G4double dLidOutTubeTopFlangeOffsetZ = dLidOutTubeTopOffsetZ;

	G4Tubs *pCryostatLidOutTubeTopFlangeTubs 		= new G4Tubs("CryostatLidOutTubeTopFlangeTube", dLidOutTubeInnerRadius, dLidOutTubeTopFlangeRadius, 0.5*dLidOutTubeTopFlangeThickness, 0.*deg, 360.*deg);
	m_pCryostatLidOutTubeTopFlangeLogicalVolume 	= new G4LogicalVolume(pCryostatLidOutTubeTopFlangeTubs, SS304LSteel, "CryostatLidOutTubeTopFlangeLogicalVolume", 0, 0, 0);
	m_pCryostatLidOutTubeTopFlangePhysicalVolume	= new G4PVPlacement(m_pRotationXMinus90,G4ThreeVector(dLidOutTubeOffsetX, dLidOutTubeTopFlangeOffsetY, dLidOutTubeOffsetZ+dLidOutTubeTopFlangeOffsetZ), m_pCryostatLidOutTubeTopFlangeLogicalVolume, "CryostatLidOutTubeTopFlange", m_pMotherLogicalVolume, false, 0);

	//-------------------------------- lid out tube top to shield ------------------------------
	const G4double dLidOutTubeTopToShieldLength 	= dCavityHalfY + (dLidOutTubeTopFlangeOffsetY-0.5*dLidOutTubeTopFlangeThickness) - dCopperSideThickness;
	const G4double dLidOutTubeTopToShieldOffsetX 	= dLidOutTubeOffsetX;
	const G4double dLidOutTubeTopToShieldOffsetY 	= dLidOutTubeTopFlangeOffsetY-0.5*dLidOutTubeTopFlangeThickness-0.5*dLidOutTubeTopToShieldLength;
	const G4double dLidOutTubeTopToShieldOffsetZ 	= dLidOutTubeOffsetZ+dLidOutTubeTopOffsetZ;

	G4Tubs *pCryostatLidOutTubeTopToShieldTube 		= new G4Tubs("CryostatLidOutTubeTopToShieldTube", dLidOutTubeInnerRadius, dLidOutTubeRadius, 0.5*dLidOutTubeTopToShieldLength, 0.*deg, 360.*deg);
	m_pCryostatLidOutTubeTopToShieldLogicalVolume 	= new G4LogicalVolume(pCryostatLidOutTubeTopToShieldTube, SS304LSteel, "CryostatLidOutTubeTopToShieldLogicalVolume", 0, 0, 0);
	m_pCryostatLidOutTubeTopToShieldPhysicalVolume 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(dLidOutTubeTopToShieldOffsetX, dLidOutTubeTopToShieldOffsetY, dLidOutTubeOffsetZ+dLidOutTubeTopOffsetZ), m_pCryostatLidOutTubeTopToShieldLogicalVolume, "CryostatLidOutTubeTopToShield", m_pMotherLogicalVolume, false, 0);

	const G4double dLidOutTubeTopToShieldOutsideLength 	= 10 *cm;
	const G4double dLidOutTubeTopToShieldOutsideOffsetX = dLidOutTubeTopToShieldOffsetX;
	const G4double dLidOutTubeTopToShieldOutsideOffsetY = -dShieldHalfY - dLidOutTubeTopToShieldOutsideLength/2;
	const G4double dLidOutTubeTopToShieldOutsideOffsetZ = dLidOutTubeTopToShieldOffsetZ-21.5*mm;

	G4Tubs *pCryostatLidOutTubeTopToShieldTubeOutside 		= new G4Tubs("CryostatLidOutTubeTopToShieldTubeOutside", dLidOutTubeInnerRadius, dLidOutTubeRadius, 0.5*dLidOutTubeTopToShieldOutsideLength, 0.*deg, 360.*deg);
	m_pCryostatLidOutTubeTopToShieldOutsideLogicalVolume 	= new G4LogicalVolume(pCryostatLidOutTubeTopToShieldTubeOutside, SS304LSteel, "CryostatLidOutTubeTopToShieldOutsideLogicalVolume", 0, 0, 0);
	m_pCryostatLidOutTubeTopToShieldOutsidePhysicalVolume 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(dLidOutTubeTopToShieldOutsideOffsetX, dLidOutTubeTopToShieldOutsideOffsetY, dLidOutTubeTopToShieldOutsideOffsetZ), m_pCryostatLidOutTubeTopToShieldOutsideLogicalVolume, "CryostatLidOutTubeTopToShieldOutside", m_pLabLogicalVolume, false, 0);

	//---------------------------------- lid vacuum ---------------------------------
//    G4Tubs *pVacuumCryostatLidTubs = new G4Tubs("VacuumCryostatLidTube", 0.*cm, dLidRadius1,
//        dVacuumLidHalfZ, 0.*deg, 360.*deg);

//    m_pVacuumCryostatLidLogicalVolume = new G4LogicalVolume(pVacuumCryostatLidTubs, Vacuum,
//        "VacuumCryostatLidLogicalVolume", 0, 0, 0);

//    m_pVacuumCryostatLidPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dVacuumLidOffsetZ),
//        m_pVacuumCryostatLidLogicalVolume, "VacuumCryostatLid", m_pMotherLogicalVolume, false, 0);


	// Cryostat Side Tube [reference P20008]
	//const G4double inside_depth = GetGeometryParameter("CavityDepth");
    G4RotationMatrix RotM90y;
	G4RotationMatrix RotP90y;
	G4RotationMatrix RotM90x;
	G4RotationMatrix RotP90x;
	G4RotationMatrix RotP0x;

	RotP0x.rotateX   (0. *deg);
	RotP90x.rotateX (90. *deg);
	RotM90x.rotateX(-90. *deg);
	RotP90y.rotateY (90. *deg);
	RotM90y.rotateY(-90. *deg);

	// outer
	G4double CST_thick	= 1.5 *mm;

	G4double CSTout_oD	= 76.20 *mm;
	G4double CSTout_iD	= CSTout_oD - 2*CST_thick;
	G4double pRtorOut	= 40. *mm;
	G4Torus* CSToutM	= new G4Torus("CSToutM", CSTout_iD*0.5, CSTout_oD*0.5, pRtorOut, 0, 90.*deg);
	//G4double CSToutL_length	= inside_depth*0.5 - pRtorOut;
	//G4double CSToutL_length	= inside_depth*0.5;
	G4double CSToutL_length	= dCavityHalfY - pRtorOut - dCopperSideThickness;
	G4Tubs*	 CSToutL	= new G4Tubs("CSToutL", CSTout_iD*0.5, CSTout_oD*0.5, CSToutL_length*0.5,  0.*deg, 360.*deg);
	G4double CSToutU_length	= 50*mm;
	G4Tubs*	 CSToutU	= new G4Tubs("CSToutU", CSTout_iD*0.5, CSTout_oD*0.5, CSToutU_length*0.5,  0.*deg, 360.*deg);

	G4ThreeVector CSToutL_V(pRtorOut,				-CSToutL_length*0.5,	0.0);
	G4ThreeVector CSToutU_V(-CSToutU_length*0.5,	pRtorOut,				0.0);

	G4Transform3D CSToutL_T(RotP90x, CSToutL_V);
	G4Transform3D CSToutU_T(RotP90y, CSToutU_V);

	G4UnionSolid* CSTout1	= new G4UnionSolid("CSTout1", CSToutM, CSToutL, CSToutL_T);
	G4UnionSolid* CSTout	= new G4UnionSolid("CSTout", CSTout1, CSToutU, CSToutU_T);

	// inner
	G4double CSTin_oD	= 50.93 *mm;
	G4double CSTin_iD	= CSTin_oD - 2*CST_thick;
	G4double pRtorIn	= 40. *mm;
	G4Torus* CSTinM	= new G4Torus("CSTinM", CSTin_iD*0.5, CSTin_oD*0.5, pRtorIn, 0, 90.*deg);
	//G4double CSTinL_length	= 50*mm;
	//G4double CSTinL_length	= inside_depth*0.5;
	G4double CSTinL_length	= dCavityHalfY - pRtorIn - dCopperSideThickness;
	G4Tubs*	 CSTinL	= new G4Tubs("CSTinL", CSTin_iD*0.5, CSTin_oD*0.5, CSTinL_length*0.5,  0.*deg, 360.*deg);
	G4double CSTinU_length	= 50*mm;
	G4Tubs*	 CSTinU	= new G4Tubs("CSTinU", CSTin_iD*0.5, CSTin_oD*0.5, CSTinU_length*0.5,  0.*deg, 360.*deg);

	G4ThreeVector CSTinL_V(pRtorIn,				-CSTinL_length*0.5,		0.0);
	G4ThreeVector CSTinU_V(-CSTinU_length*0.5,	pRtorIn,				0.0);

	G4Transform3D CSTinL_T(RotP90x, CSTinL_V);
	G4Transform3D CSTinU_T(RotP90y, CSTinU_V);

	G4UnionSolid* CSTin1	= new G4UnionSolid("CSTin1",	CSTinM, CSTinL, CSTinL_T);
	G4UnionSolid* CSTin		= new G4UnionSolid("CSTin",		CSTin1, CSTinU, CSTinU_T);

	G4LogicalVolume* CryoSideTubeOut_log	= new G4LogicalVolume(CSTout,	SS304LSteel,	"CryoSideTubeOut_log");
	G4LogicalVolume* CryoSideTubeIn_log		= new G4LogicalVolume(CSTin,	SS304LSteel,	"CryoSideTubeIn_log");

	G4double CSTout_X = dVesselOuterWallRadius + CSToutU_length;
	G4double CSTout_Y = -pRtorOut;
	G4double CSTout_Z = dVesselOuterWallOffsetZ + dVesselOuterWallHalfZ - 69.29*mm;
	G4PVPlacement* CryoSideTubeOut_phys	= new G4PVPlacement(0, G4ThreeVector(CSTout_X, CSTout_Y, CSTout_Z), CryoSideTubeOut_log, "CryoSideTubeOut_phys", m_pMotherLogicalVolume, false, 0);
	G4double CSTin_X = dVesselOuterWallRadius + CSTinU_length;
	G4double CSTin_Y = -pRtorIn;
	G4double CSTin_Z = dVesselOuterWallOffsetZ + dVesselOuterWallHalfZ - 69.29*mm;
	G4PVPlacement* CryoSideTubeIn_phys	= new G4PVPlacement(0, G4ThreeVector(CSTin_X, CSTin_Y, CSTin_Z), CryoSideTubeIn_log, "CryoSideTubeIn_phys", m_pMotherLogicalVolume, false, 0);
	
	// PORCUPINES
	G4double Porcupine_iD 		= 0.0;
	G4double Porcupine_oD 		= 15 *cm;
	G4double Porcupine_thick 	= 10 *cm;
	G4Tubs*	 Porcupine	= new G4Tubs("Porcupine", Porcupine_iD*0.5, Porcupine_oD*0.5, Porcupine_thick*0.5,  0.*deg, 360.*deg);
	G4LogicalVolume* PorcupineLeft_log	= new G4LogicalVolume(Porcupine,	SS304LSteel,	"PorcupineLeft_log");
	G4LogicalVolume* PorcupineRight_log	= new G4LogicalVolume(Porcupine,	SS304LSteel,	"PorcupineRight_log");

	const G4double PorcupineLeftOffsetX = dLidInTubeTopToShieldOutsideOffsetX;
	const G4double PorcupineLeftOffsetY = dLidInTubeTopToShieldOutsideOffsetY-dLidInTubeTopToShieldOutsideLength/2-Porcupine_thick/2;
	const G4double PorcupineLeftOffsetZ = dLidInTubeTopToShieldOutsideOffsetZ;

	const G4double PorcupineRightOffsetX = dLidOutTubeTopToShieldOutsideOffsetX;
	const G4double PorcupineRightOffsetY = dLidOutTubeTopToShieldOutsideOffsetY-dLidOutTubeTopToShieldOutsideLength/2-Porcupine_thick/2;
	const G4double PorcupineRightOffsetZ = dLidOutTubeTopToShieldOutsideOffsetZ;

	G4PVPlacement *PorcupineLeft_phys 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(PorcupineLeftOffsetX, PorcupineLeftOffsetY, PorcupineLeftOffsetZ), PorcupineLeft_log, "PorcupineLeft_phys", m_pLabLogicalVolume, false, 0);
	G4PVPlacement *PorcupineRight_phys 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(PorcupineRightOffsetX, PorcupineRightOffsetY, PorcupineRightOffsetZ), PorcupineRight_log, "PorcupineRight_phys", m_pLabLogicalVolume, false, 0);
	
	// PTR
	G4double PTR_thick 	= 2.0 *mm;
	G4double PTR_oD		= 23.0*cm;
	G4double PTR_iD		= PTR_oD - 2*PTR_thick;
	G4double PTR_height	= 30*cm;

	G4double PTRflange_thick 	= 2.0 *cm;
	G4double PTRflange_iD 		= 0.0;
	G4double PTRflange_oD 		= PTR_oD + 3.0*cm;

	G4double PTRpipe_thick 	= CST_thick;
	G4double PTRpipe_oD		= CSTout_oD;
	G4double PTRpipe_iD		= CSTout_iD;
	G4double PTRpipe_length	= 20.0 *cm;

	G4double PTRangle_oD	= PTRpipe_oD;
	G4double PTRangle_iD	= PTRpipe_iD;
	G4double PTRangle_Rtor	= 50. *mm;


	G4Tubs*	 PTRpipe	= new G4Tubs("PTRpipe", PTRpipe_iD*0.5, PTRpipe_oD*0.5, PTRpipe_length*0.5,  0.*deg, 360.*deg);
	G4Torus* PTRangle	= new G4Torus("PTRangle", PTRangle_iD*0.5, PTRangle_oD*0.5, PTRangle_Rtor, 0, 90.*deg);
	G4Tubs*	 PTRflange	= new G4Tubs("PTRflange", PTRflange_iD*0.5, PTRflange_oD*0.5, PTRflange_thick*0.5,  0.*deg, 360.*deg);
	G4Tubs*	 PTR		= new G4Tubs("PTR", PTR_iD*0.5, PTR_oD*0.5, PTR_height*0.5,  0.*deg, 360.*deg);

	G4LogicalVolume* PTRpipe_log	= new G4LogicalVolume(PTRpipe,		SS304LSteel,	"PTRpipe_log");
	G4LogicalVolume* PTRangle_log	= new G4LogicalVolume(PTRangle,		SS304LSteel,	"PTRangle_log");
	G4LogicalVolume* PTRflange_log	= new G4LogicalVolume(PTRflange,	SS304LSteel,	"PTRflange_log");
	G4LogicalVolume* PTR_log		= new G4LogicalVolume(PTR,			SS304LSteel,	"PTR_log");
	
	const G4double PTRpipe_OffsetX = CSTout_X;
	const G4double PTRpipe_OffsetY = -dShieldHalfY - PTRpipe_length/2;
	const G4double PTRpipe_OffsetZ = CSTout_Z - 21.5*mm;

	const G4double PTRangle_OffsetX = PTRpipe_OffsetX;
	const G4double PTRangle_OffsetY = PTRpipe_OffsetY - PTRpipe_length/2;
	const G4double PTRangle_OffsetZ = PTRpipe_OffsetZ + PTRangle_Rtor;

	const G4double PTRflange_OffsetX = PTRpipe_OffsetX;
	const G4double PTRflange_OffsetY = PTRpipe_OffsetY - PTRpipe_length/2 - PTRangle_Rtor;
	const G4double PTRflange_OffsetZ = PTRpipe_OffsetZ + PTRangle_Rtor + PTRflange_thick/2;

	const G4double PTR_OffsetX = PTRflange_OffsetX;
	const G4double PTR_OffsetY = PTRflange_OffsetY;
	const G4double PTR_OffsetZ = PTRflange_OffsetZ + PTRflange_thick/2 + PTR_height/2;
	
	G4PVPlacement *PTRpipe_phys 	= new G4PVPlacement(m_pRotationXMinus90, G4ThreeVector(PTRpipe_OffsetX, PTRpipe_OffsetY, PTRpipe_OffsetZ), PTRpipe_log, "PTRpipe_phys", m_pLabLogicalVolume, false, 0);
	G4PVPlacement *PTRangle_phys 	= new G4PVPlacement(m_pRotationXPlus90YMinus90, G4ThreeVector(PTRangle_OffsetX, PTRangle_OffsetY, PTRangle_OffsetZ), PTRangle_log, "PTRangle_phys", m_pLabLogicalVolume, false, 0);
	G4PVPlacement *PTRflange_phys 	= new G4PVPlacement(0, G4ThreeVector(PTRflange_OffsetX, PTRflange_OffsetY, PTRflange_OffsetZ), PTRflange_log, "PTRflange_phys", m_pLabLogicalVolume, false, 0);
	G4PVPlacement *PTR_phys 		= new G4PVPlacement(0, G4ThreeVector(PTR_OffsetX, PTR_OffsetY, PTR_OffsetZ), PTR_log, "PTR_phys", m_pLabLogicalVolume, false, 0);
	
	
	//================================== attributes =================================
	//G4Colour hStainlessSteelColor(0.500, 0.500, 0.500, 0.1);
	//G4Colour hVacuumColor(0.600, 0.200, 0.250, 0.1);
	//G4Colour hSteelColour(0.0, 0.0, 1.0);		// xlblue
	G4Colour hSteelColour(0.0, 0.0, 0.0);
	
	G4VisAttributes *pCryostatVesselFlangeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatVesselFlangeVisAtt->SetVisibility(true);	//f
	m_pCryostatVesselFlangeLogicalVolume->SetVisAttributes(pCryostatVesselFlangeVisAtt);

	G4VisAttributes *pCryostatVesselInnerWallVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatVesselInnerWallVisAtt->SetVisibility(true); //f
	m_pCryostatVesselInnerWallLogicalVolume->SetVisAttributes(pCryostatVesselInnerWallVisAtt);

	G4VisAttributes *pCryostatVesselInnerWallDomeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatVesselInnerWallDomeVisAtt->SetVisibility(true);	//f
	m_pCryostatVesselInnerWallDomeLogicalVolume->SetVisAttributes(pCryostatVesselInnerWallDomeVisAtt);

	G4VisAttributes *pCryostatVesselOuterWallVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatVesselOuterWallVisAtt->SetVisibility(true);	//f
	m_pCryostatVesselOuterWallLogicalVolume->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		CryoSideTubeOut_log	->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		CryoSideTubeIn_log	->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PorcupineLeft_log	->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PorcupineRight_log	->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PTRpipe_log			->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PTRangle_log		->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PTRflange_log		->SetVisAttributes(pCryostatVesselOuterWallVisAtt);
		PTR_log				->SetVisAttributes(pCryostatVesselOuterWallVisAtt);

	G4VisAttributes *pCryostatVesselOuterWallDomeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatVesselOuterWallDomeVisAtt->SetVisibility(true);	//f
	m_pCryostatVesselOuterWallDomeLogicalVolume			->SetVisAttributes(pCryostatVesselOuterWallDomeVisAtt);

	m_pVacuumCryostatVesselLogicalVolume				->SetVisAttributes(G4VisAttributes::Invisible);

	G4VisAttributes *pCryostatLidFlangeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidFlangeVisAtt->SetVisibility(true);
	m_pCryostatLidFlangeLogicalVolume					->SetVisAttributes(pCryostatLidFlangeVisAtt);

	G4VisAttributes *pCryostatLidInnerWallVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidInnerWallVisAtt->SetVisibility(true);
	m_pCryostatLidInnerWallLogicalVolume				->SetVisAttributes(pCryostatLidInnerWallVisAtt);

	G4VisAttributes *pCryostatLidInnerWallDomeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidInnerWallDomeVisAtt->SetVisibility(true);
	m_pCryostatLidInnerWallDomeLogicalVolume			->SetVisAttributes(pCryostatLidInnerWallDomeVisAtt);

	G4VisAttributes *pCryostatLidOuterWallVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidOuterWallVisAtt->SetVisibility(true);
	m_pCryostatLidOuterWallLogicalVolume				->SetVisAttributes(pCryostatLidOuterWallVisAtt);

	G4VisAttributes *pCryostatLidOuterWallDomeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidOuterWallDomeVisAtt->SetVisibility(true);
	m_pCryostatLidOuterWallDomeLogicalVolume			->SetVisAttributes(pCryostatLidOuterWallVisAtt);

	G4VisAttributes *pCryostatLidInTubeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidInTubeVisAtt->SetVisibility(true);
	m_pCryostatLidInTubeLogicalVolume					->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeInsideVetoGXeLogicalVolume		->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeInsideVetoLXeLogicalVolume		->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeFlangeLogicalVolume				->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeElbowLogicalVolume				->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeTopLogicalVolume				->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeTopFlangeLogicalVolume			->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeTopToShieldLogicalVolume		->SetVisAttributes(pCryostatLidInTubeVisAtt);
	m_pCryostatLidInTubeTopToShieldOutsideLogicalVolume	->SetVisAttributes(pCryostatLidInTubeVisAtt);

	G4VisAttributes *pCryostatLidInTubeJacketVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidInTubeJacketVisAtt->SetVisibility(true);
	m_pCryostatLidInTubeJacketLogicalVolume				->SetVisAttributes(pCryostatLidInTubeJacketVisAtt);
	m_pCryostatLidInTubeJacketCoverLogicalVolume		->SetVisAttributes(pCryostatLidInTubeJacketVisAtt);

	G4VisAttributes *pCryostatLidInTubeOuterJacketVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidInTubeOuterJacketVisAtt->SetVisibility(true);
	m_pCryostatLidInTubeOuterJacketLogicalVolume		->SetVisAttributes(pCryostatLidInTubeOuterJacketVisAtt);
	m_pCryostatLidInTubeOuterJacketCoverLogicalVolume	->SetVisAttributes(pCryostatLidInTubeOuterJacketVisAtt);

	G4VisAttributes *pCryostatLidOutTubeVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidOutTubeVisAtt->SetVisibility(true);
	m_pCryostatLidOutTubeLogicalVolume					->SetVisAttributes(pCryostatLidOutTubeVisAtt);
	m_pCryostatLidOutTubeElbowLogicalVolume				->SetVisAttributes(pCryostatLidOutTubeVisAtt);
	m_pCryostatLidOutTubeTopLogicalVolume				->SetVisAttributes(pCryostatLidOutTubeVisAtt);
	m_pCryostatLidOutTubeTopFlangeLogicalVolume			->SetVisAttributes(pCryostatLidOutTubeVisAtt);
	m_pCryostatLidOutTubeTopToShieldLogicalVolume		->SetVisAttributes(pCryostatLidOutTubeVisAtt);
	m_pCryostatLidOutTubeTopToShieldOutsideLogicalVolume->SetVisAttributes(pCryostatLidOutTubeVisAtt);

	G4VisAttributes *pCryostatLidOutTubeJacketVisAtt = new G4VisAttributes(hSteelColour);
	pCryostatLidOutTubeJacketVisAtt->SetVisibility(true);
	m_pCryostatLidOutTubeJacketLogicalVolume		->SetVisAttributes(pCryostatLidOutTubeJacketVisAtt);
	m_pCryostatLidOutTubeJacketCoverLogicalVolume	->SetVisAttributes(pCryostatLidOutTubeJacketVisAtt);
//    G4VisAttributes *pVacuumCryostatLidVisAtt = new G4VisAttributes(hVacuumColor);
//    pVacuumCryostatLidVisAtt->SetVisibility(false);
//    m_pVacuumCryostatLidLogicalVolume->SetVisAttributes(pVacuumCryostatLidVisAtt);
}

void
Xenon100DetectorConstruction::PrintGeometryInformation()
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< xenon >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeTotalMass = m_pLXeLogicalVolume->GetMass(false, false)/kg;
	G4cout
		<< "Total LXe Mass: " << dLXeTotalMass << " kg" << G4endl;

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< pmt arrays >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");

	G4cout
		<< "Pmts: " << G4endl
		<< "  NbTopPmts: " << iNbTopPmts << " NbBottomPmts: " << iNbBottomPmts << G4endl
		<< "  NbTopVetoPmts: " << iNbTopVetoPmts << " NbBottomVetoPmts: " << iNbBottomVetoPmts << G4endl
		<< "  Total: " << iNbTopPmts + iNbBottomPmts + iNbTopVetoPmts + iNbBottomVetoPmts << G4endl;

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cryostat >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	
	//=================================== vessel ====================================
	G4double dCryostatVesselMass = 0.;
	dCryostatVesselMass += m_pCryostatVesselFlangeLogicalVolume			->GetMass(false, false)/kg;
	dCryostatVesselMass += m_pCryostatVesselInnerWallLogicalVolume		->GetMass(false, false)/kg;
	dCryostatVesselMass += m_pCryostatVesselInnerWallDomeLogicalVolume	->GetMass(false, false)/kg;
	dCryostatVesselMass += m_pCryostatVesselOuterWallLogicalVolume		->GetMass(false, false)/kg;
	dCryostatVesselMass += m_pCryostatVesselOuterWallDomeLogicalVolume	->GetMass(false, false)/kg;
	//===================================== lid =====================================
	G4double dCryostatLidMass = 0.;
	dCryostatLidMass += m_pCryostatLidFlangeLogicalVolume					->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInnerWallLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInnerWallDomeLogicalVolume			->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOuterWallLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOuterWallDomeLogicalVolume			->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeLogicalVolume					->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeInsideVetoGXeLogicalVolume		->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeInsideVetoLXeLogicalVolume		->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeFlangeLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeJacketLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeJacketCoverLogicalVolume		->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeOuterJacketLogicalVolume		->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeOuterJacketCoverLogicalVolume	->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeElbowLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidInTubeTopLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOutTubeLogicalVolume					->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOutTubeJacketLogicalVolume			->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOutTubeJacketCoverLogicalVolume		->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOutTubeElbowLogicalVolume				->GetMass(false, false)/kg;
	dCryostatLidMass += m_pCryostatLidOutTubeTopLogicalVolume				->GetMass(false, false)/kg;
	//=================================== vessel =======================================
	G4double dBellMass = 0.;
	dBellMass += m_pBellLogicalVolume->GetMass(false, false)/kg;
	dBellMass += m_pBellPmtSupportRingLogicalVolume->GetMass(false, false)/kg;

	//==================================== total =======================================
	const G4double dCryostatTotalMass = dCryostatVesselMass+dCryostatLidMass+dBellMass;

	G4cout
		<< "Cryostat Mass:" << G4endl
		<< "  CryostatVessel: " << dCryostatVesselMass << " kg" << G4endl
		<< "  CryostatLid:    " << dCryostatLidMass << " kg" << G4endl
		<< "  Bell:           " << dBellMass << " kg" << G4endl
		<< "  TOTAL:          " << dCryostatTotalMass << " kg" << G4endl;
		
	const G4double dCryostatVesselInnerWallMass		= m_pCryostatVesselInnerWallLogicalVolume		->GetMass(false, false)/kg;
	const G4double dCryostatVesselInnerWallDomeMass	= m_pCryostatVesselInnerWallDomeLogicalVolume	->GetMass(false, false)/kg;
	const G4double dCryostatVesselOuterWallMass		= m_pCryostatVesselOuterWallLogicalVolume		->GetMass(false, false)/kg;
	const G4double dCryostatVesselOuterWallDomeMass	= m_pCryostatVesselOuterWallDomeLogicalVolume	->GetMass(false, false)/kg;
	const G4double dCryostatVesselFlangeMass		= m_pCryostatVesselFlangeLogicalVolume			->GetMass(false, false)/kg;

	
	G4cout
		<< "   Vessel Inner Wall      = " << dCryostatVesselInnerWallMass		<< " kg" << G4endl
		<< "   Vessel Inner Wall Dome = " << dCryostatVesselInnerWallDomeMass	<< " kg" << G4endl
		<< "   Vessel Outer Wall      = " << dCryostatVesselOuterWallMass		<< " kg" << G4endl
		<< "   Vessel Outer Wall Dome = " << dCryostatVesselOuterWallDomeMass	<< " kg" << G4endl
		<< "   Vessel Flange          = "	<< dCryostatVesselFlangeMass		<< " kg" << G4endl;
	//==================================================================================
	
	G4double dTeflonTotalMass = 0.;
	dTeflonTotalMass += m_pTeflonSideVetoLiningLogicalVolume	->GetMass(false, false)/kg;
	dTeflonTotalMass += m_pTeflonPanelLogicalVolume				->GetMass(false, false)/kg;
	dTeflonTotalMass += m_pTopPmtTeflonHolderLogicalVolume		->GetMass(false, false)/kg;
	dTeflonTotalMass += m_pTeflonRodsLogicalVolume				->GetMass(false, false)/kg;
	dTeflonTotalMass += m_pLowerTeflonPanelLogicalVolume		->GetMass(false, false)/kg;
	
	G4cout
		<<" "<< G4endl
		<<"Teflon Mass: " << dTeflonTotalMass << " kg" << G4endl
		<<" "<< G4endl;

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< shield >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    const G4double dShieldPolishLeadMass = m_pShieldPolishLeadLogicalVolume		->GetMass(false, false)/kg;
    const G4double dShieldFrenchLeadMass = m_pShieldFrenchLeadLogicalVolume		->GetMass(false, false)/kg;
    //const G4double dShieldInnerLeadMass = m_pShieldInnerLeadLogicalVolume		->GetMass(false, false)/kg;
    const G4double dShieldPolyethyleneMass = m_pShieldPolyethyleneLogicalVolume	->GetMass(false, false)/kg;
	const G4double dShieldCopperMass = m_pShieldCopperLogicalVolume				->GetMass(false, false)/kg;
	const G4double dShieldCavityMass = m_pShieldCavityLogicalVolume				->GetMass(false, false)/kg;
	
    G4cout
        << "Shield:" << G4endl
        << "  ShieldPolishLead: " << dShieldPolishLeadMass << " kg" << G4endl
        << "  ShieldFrenchLead: " << dShieldFrenchLeadMass << " kg" << G4endl
       // << "  ShieldInnerLead (French, 3mm thick) : " << dShieldInnerLeadMass << " kg" << G4endl
        << "  ShieldPolyethylene: " << dShieldPolyethyleneMass << " kg" << G4endl
        << "  ShieldCopper: " << dShieldCopperMass << " kg" << G4endl
        << "  ShieldCavity mass: " << dShieldCavityMass << " kg" << G4endl;

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< copper >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	G4double dCopperTotalMass = 0.;
	dCopperTotalMass += m_pTopVetoAngleLogicalVolume		->GetMass(false, false)/kg;
	dCopperTotalMass += m_pUpperSideVetoAngleLogicalVolume	->GetMass(false, false)/kg;
	dCopperTotalMass += m_pTopPlateLogicalVolume			->GetMass(false, false)/kg;
	dCopperTotalMass += m_pBottomPlateLogicalVolume			->GetMass(false, false)/kg;
	dCopperTotalMass += m_pUpperBasePlateLogicalVolume		->GetMass(false, false)/kg;
	dCopperTotalMass += m_pLowerBasePlateLogicalVolume		->GetMass(false, false)/kg;
	dCopperTotalMass += m_pBottomPmtPlateLogicalVolume		->GetMass(false, false)/kg;
	dCopperTotalMass += m_pBottomVetoAngleLogicalVolume		->GetMass(false, false)/kg;
	dCopperTotalMass += m_pLowerSideVetoAngleLogicalVolume	->GetMass(false, false)/kg;	
	G4cout
		<<" "<< G4endl
		<<"Copper Mass: " << dCopperTotalMass << " kg" << G4endl
		<<" "<< G4endl;
	
	G4double dSBarsTotalMass = 0.;
	dSBarsTotalMass += m_pSBarSideLeftLogicalVolume		->GetMass(false, false)/kg;
	dSBarsTotalMass += m_pSBarSideRightLogicalVolume	->GetMass(false, false)/kg;
	dSBarsTotalMass += m_pSBarFrontLogicalVolume		->GetMass(false, false)/kg;	
	G4cout
		<< G4endl
		<<"Support Bars Mass: "<< dSBarsTotalMass << " kg" << G4endl;

	G4double dResistorChainMass = m_pResistorChainLogicalVolume ->GetMass(false, false)/kg;
	G4cout
		<<"Resistors Mass: "<< dResistorChainMass << " kg" << G4endl;

	//G4double dCableMass = m_pCableLogicalVolume ->GetMass(false, false)/kg;
	//G4cout
	//	<<"Cable Mass: "<< dCableMass << " kg" << G4endl;

	G4double dCathodeGridMeshMass = m_pCathodeGridMeshLogicalVolume ->GetMass(false, false)/kg;
	G4double dCathodeGridRingMass = m_pCathodeGridRingLogicalVolume ->GetMass(false, false)/kg;
	G4cout
		<< G4endl
		<<"Cathode Grid Mesh Mass: "<< dCathodeGridMeshMass << G4endl
		<<"Cathode Grid Ring Mass: "<< dCathodeGridRingMass << G4endl;
	
	G4double dTopGridsMass = 0.;
	dTopGridsMass 	+= m_pTopGridRingLogicalVolume 		->GetMass(false, false)/kg;
	dTopGridsMass 	+= m_pAnodeGridRingLogicalVolume 	->GetMass(false, false)/kg;
	dTopGridsMass 	+= m_pBottomGridRingLogicalVolume	->GetMass(false, false)/kg;
	G4cout
		<< G4endl
		<<"Top Grids Mass: "<< dTopGridsMass << G4endl;

	G4double dLeadBrickMass = m_pLeadBrickLogicalVolume ->GetMass(false, false)/kg;
	G4cout
		<< G4endl
		<<"Lead Brick Mass: "<< dLeadBrickMass << G4endl;

	G4double dBrickHolderMass = 0.;
	dBrickHolderMass 	+= m_pLongHolderLogicalVolume 	->GetMass(false, false)/kg;
	dBrickHolderMass 	+= m_pShortHolder1LogicalVolume ->GetMass(false, false)/kg;
	dBrickHolderMass 	+= m_pShortHolder2LogicalVolume	->GetMass(false, false)/kg;
	dBrickHolderMass 	+= m_pShortHolder3LogicalVolume	->GetMass(false, false)/kg;
	G4cout
		<< G4endl
		<<"Copper Holder Mass (for Lead Brick): "<< dBrickHolderMass << G4endl;

	//G4double dWireMass = m_pWireLogicalVolume ->GetMass(false, false)/g;
	//G4cout
	//	<< G4endl
	//	<<"Wire Mass: "<< dWireMass <<" g"<< G4endl;
	
}


////////////////////////////////////////////////////////////////////////////////////
G4ThreeVector
Xenon100DetectorConstruction::GetPmtPosition(G4int iPmtNb, PmtPart ePmtPart)
{
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");

	G4ThreeVector hPos;

	if(iPmtNb < iNbTopPmts)
		hPos = GetPmtPositionTopArray(iPmtNb, ePmtPart);
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts)
		hPos = GetPmtPositionBottomArray(iPmtNb, ePmtPart);
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts)
		hPos = GetPmtPositionTopVetoArray(iPmtNb, ePmtPart);
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts)
		hPos = GetPmtPositionBottomVetoArray(iPmtNb, ePmtPart);

	return hPos;
}

G4ThreeVector
Xenon100DetectorConstruction::GetPmtPositionTopArray(G4int iPmtNb, PmtPart ePmtPart)
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bottom array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dGXeHeight = GetGeometryParameter("GXeHeight");

	const G4double dTopPmtTeflonHolderCutDepth = GetGeometryParameter("TopPmtTeflonHolderCutDepth");
	const G4double dTopPmtTeflonHolderToBell = GetGeometryParameter("TopPmtTeflonHolderToBell");

	const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
	const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
	const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

	const G4double dPmtToPmtBase = GetGeometryParameter("PmtToPmtBase");

	const G4double dPmtSixthRingRadius = GetGeometryParameter("PmtSixthRingRadius");
	const G4double dPmtFifthRingRadius = GetGeometryParameter("PmtFifthRingRadius");
	const G4double dPmtFourthRingRadius = GetGeometryParameter("PmtFourthRingRadius");
	const G4double dPmtThirdRingRadius = GetGeometryParameter("PmtThirdRingRadius");

	vector<G4int> hTopPmtsPerRing;
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsSixthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFifthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFourthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsThirdRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsSecondRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFirstRing"));

	//================================ bottom array ================================-
	const G4double dGXeHalfZ = 0.5*dGXeHeight;

	const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
	const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
	const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

	const G4double dBottomPmtWindowOffsetZ = dGXeHalfZ-dTopPmtTeflonHolderToBell-dTopPmtTeflonHolderCutDepth+dPmtWindowHalfZ;
	const G4double dBottomPmtCasingOffsetZ = dBottomPmtWindowOffsetZ+dPmtWindowHalfZ+dPmtCasingHalfZ;
	const G4double dBottomPmtBaseOffsetZ = dBottomPmtCasingOffsetZ+dPmtCasingHalfZ+dPmtToPmtBase+dPmtBaseHalfZ;

	G4int iRing = 6;
	G4int iTotal = hTopPmtsPerRing[0];

	while(iPmtNb > iTotal-1)
		iTotal += hTopPmtsPerRing[6-(--iRing)];

	G4int iOffsetPmtNb = iTotal-hTopPmtsPerRing[6-iRing];

	const G4double dPmtDistance = 30.*mm;

	G4ThreeVector hPos;

	switch(iRing)
	{
		case 6:
			hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iOffsetPmtNb, hTopPmtsPerRing[6-iRing], dPmtSixthRingRadius);
			break;

		case 5:
			hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iOffsetPmtNb, hTopPmtsPerRing[6-iRing], dPmtFifthRingRadius);
			break;

		case 4:
			hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iOffsetPmtNb, hTopPmtsPerRing[6-iRing], dPmtFourthRingRadius);
			break;

		case 3:
			hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iOffsetPmtNb, hTopPmtsPerRing[6-iRing], dPmtThirdRingRadius);
			break;

		case 2:
			switch(iPmtNb)
			{
				case 88:
				case 92:
					hPos.setX((iPmtNb == 88)?(-1.5*dPmtDistance):(1.5*dPmtDistance));
					hPos.setY(0.*mm);
					break;

				case 90:
				case 94:
					hPos.setX(0.*mm);
					hPos.setY((iPmtNb == 90)?(1.5*dPmtDistance):(-1.5*dPmtDistance));
					break;

				case 95:
				case 89:
					hPos.setX(-dPmtDistance);
					hPos.setY((iPmtNb == 89)?(dPmtDistance):(-dPmtDistance));
					break;

				case 93:
				case 91:
					hPos.setX(dPmtDistance);
					hPos.setY((iPmtNb == 91)?(dPmtDistance):(-dPmtDistance));
					break;
			}
			break;

		case 1:
			hPos.setX((iPmtNb == 96)?(-0.5*dPmtDistance):(0.5*dPmtDistance));
			hPos.setY(0.*mm);
			break;
	}

	switch(ePmtPart)
	{
		case PMT_WINDOW:
			hPos.setZ(dBottomPmtWindowOffsetZ);
			break;

		case PMT_CASING:
			hPos.setZ(dBottomPmtCasingOffsetZ);
			break;

		case PMT_BASE:
			hPos.setZ(dBottomPmtBaseOffsetZ);
			break;
	}

	return hPos;
}

G4ThreeVector
Xenon100DetectorConstruction::GetPmtPositionBottomArray(G4int iPmtNb, PmtPart ePmtPart)
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bottom array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");

	const G4double dBellHeight = GetGeometryParameter("BellHeight");
	const G4double dTopPlateThickness = GetGeometryParameter("TopPlateThickness");
	const G4double dTeflonPanelHeight = GetGeometryParameter("TeflonPanelHeight");
	const G4double dBottomPlateThickness = GetGeometryParameter("BottomPlateThickness");
	const G4double dUpperBasePlateThickness = GetGeometryParameter("UpperBasePlateThickness");
	const G4double dLowerTeflonPanelHeight = GetGeometryParameter("LowerTeflonPanelHeight");
	const G4double dLowerBasePlateThickness = GetGeometryParameter("LowerBasePlateThickness");

	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");

	const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
	const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
	const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

	const G4double dPmtToPmtBase= GetGeometryParameter("PmtToPmtBase");
	const G4double dPmtBaseSpacerHeight= GetGeometryParameter("PmtBaseSpacerHeight");

	vector<G4int> hBottomPmtsPerRow;
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsFirstRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsSecondRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsThirdRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsFourthRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsFifthRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsSixthRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsSeventhRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsEighthRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsNinthRow"));
	hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsTenthRow"));

	//================================ bottom array ================================-
	const G4double dLXeHalfZ = 0.5*dLXeHeight;

	const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
	const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
	const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

	const G4double dBottomPmtBaseOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight-dTopPlateThickness-dTeflonPanelHeight-dBottomPlateThickness-dUpperBasePlateThickness-dLowerTeflonPanelHeight-dLowerBasePlateThickness+dPmtBaseSpacerHeight+dPmtBaseHalfZ;
	const G4double dBottomPmtCasingOffsetZ = dBottomPmtBaseOffsetZ+dPmtToPmtBase+dPmtBaseHalfZ+dPmtCasingHalfZ;
	const G4double dBottomPmtWindowOffsetZ = dBottomPmtCasingOffsetZ+dPmtCasingHalfZ+dPmtWindowHalfZ;

	G4ThreeVector hPos = ComputeXYPmtPositionForGridPattern(iPmtNb-iNbTopPmts, hBottomPmtsPerRow);

	switch(ePmtPart)
	{
		case PMT_WINDOW:
			hPos.setZ(dBottomPmtWindowOffsetZ);
			break;

		case PMT_CASING:
			hPos.setZ(dBottomPmtCasingOffsetZ);
			break;

		case PMT_BASE:
			hPos.setZ(dBottomPmtBaseOffsetZ);
			break;
	}

	return hPos;
}

G4ThreeVector
Xenon100DetectorConstruction::GetPmtPositionTopVetoArray(G4int iPmtNb, PmtPart ePmtPart)
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< top veto array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");

	const G4double dUpperSideVetoAngleThickness = GetGeometryParameter("UpperSideVetoAngleThickness");
	const G4double dUpperSideVetoAngleHeight = GetGeometryParameter("UpperSideVetoAngleHeight");
	const G4double dUpperSideVetoAngleLength = GetGeometryParameter("UpperSideVetoAngleLength");

	const G4double dTopVetoAngleThickness = GetGeometryParameter("TopVetoAngleThickness");
	const G4double dTopVetoAngleLength = GetGeometryParameter("TopVetoAngleLength");
	const G4double dTopVetoAngleLowerHeight = GetGeometryParameter("TopVetoAngleLowerHeight");
	const G4double dTopVetoAngleUpperHeight = GetGeometryParameter("TopVetoAngleUpperHeight");

	const G4double dBellPmtSupportRingRadius = GetGeometryParameter("BellPmtSupportRingRadius");

	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	
	const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
	const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
	const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

	const G4double dPmtToPmtBase = GetGeometryParameter("PmtToPmtBase");
	const G4double dPmtBaseSpacerHeight = GetGeometryParameter("PmtBaseSpacerHeight");
	
	//=============================== top veto array ================================
	const G4double dLXeHalfZ = 0.5*dLXeHeight;

	const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
	const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
	const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

	const G4double dUpperSideVetoPmtOffsetX = dBellPmtSupportRingRadius+0.5*(dUpperSideVetoAngleLength+dUpperSideVetoAngleThickness);

	const G4double dUpperSideVetoPmtBaseOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))+dUpperSideVetoAngleHeight-dUpperSideVetoAngleThickness-dPmtBaseHalfZ-dPmtBaseSpacerHeight;
	const G4double dUpperSideVetoPmtCasingOffsetZ = dUpperSideVetoPmtBaseOffsetZ-dPmtBaseHalfZ-dPmtToPmtBase-dPmtCasingHalfZ;
	const G4double dUpperSideVetoPmtWindowOffsetZ = dUpperSideVetoPmtCasingOffsetZ-dPmtCasingHalfZ-dPmtWindowHalfZ;

	const G4double dTopVetoPmtOffsetX = dBellPmtSupportRingRadius+dTopVetoAngleLength-dTopVetoAngleThickness;

	const G4double dTopVetoPmtBaseOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))+dTopVetoAngleLowerHeight+dTopVetoAngleUpperHeight-dTopVetoAngleThickness-0.5*(dTopVetoAngleUpperHeight-dTopVetoAngleThickness);

	const G4double dTopVetoPmtBaseOffsetX = dTopVetoPmtOffsetX-dPmtBaseSpacerHeight-dPmtBaseHalfZ;
	const G4double dTopVetoPmtCasingOffsetX = dTopVetoPmtBaseOffsetX-dPmtBaseHalfZ-dPmtToPmtBase-dPmtCasingHalfZ;
	const G4double dTopVetoPmtWindowOffsetX = dTopVetoPmtCasingOffsetX-dPmtCasingHalfZ-dPmtWindowHalfZ;

	G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts, iNbTopVetoPmts, dUpperSideVetoPmtOffsetX);
	G4ThreeVector hPosBase = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts, iNbTopVetoPmts, dTopVetoPmtBaseOffsetX);
	G4ThreeVector hPosCasing = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts, iNbTopVetoPmts, dTopVetoPmtCasingOffsetX);
	G4ThreeVector hPosWindow = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts, iNbTopVetoPmts, dTopVetoPmtWindowOffsetX);

	if(pow(-1, iPmtNb) == -1)
	{
		switch(ePmtPart)
		{
			case PMT_BASE:
				hPos.set(hPosBase.x(), hPosBase.y(), dTopVetoPmtBaseOffsetZ);
				break;

			case PMT_CASING:
				hPos.set(hPosCasing.x(), hPosCasing.y(), dTopVetoPmtBaseOffsetZ);
				break;

			case PMT_WINDOW:
				hPos.set(hPosWindow.x(), hPosWindow.y(), dTopVetoPmtBaseOffsetZ);
				break;
		}
	}
	else
	{
		switch(ePmtPart)
		{
			case PMT_BASE:
				hPos.setZ(dUpperSideVetoPmtBaseOffsetZ);
				break;

			case PMT_CASING:
				hPos.setZ(dUpperSideVetoPmtCasingOffsetZ);
				break;

			case PMT_WINDOW:
				hPos.setZ(dUpperSideVetoPmtWindowOffsetZ);
				break;
		}
	}

	return hPos;
}

G4ThreeVector
Xenon100DetectorConstruction::GetPmtPositionBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart)
{
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bottom veto array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

	const G4double dLidInnerWallDomeHeight = GetGeometryParameter("CryostatLidInnerWallDomeHeight");
	const G4double dLidInTubeHeightInside = GetGeometryParameter("CryostatLidInTubeHeightInside");	
	const G4double dLidInTubeOffsetDifference = GetGeometryParameter("CryostatLidInTubeOffsetDifference");

	const G4double dBellHeight = GetGeometryParameter("BellHeight");
	const G4double dTopPlateThickness = GetGeometryParameter("TopPlateThickness");
	const G4double dTeflonPanelHeight = GetGeometryParameter("TeflonPanelHeight");
	const G4double dBottomPlateThickness = GetGeometryParameter("BottomPlateThickness");
	const G4double dUpperBasePlateThickness = GetGeometryParameter("UpperBasePlateThickness");
	const G4double dLowerTeflonPanelHeight = GetGeometryParameter("LowerTeflonPanelHeight");
	const G4double dLowerBasePlateThickness = GetGeometryParameter("LowerBasePlateThickness");

	const G4double dBottomVetoAngleThickness = GetGeometryParameter("BottomVetoAngleThickness");
	const G4double dBottomVetoAngleHeight = GetGeometryParameter("BottomVetoAngleHeight");
	const G4double dBottomVetoAngleLength = GetGeometryParameter("BottomVetoAngleLength");

	const G4double dLowerSideVetoAngleThickness = GetGeometryParameter("LowerSideVetoAngleThickness");
	const G4double dLowerSideVetoAngleHeight = GetGeometryParameter("LowerSideVetoAngleHeight");
	const G4double dLowerSideVetoAngleShortLength = GetGeometryParameter("LowerSideVetoAngleShortLength");
	const G4double dLowerSideVetoAngleLongLength = GetGeometryParameter("LowerSideVetoAngleLongLength");

	const G4double dLowerBasePlateSupportRingRadius = GetGeometryParameter("LowerBasePlateSupportRingRadius");

	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");
	
	const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
	const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
	const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

	const G4double dPmtToPmtBase = GetGeometryParameter("PmtToPmtBase");
	const G4double dPmtBaseSpacerHeight = GetGeometryParameter("PmtBaseSpacerHeight");
	
	//============================= bottom veto array ===============================
	const G4double dLXeHalfZ = 0.5*dLXeHeight;

	const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
	const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
	const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

	const G4double dLowerSideVetoPmtOffsetX = dLowerBasePlateSupportRingRadius+dLowerSideVetoAngleShortLength+0.5*(dLowerSideVetoAngleLongLength-dLowerSideVetoAngleThickness);

	const G4double dLowerSideVetoPmtBaseOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight-dTopPlateThickness-dTeflonPanelHeight-dBottomPlateThickness-dUpperBasePlateThickness-dLowerTeflonPanelHeight-dLowerBasePlateThickness-dLowerSideVetoAngleHeight+dLowerSideVetoAngleThickness+dPmtBaseSpacerHeight+dPmtBaseHalfZ;

	const G4double dLowerSideVetoPmtCasingOffsetZ = dLowerSideVetoPmtBaseOffsetZ+dPmtBaseHalfZ+dPmtToPmtBase+dPmtCasingHalfZ;
	const G4double dLowerSideVetoPmtWindowOffsetZ = dLowerSideVetoPmtCasingOffsetZ+dPmtCasingHalfZ+dPmtWindowHalfZ;

	const G4double dBottomVetoPmtOffsetX = dLowerBasePlateSupportRingRadius+dBottomVetoAngleLength-dBottomVetoAngleThickness;

	const G4double dBottomVetoPmtBaseOffsetZ = dLXeHalfZ-(dLidInTubeHeightInside+(dLidInnerWallDomeHeight-dLidInTubeOffsetDifference))-dBellHeight-dTopPlateThickness-dTeflonPanelHeight-dBottomPlateThickness-dUpperBasePlateThickness-dLowerTeflonPanelHeight-dLowerBasePlateThickness-0.5*(dBottomVetoAngleHeight+dBottomVetoAngleThickness);

	const G4double dBottomVetoPmtBaseOffsetX = dBottomVetoPmtOffsetX-dPmtBaseSpacerHeight-dPmtBaseHalfZ;
	const G4double dBottomVetoPmtCasingOffsetX = dBottomVetoPmtBaseOffsetX-dPmtBaseHalfZ-dPmtToPmtBase-dPmtCasingHalfZ;
	const G4double dBottomVetoPmtWindowOffsetX = dBottomVetoPmtCasingOffsetX-dPmtCasingHalfZ-dPmtWindowHalfZ;

	G4ThreeVector hPos = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts-iNbBottomVetoPmts, iNbTopVetoPmts, dLowerSideVetoPmtOffsetX);
	G4ThreeVector hPosBase = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts-iNbBottomVetoPmts, iNbTopVetoPmts, dBottomVetoPmtBaseOffsetX);
	G4ThreeVector hPosCasing = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts-iNbBottomVetoPmts, iNbTopVetoPmts, dBottomVetoPmtCasingOffsetX);
	G4ThreeVector hPosWindow = ComputeXYPmtPositionForRingPattern(iPmtNb-iNbTopPmts-iNbBottomPmts-iNbBottomVetoPmts, iNbTopVetoPmts, dBottomVetoPmtWindowOffsetX);

	if(pow(-1, iPmtNb) == -1)
	{
		switch(ePmtPart)
		{
			case PMT_WINDOW:
				hPos.set(hPosWindow.x(), hPosWindow.y(), dBottomVetoPmtBaseOffsetZ);
				break;

			case PMT_CASING:
				hPos.set(hPosCasing.x(), hPosCasing.y(), dBottomVetoPmtBaseOffsetZ);
				break;

			case PMT_BASE:
				hPos.set(hPosBase.x(), hPosBase.y(), dBottomVetoPmtBaseOffsetZ);
				break;
		}
	}
	else
	{
		switch(ePmtPart)
		{
			case PMT_WINDOW:
				hPos.setZ(dLowerSideVetoPmtWindowOffsetZ);
				break;

			case PMT_CASING:
				hPos.setZ(dLowerSideVetoPmtCasingOffsetZ);
				break;

			case PMT_BASE:
				hPos.setZ(dLowerSideVetoPmtBaseOffsetZ);
				break;
		}
	}

	return hPos;
}

G4RotationMatrix *
Xenon100DetectorConstruction::GetPmtRotation(G4int iPmtNb, PmtPart ePmtPart)
{
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");

	G4RotationMatrix *pRotationMatrix;

	if(iPmtNb < iNbTopPmts)
		pRotationMatrix = GetPmtRotationTopArray(iPmtNb, ePmtPart);
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts)
		pRotationMatrix = m_pRotationX180;
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts)
		pRotationMatrix = GetPmtRotationTopVetoArray(iPmtNb, ePmtPart);
	else if(iPmtNb < iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts)
		pRotationMatrix = GetPmtRotationBottomVetoArray(iPmtNb, ePmtPart);

	return pRotationMatrix;
}

G4RotationMatrix *
Xenon100DetectorConstruction::GetPmtRotationTopArray(G4int iPmtNb, PmtPart ePmtPart)
{
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");

	vector<G4int> hTopPmtsPerRing;
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsSixthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFifthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFourthRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsThirdRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsSecondRing"));
	hTopPmtsPerRing.push_back((G4int) GetGeometryParameter("NbTopPmtsFirstRing"));

	G4int iRing = 6;
	G4int iTotal = hTopPmtsPerRing[0];

	while(iPmtNb > iTotal-1)
		iTotal += hTopPmtsPerRing[6-(--iRing)];

	G4int iOffsetPmtNb = iTotal-hTopPmtsPerRing[6-iRing];

	G4double dPmtAngularSpacing = 360.*deg/iNbTopPmts;

	switch(iRing)
	{
		case 6:
			dPmtAngularSpacing = 360.*deg/hTopPmtsPerRing[0];
			break;

		case 5:
			dPmtAngularSpacing = 360.*deg/hTopPmtsPerRing[1];
			break;

		case 4:
			dPmtAngularSpacing = 360.*deg/hTopPmtsPerRing[2];
			break;

		case 3:
			dPmtAngularSpacing = 360.*deg/hTopPmtsPerRing[3];
			break;

		case 2:
		case 1:
			dPmtAngularSpacing = 0.*deg;
			break;
	}

	G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();

	pRotationMatrix->rotateZ((iPmtNb-iOffsetPmtNb)*dPmtAngularSpacing);

	return pRotationMatrix;
}

G4RotationMatrix *
Xenon100DetectorConstruction::GetPmtRotationTopVetoArray(G4int iPmtNb, PmtPart ePmtPart)
{
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");

	G4double dPmtAngularSpacing = 360.*deg/iNbTopVetoPmts;

	G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();

	pRotationMatrix->rotateZ((iPmtNb-iNbTopPmts-iNbBottomPmts)*dPmtAngularSpacing-90.*deg);
	if(pow(-1, iPmtNb) == -1)
		pRotationMatrix->rotateX(90.*deg);

	return pRotationMatrix;
}

G4RotationMatrix *
Xenon100DetectorConstruction::GetPmtRotationBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart)
{
	const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");
	const G4int iNbTopVetoPmts = (G4int) GetGeometryParameter("NbTopVetoPmts");
	const G4int iNbBottomVetoPmts = (G4int) GetGeometryParameter("NbBottomVetoPmts");

	G4double dPmtAngularSpacing = 360.*deg/iNbBottomVetoPmts;

	G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();

	pRotationMatrix->rotateZ((iPmtNb-iNbTopPmts-iNbBottomPmts-iNbTopVetoPmts)*dPmtAngularSpacing-90.*deg);
	pRotationMatrix->rotateX(180.*deg);
	if(pow(-1, iPmtNb) == -1)
		pRotationMatrix->rotateX(-90.*deg);

	return pRotationMatrix;
}

vector<G4int>
Xenon100DetectorConstruction::ComputePmtPartition(G4double dRadius, G4double dHeightFraction)
{
	G4double dPmtWidth = Xenon100DetectorConstruction::GetGeometryParameter("PmtWidth");
	G4double dPmtSpacing = Xenon100DetectorConstruction::GetGeometryParameter("PmtSpacing");
	G4int iNbRows = (G4int) ceil(2*dRadius/(dPmtWidth+dPmtSpacing));

	vector<G4int> hPmtsPerRow(iNbRows);

	G4int iStartRow = (iNbRows%2)?(iNbRows/2):(iNbRows/2-1);
	G4int iCurrentRow = iStartRow;
	G4int iMirrorCurrentRow = iNbRows/2;

	hPmtsPerRow[iCurrentRow--] = iNbRows;
	hPmtsPerRow[iMirrorCurrentRow++] = iNbRows;

	while(iCurrentRow >= 0)
	{
		G4double dY = (iStartRow - iCurrentRow + dHeightFraction)*(dPmtWidth+dPmtSpacing);
		G4double dWidth = 2.*sqrt(dRadius*dRadius - dY*dY);
		G4int iNbPmtsInRow = (G4int) ceil(dWidth/(dPmtWidth+dPmtSpacing));

		hPmtsPerRow[iCurrentRow--] = iNbPmtsInRow;
		hPmtsPerRow[iMirrorCurrentRow++] = iNbPmtsInRow;
	}

	return(hPmtsPerRow);
}

G4ThreeVector
Xenon100DetectorConstruction::ComputeXYPmtPositionForGridPattern(G4int iPmtNb, const vector<G4int> &hPmtsPerRow)
{
	G4double dPmtWidth = Xenon100DetectorConstruction::GetGeometryParameter("PmtWidth");
	G4double dPmtSpacing = Xenon100DetectorConstruction::GetGeometryParameter("PmtSpacing");

	G4int iRow = 0, iColumn;
	G4int iTotal = hPmtsPerRow[0];

	while(iPmtNb > iTotal-1)
		iTotal += hPmtsPerRow[++iRow];
	
	iColumn = hPmtsPerRow[iRow]-(iTotal-iPmtNb);

	G4double dX = -(hPmtsPerRow[iRow]-1)*(dPmtWidth+dPmtSpacing)/2. + iColumn*(dPmtWidth+dPmtSpacing);
	G4double dY = (hPmtsPerRow.size()-1)*(dPmtWidth+dPmtSpacing)/2. - iRow*(dPmtWidth+dPmtSpacing);

	return G4ThreeVector(dX, dY, 0.);
}

G4ThreeVector
Xenon100DetectorConstruction::ComputeXYPmtPositionForRingPattern(G4int iPmtNb, G4int iNbPmts, G4double dRingRadius)
{
	G4double dPmtAngularSpacing = 2.*M_PI/iNbPmts;

	G4double dX = -dRingRadius*cos(iPmtNb*dPmtAngularSpacing);
	G4double dY = dRingRadius*sin(iPmtNb*dPmtAngularSpacing);

	return G4ThreeVector(dX, dY, 0.);
}

