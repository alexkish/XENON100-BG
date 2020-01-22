#ifndef __XENON10PDETECTORCONSTRUCTION_H__
#define __XENON10PDETECTORCONSTRUCTION_H__

#include <globals.hh>

#include <vector>
#include <map>

using std::vector;
using std::map;

class G4LogicalVolume;
class G4VPhysicalVolume;

#include <G4VUserDetectorConstruction.hh>

class Xenon100DetectorConstruction: public G4VUserDetectorConstruction
{
public:
	Xenon100DetectorConstruction();
	~Xenon100DetectorConstruction();

	G4VPhysicalVolume* Construct();

	static G4double GetGeometryParameter(const char *szParameter);

private:
	//void DefineColours();
	void DefineMaterials();
	void DefineGeometryParameters();

	void ConstructLaboratory();
	void ConstructShield();
	void ConstructXenon();
	void ConstructBell();
	void ConstructFieldCage();
	void ConstructTPC();
	void ConstructPmtSupports();
	void ConstructPmtArrays();
	void ConstructCryostat();

	void PrintGeometryInformation();

	typedef enum {PMT_CASING, PMT_WINDOW, PMT_BASE} PmtPart;

	G4ThreeVector GetPmtPosition(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionTopArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionBottomArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionTopVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotation(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationTopArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationTopVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart);

	vector<G4int> ComputePmtPartition(G4double dRadius, G4double dHeightFraction);
	G4ThreeVector ComputeXYPmtPositionForGridPattern(G4int iPmtNb, const vector<G4int> &hPmtsPerRow);
	G4ThreeVector ComputeXYPmtPositionForRingPattern(G4int iPmtNb, G4int iNbPmts, G4double dRingRadius);

private:
	// rotation matrices
	G4RotationMatrix *m_pRotationXPlus225;
	G4RotationMatrix *m_pRotationXMinus225;
	G4RotationMatrix *m_pRotationXPlus45;
	G4RotationMatrix *m_pRotationXMinus45;
	G4RotationMatrix *m_pRotationXPlus90;
	G4RotationMatrix *m_pRotationXMinus90;
	G4RotationMatrix *m_pRotationYPlus90;
	G4RotationMatrix *m_pRotationYMinus90;
	G4RotationMatrix *m_pRotationX180;
	G4RotationMatrix *m_pRotationZTeflonLedBlock;
	G4RotationMatrix *m_pRotationXPlus90YMinus90;
	// logical volumes
	G4LogicalVolume *m_pMotherLogicalVolume;

	G4LogicalVolume *m_pLabLogicalVolume;

	G4LogicalVolume *m_pShieldPolishLeadLogicalVolume;
	G4LogicalVolume *m_pShieldFrenchLeadLogicalVolume;
	//G4LogicalVolume *m_pShieldInnerLeadLogicalVolume;
	G4LogicalVolume *m_pShieldPolyethyleneLogicalVolume;
	G4LogicalVolume *m_pPolyethyleneSlabLogicalVolume;
	G4LogicalVolume *m_pShieldCopperLogicalVolume;
	G4LogicalVolume *m_pWaterTopLogicalVolume;
	G4LogicalVolume *m_pWaterBackLogicalVolume;
	G4LogicalVolume *m_pWaterLeftLogicalVolume;
	G4LogicalVolume *m_pWaterRightLogicalVolume;
	G4LogicalVolume *m_pShieldCavityLogicalVolume;
	
	G4LogicalVolume *m_pSBarSideLeftLogicalVolume;
	G4LogicalVolume *m_pSBarSideRightLogicalVolume;
	G4LogicalVolume *m_pSBarFrontLogicalVolume;
	
	G4LogicalVolume *m_pLongHolderLogicalVolume;
	G4LogicalVolume *m_pShortHolder1LogicalVolume;
	G4LogicalVolume *m_pShortHolder2LogicalVolume;
	G4LogicalVolume *m_pShortHolder3LogicalVolume;
	G4LogicalVolume *m_pLeadBrickLogicalVolume;
	G4LogicalVolume *m_pSourcePipeLogicalVolume;

	//G4LogicalVolume *m_pWireLogicalVolume;
	//G4LogicalVolume *m_pCollimatorLogicalVolume;
	G4LogicalVolume *m_pCalibrationSourceLogicalVolume;
	
	G4LogicalVolume *m_pLXeLogicalVolume;
	G4LogicalVolume *m_pGXeLogicalVolume;
	G4LogicalVolume *m_pVetoGXeLogicalVolume;
	G4LogicalVolume *m_pLidInTubeGXeLogicalVolume;

	G4LogicalVolume *m_pBellLogicalVolume;	
	G4LogicalVolume *m_pTopPmtTeflonHolderLogicalVolume;
	G4LogicalVolume *m_pBellPmtSupportRingLogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInLXeLogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInVetoGXe1LogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInVetoGXe2LogicalVolume;
	G4LogicalVolume *m_pTopVetoAngleLogicalVolume;
	G4LogicalVolume *m_pUpperSideVetoAngleLogicalVolume;

	G4LogicalVolume *m_pTopGridMeshLogicalVolume;
	G4LogicalVolume *m_pTopGridRingLogicalVolume;
	G4LogicalVolume *m_pAnodeGridMeshLogicalVolume;
	G4LogicalVolume *m_pAnodeGridRingLogicalVolume;
	G4LogicalVolume *m_pBottomGridMeshLogicalVolume;
	G4LogicalVolume *m_pBottomGridRingLogicalVolume;
	G4LogicalVolume *m_pCathodeGridMeshLogicalVolume;
	G4LogicalVolume *m_pCathodeGridRingLogicalVolume;
	G4LogicalVolume *m_pScreenMeshLogicalVolume;

	G4LogicalVolume *m_pGridHolderInLXeLogicalVolume;
	G4LogicalVolume *m_pGridHolderInGXeLogicalVolume;
	G4LogicalVolume *m_pCopperScreenLogicalVolume;

	G4LogicalVolume *m_pTopPlateLogicalVolume;	
	G4LogicalVolume *m_pTeflonPanelLogicalVolume;
	G4LogicalVolume *m_pTeflonRodsLogicalVolume;
	G4LogicalVolume *m_pBottomPlateLogicalVolume;	
	G4LogicalVolume *m_pUpperBasePlateLogicalVolume;	
	G4LogicalVolume *m_pLowerTeflonPanelLogicalVolume;
	G4LogicalVolume *m_pLowerBasePlateLogicalVolume;	
	G4LogicalVolume *m_pBottomPmtPlateLogicalVolume;	
	G4LogicalVolume *m_pBottomVetoAngleLogicalVolume;
	G4LogicalVolume *m_pLowerSideVetoAngleLogicalVolume;
	G4LogicalVolume *m_pTeflonSideVetoLiningLogicalVolume;

	G4LogicalVolume *m_pTeflonLedBlockLogicalVolume;
	G4LogicalVolume *m_pTeflonBottomReflectorTopLogicalVolume;
	G4LogicalVolume *m_pTeflonBottomReflectorBottomLogicalVolume;
	G4LogicalVolume *m_pBottomReflectorHolderLogicalVolume;

	G4LogicalVolume *m_pPmtWindowLogicalVolume;
	G4LogicalVolume *m_pPmtCasingLogicalVolume;
	G4LogicalVolume *m_pPmtInteriorLogicalVolume;
	G4LogicalVolume *m_pPmtPhotoCathodeLogicalVolume;
	G4LogicalVolume *m_pPmtBaseLogicalVolume;

	G4LogicalVolume *m_pCryostatVesselFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselInnerWallLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselInnerWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselOuterWallLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselOuterWallDomeLogicalVolume;
	G4LogicalVolume *m_pVacuumCryostatVesselLogicalVolume;

	G4LogicalVolume *m_pCryostatLidFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInnerWallLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInnerWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOuterWallLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOuterWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeInsideVetoGXeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeInsideVetoLXeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeElbowLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopToShieldLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopToShieldOutsideLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeJacketCoverLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeOuterJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeOuterJacketCoverLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeElbowLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopToShieldLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopToShieldOutsideLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeJacketCoverLogicalVolume;

	G4LogicalVolume *m_pCryostatLidLogicalVolume;
	G4LogicalVolume *m_pVacuumCryostatLidLogicalVolume;

	G4LogicalVolume *CryoSideTubeOut_log;
	G4LogicalVolume *CryoSideTubeIn_log;
	G4LogicalVolume *PorcupineLeft_log;
	G4LogicalVolume *PorcupineRight_log;
	G4LogicalVolume *PTRpipe_log;
	G4LogicalVolume *PTRangle_log;
	G4LogicalVolume *PTRflange_log;
	G4LogicalVolume *PTR_log;
	
	G4LogicalVolume *m_pResistorChainLogicalVolume;
	G4LogicalVolume *m_pCableLogicalVolume;
	
	
	
	// physical volumes
	G4VPhysicalVolume *m_pLabPhysicalVolume;

	G4VPhysicalVolume *m_pShieldPolishLeadPhysicalVolume;
	G4VPhysicalVolume *m_pShieldFrenchLeadPhysicalVolume;
	//G4VPhysicalVolume *m_pShieldInnerLeadPhysicalVolume;
	G4VPhysicalVolume *m_pShieldPolyethylenePhysicalVolume;
	G4VPhysicalVolume *m_pPolyethyleneSlabPhysicalVolume;
	G4VPhysicalVolume *m_pShieldCopperPhysicalVolume;
	G4VPhysicalVolume *m_pWaterTopPhysicalVolume;
	G4VPhysicalVolume *m_pWaterBackPhysicalVolume;
	G4VPhysicalVolume *m_pWaterLeftPhysicalVolume;
	G4VPhysicalVolume *m_pWaterRightPhysicalVolume;
	G4VPhysicalVolume *m_pShieldCavityPhysicalVolume;
	
	G4VPhysicalVolume *m_pSBarSideLeftPhysicalVolume;
	G4VPhysicalVolume *m_pSBarSideRightPhysicalVolume;
	G4VPhysicalVolume *m_pSBarFrontPhysicalVolume;

	G4VPhysicalVolume *m_pLongHolderPhysicalVolume;
	G4VPhysicalVolume *m_pShortHolder1PhysicalVolume;
	G4VPhysicalVolume *m_pShortHolder2PhysicalVolume;
	G4VPhysicalVolume *m_pShortHolder3PhysicalVolume;
	G4VPhysicalVolume *m_pLeadBrickPhysicalVolume;
	G4VPhysicalVolume *m_pSourcePipePhysicalVolume;
	
	//G4VPhysicalVolume *m_pWirePhysicalVolume;
	//G4VPhysicalVolume *m_pCollimatorPhysicalVolume;

	G4VPhysicalVolume *m_pCalibrationSourcePhysicalVolume;

	G4VPhysicalVolume *m_pLXePhysicalVolume;
	G4VPhysicalVolume *m_pGXePhysicalVolume;
	G4VPhysicalVolume *m_pVetoGXePhysicalVolume;
	G4VPhysicalVolume *m_pLidInTubeGXePhysicalVolume;

	G4VPhysicalVolume *m_pBellPhysicalVolume;
	G4VPhysicalVolume *m_pTopPmtTeflonHolderPhysicalVolume;
	G4VPhysicalVolume *m_pBellPmtSupportRingPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hTopVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hUpperSideVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hBellSupportCylinderInLXePhysicalVolumes;
    G4VPhysicalVolume *m_pBellSupportCylinderInVetoGXe1PhysicalVolume;
    G4VPhysicalVolume *m_pBellSupportCylinderInVetoGXe2PhysicalVolume;
	
	G4VPhysicalVolume *m_pTopGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pTopGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pBottomGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pBottomGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridRingPhysicalVolume;	
	G4VPhysicalVolume *m_pScreenMeshPhysicalVolume;

	G4VPhysicalVolume *m_pGridHolderInLXePhysicalVolume;
	G4VPhysicalVolume *m_pGridHolderInGXePhysicalVolume;
	G4VPhysicalVolume *m_pCopperScreenPhysicalVolume;

	G4VPhysicalVolume *m_pTopPlatePhysicalVolume;
	G4VPhysicalVolume *m_pTeflonPanelPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hTeflonRodsPhysicalVolumes;
	G4VPhysicalVolume *m_pBottomPlatePhysicalVolume;
	G4VPhysicalVolume *m_pUpperBasePlatePhysicalVolume;
	G4VPhysicalVolume *m_pLowerTeflonPanelPhysicalVolume;
	G4VPhysicalVolume *m_pLowerBasePlatePhysicalVolume;
	G4VPhysicalVolume *m_pBottomPmtPlatePhysicalVolume;
	vector<G4VPhysicalVolume *> m_hBottomVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hLowerSideVetoAnglePhysicalVolumes;
	G4VPhysicalVolume *m_pTeflonSideVetoLiningPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hTeflonLedBlockPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hBottomReflectorHolderPhysicalVolumes;
	G4VPhysicalVolume *m_pTeflonBottomReflectorTopPhysicalVolume;
	G4VPhysicalVolume *m_pTeflonBottomReflectorBottomPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hPmtWindowPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hPmtCasingPhysicalVolumes;
	G4VPhysicalVolume *m_pPmtInteriorPhysicalVolume;
	G4VPhysicalVolume *m_pPmtPhotoCathodePhysicalVolume;
	vector<G4VPhysicalVolume *> m_hPmtBasePhysicalVolumes;

	G4VPhysicalVolume *m_pCryostatVesselFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselInnerWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselInnerWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselOuterWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselOuterWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pVacuumCryostatVesselPhysicalVolume;

	G4VPhysicalVolume *m_pCryostatLidFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInnerWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInnerWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOuterWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOuterWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeInsideVetoGXePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeInsideVetoLXePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeElbowPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopToShieldPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopToShieldOutsidePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeJacketCoverPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeOuterJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeOuterJacketCoverPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeElbowPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopToShieldPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopToShieldOutsidePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeJacketCoverPhysicalVolume;

	G4VPhysicalVolume *m_pCryostatLidPhysicalVolume;
	G4VPhysicalVolume *m_pVacuumCryostatLidPhysicalVolume;

	G4VPhysicalVolume *CryoSideTubeOut_phys;
	G4VPhysicalVolume *CryoSideTubeIn_phys;
	G4VPhysicalVolume *PorcupineLeft_phys;
	G4VPhysicalVolume *PorcupineRight_phys;
	G4VPhysicalVolume *PTRpipe_phys;
	G4VPhysicalVolume *PTRangle_phys;
	G4VPhysicalVolume *PTRflange_phys;
	G4VPhysicalVolume *PTR_phys;
		
	G4VPhysicalVolume *m_pResistorChainPhysicalVolume;
	G4VPhysicalVolume *m_pCablePhysicalVolume;

	static map<G4String, G4double> m_hGeometryParameters;
};

#endif // __XENON10PDETECTORCONSTRUCTION_H__

