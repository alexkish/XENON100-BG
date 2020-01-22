#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>

#include <numeric>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>

#include "Xenon100DetectorConstruction.hh"
#include "Xenon100LXeHit.hh"
#include "Xenon100PmtHit.hh"
#include "Xenon100PrimaryGeneratorAction.hh"
#include "Xenon100EventData.hh"

#include "Xenon100AnalysisManager.hh"

Xenon100AnalysisManager::Xenon100AnalysisManager(Xenon100PrimaryGeneratorAction *pPrimaryGeneratorAction)
{
	m_iLXeHitsCollectionID = -1;
	m_iPmtHitsCollectionID = -1;

	m_hDataFilename = "events.root";

	m_pPrimaryGeneratorAction = pPrimaryGeneratorAction;

	m_pEventData = new Xenon100EventData();
}

Xenon100AnalysisManager::~Xenon100AnalysisManager()
{
}

void
Xenon100AnalysisManager::BeginOfRun(const G4Run *pRun)
{
	m_pTreeFile = new TFile(m_hDataFilename.c_str(), "RECREATE", "File containing event data for Xenon100");
	m_pTree = new TTree("t1", "Tree containing event data for Xenon100");

	gROOT->ProcessLine("#include <vector>");

	m_pTree->Branch("eventid", &m_pEventData->m_iEventId, "eventid/I");
	m_pTree->Branch("ntpmthits", &m_pEventData->m_iNbTopPmtHits, "ntpmthits/I");
	m_pTree->Branch("nbpmthits", &m_pEventData->m_iNbBottomPmtHits, "nbpmthits/I");
	m_pTree->Branch("pmthits", "vector<int>", &m_pEventData->m_pPmtHits);
	m_pTree->Branch("etot", &m_pEventData->m_fTotalEnergyDeposited, "etot/F");
	m_pTree->Branch("nsteps", &m_pEventData->m_iNbSteps, "nsteps/I");
	
	m_pTree->Branch("trackid", "vector<int>", &m_pEventData->m_pTrackId);
	m_pTree->Branch("type", "vector<string>", &m_pEventData->m_pParticleType);
	m_pTree->Branch("parentid", "vector<int>", &m_pEventData->m_pParentId);
	m_pTree->Branch("parenttype", "vector<string>", &m_pEventData->m_pParentType);
	m_pTree->Branch("creaproc", "vector<string>", &m_pEventData->m_pCreatorProcess);
	m_pTree->Branch("edproc", "vector<string>", &m_pEventData->m_pDepositingProcess);
	m_pTree->Branch("xp", "vector<float>", &m_pEventData->m_pX);
	m_pTree->Branch("yp", "vector<float>", &m_pEventData->m_pY);
	m_pTree->Branch("zp", "vector<float>", &m_pEventData->m_pZ);
	m_pTree->Branch("ed", "vector<float>", &m_pEventData->m_pEnergyDeposited);
	m_pTree->Branch("time", "vector<float>", &m_pEventData->m_pTime);

	m_pTree->Branch("type_pri", "vector<string>", &m_pEventData->m_pPrimaryParticleType);
	m_pTree->Branch("xp_pri", &m_pEventData->m_fPrimaryX, 	"xp_pri/F");
	m_pTree->Branch("yp_pri", &m_pEventData->m_fPrimaryY, 	"yp_pri/F");
	m_pTree->Branch("zp_pri", &m_pEventData->m_fPrimaryZ, 	"zp_pri/F");
	m_pTree->Branch("e_pri",  &m_pEventData->m_fPrimaryE,	"e_pri/F");

	// switch to a new file with a certain number of events
	//m_pTree->SetMaxTreeSize(10e9);
	//m_pTree->AutoSave();
	// write everything to one file, do not switch
	m_pTree->SetMaxTreeSize(10737418240LL); // 10G bytes, don't split file automatically
	//m_pTree->SetAutoSave(1073741824L); 		// 1G bytes, esentially switches off autosave

	m_pNbEventsToSimulateParameter = new TParameter<int>("nbevents", m_iNbEventsToSimulate);
	m_pNbEventsToSimulateParameter->Write();

	// outfile for nuclei
    //G4String fileISO = "xe100_AmBe.iso";
    //std::ofstream foutIso(fileISO);
    //foutIso.close();

}

void
Xenon100AnalysisManager::EndOfRun(const G4Run *pRun)
{
	m_pTreeFile->Write();
	m_pTreeFile->Close();
}

void
Xenon100AnalysisManager::BeginOfEvent(const G4Event *pEvent)
{
	if(m_iLXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		m_iLXeHitsCollectionID = pSDManager->GetCollectionID("LXeHitsCollection");
	} 

	if(m_iPmtHitsCollectionID == -1)
	{
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		m_iPmtHitsCollectionID = pSDManager->GetCollectionID("PmtHitsCollection");
	}
}

void
Xenon100AnalysisManager::EndOfEvent(const G4Event *pEvent)
{
	G4HCofThisEvent* pHCofThisEvent = pEvent->GetHCofThisEvent();
	Xenon100LXeHitsCollection* pLXeHitsCollection = 0;
	Xenon100PmtHitsCollection* pPmtHitsCollection = 0;

	G4int iNbLXeHits = 0, iNbPmtHits = 0;
	
	if(pHCofThisEvent)
	{
		if(m_iLXeHitsCollectionID != -1)
		{
			pLXeHitsCollection = (Xenon100LXeHitsCollection *)(pHCofThisEvent->GetHC(m_iLXeHitsCollectionID));
			iNbLXeHits = (pLXeHitsCollection)?(pLXeHitsCollection->entries()):(0);
		}

		if(m_iPmtHitsCollectionID != -1)
		{
			pPmtHitsCollection = (Xenon100PmtHitsCollection *)(pHCofThisEvent->GetHC(m_iPmtHitsCollectionID));
			iNbPmtHits = (pPmtHitsCollection)?(pPmtHitsCollection->entries()):(0);
		}
	}

	if(iNbLXeHits || iNbPmtHits)
	{
		m_pEventData->m_iEventId = pEvent->GetEventID();

		m_pEventData->m_pPrimaryParticleType->push_back(m_pPrimaryGeneratorAction->GetParticleTypeOfPrimary());

		m_pEventData->m_fPrimaryE = m_pPrimaryGeneratorAction->GetEnergyOfPrimary();

		m_pEventData->m_fPrimaryX = m_pPrimaryGeneratorAction->GetPositionOfPrimary().x();
		m_pEventData->m_fPrimaryY = m_pPrimaryGeneratorAction->GetPositionOfPrimary().y();
		m_pEventData->m_fPrimaryZ = m_pPrimaryGeneratorAction->GetPositionOfPrimary().z();

		G4int iNbSteps = 0;
		G4float fTotalEnergyDeposited = 0.;

		// LXe hits
		for(G4int i=0; i<iNbLXeHits; i++)
		{
			Xenon100LXeHit *pHit = (*pLXeHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			//if(pHit->GetParticleType()!="opticalphoton" && pHit->GetParticleType()!="gamma" && pHit->GetParticleType()!="e-" && pHit->GetParticleType()!="e+")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				fTotalEnergyDeposited += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}
		};

		m_pEventData->m_iNbSteps = iNbSteps;
		m_pEventData->m_fTotalEnergyDeposited = fTotalEnergyDeposited;

		G4int iNbTopPmts = (G4int) Xenon100DetectorConstruction::GetGeometryParameter("NbTopPmts");
		G4int iNbBottomPmts = (G4int) Xenon100DetectorConstruction::GetGeometryParameter("NbBottomPmts");
		G4int iNbTopVetoPmts = (G4int) Xenon100DetectorConstruction::GetGeometryParameter("NbTopVetoPmts");
		G4int iNbBottomVetoPmts = (G4int) Xenon100DetectorConstruction::GetGeometryParameter("NbBottomVetoPmts");

		m_pEventData->m_pPmtHits->resize(iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts, 0);

		// Pmt hits
		for(G4int i=0; i<iNbPmtHits; i++)
			(*(m_pEventData->m_pPmtHits))[(*pPmtHitsCollection)[i]->GetPmtNb()]++;

		m_pEventData->m_iNbTopPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin(), m_pEventData->m_pPmtHits->begin()+iNbTopPmts, 0);
		m_pEventData->m_iNbBottomPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts, m_pEventData->m_pPmtHits->end(), 0);

//        if((fTotalEnergyDeposited > 0. || iNbPmtHits > 0) && !FilterEvent(m_pEventData))
//		if(fTotalEnergyDeposited > 0. || iNbPmtHits > 0)
		if(fTotalEnergyDeposited > 0. || iNbPmtHits > 0)
			m_pTree->Fill();

		m_pEventData->Clear();
	}
}

void
Xenon100AnalysisManager::Step(const G4Step *pStep)
{

}

G4bool
Xenon100AnalysisManager::FilterEvent(Xenon100EventData *pEventData)
{
	G4double dEnergyDepositedSensitiveRegion = 0.;

	vector<float> *pX = pEventData->m_pX;
	vector<float> *pY = pEventData->m_pY;
	vector<float> *pZ = pEventData->m_pZ;
	vector<float> *pEnergyDeposited = pEventData->m_pEnergyDeposited;

	const G4double dDriftLength = Xenon100DetectorConstruction::GetGeometryParameter("DriftLength");
	const G4double dRadius = Xenon100DetectorConstruction::GetGeometryParameter("TeflonCylinderInnerRadius");

	for(G4int i=0; i<pEnergyDeposited->size(); i++)
	{
		if((*pZ)[i] < 0.-21.5 && (*pZ)[i] > -dDriftLength-21.5 && std::sqrt((*pX)[i]*(*pX)[i] + (*pY)[i]*(*pY)[i]) < dRadius)
			dEnergyDepositedSensitiveRegion += (*pEnergyDeposited)[i];
	}

//    if(dEnergyDepositedSensitiveRegion > 0. && dEnergyDepositedSensitiveRegion < 100.)
	if(dEnergyDepositedSensitiveRegion > 0.)
		return false;
	else
		return true;
}

	
