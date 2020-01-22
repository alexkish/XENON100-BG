#ifndef __XENON10PPMTSENSITIVEDETECTOR_H__
#define __XENON10PPMTSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "Xenon100PmtHit.hh"

class G4Step;
class G4HCofThisEvent;

class Xenon100PmtSensitiveDetector: public G4VSensitiveDetector
{
public:
	Xenon100PmtSensitiveDetector(G4String hName);
	~Xenon100PmtSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	Xenon100PmtHitsCollection* m_pPmtHitsCollection;
};

#endif // __XENON10PPMTSENSITIVEDETECTOR_H__

