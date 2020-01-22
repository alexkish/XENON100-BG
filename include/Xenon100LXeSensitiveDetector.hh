#ifndef __XENON10PLXESENSITIVEDETECTOR_H__
#define __XENON10PLXESENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "Xenon100LXeHit.hh"

class G4Step;
class G4HCofThisEvent;

class Xenon100LXeSensitiveDetector: public G4VSensitiveDetector
{
public:
	Xenon100LXeSensitiveDetector(G4String hName);
	~Xenon100LXeSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	Xenon100LXeHitsCollection* m_pLXeHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif // __XENON10PLXESENSITIVEDETECTOR_H__

