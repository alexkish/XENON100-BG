#ifndef __XENON10PEVENTACTION_H__
#define __XENON10PEVENTACTION_H__

#include <G4UserEventAction.hh>

#include "Xenon100AnalysisManager.hh"

class G4Event;

class Xenon100EventAction : public G4UserEventAction
{
public:
	Xenon100EventAction(Xenon100AnalysisManager *pAnalysisManager = 0);
	~Xenon100EventAction();

public:
	void BeginOfEventAction(const G4Event *pEvent);
	void EndOfEventAction(const G4Event *pEvent);

private:
	Xenon100AnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PEVENTACTION_H__

