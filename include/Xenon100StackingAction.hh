#ifndef __XENON10PSTACKINGACTION_H__
#define __XENON10PSTACKINGACTION_H__

#include <globals.hh>
#include <G4UserStackingAction.hh>

class Xenon100AnalysisManager;

class Xenon100StackingAction: public G4UserStackingAction
{
public:
	Xenon100StackingAction(Xenon100AnalysisManager *pAnalysisManager=0);
	~Xenon100StackingAction();
  
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	virtual void NewStage();
	virtual void PrepareNewEvent();

private:
	Xenon100AnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PSTACKINGACTION_H__

