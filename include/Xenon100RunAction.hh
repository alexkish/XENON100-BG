#ifndef __XENON10PRUNACTION_H__
#define __XENON10PRUNACTION_H__

#include <G4UserRunAction.hh>

class G4Run;

class Xenon100AnalysisManager;

class Xenon100RunAction: public G4UserRunAction
{
public:
	Xenon100RunAction(Xenon100AnalysisManager *pAnalysisManager=0);
	~Xenon100RunAction();

public:
	void BeginOfRunAction(const G4Run *pRun);
	void EndOfRunAction(const G4Run *pRun);

private:
	Xenon100AnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PRUNACTION_H__

