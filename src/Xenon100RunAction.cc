#include <Randomize.hh>

#include <sys/time.h>

#include "Xenon100AnalysisManager.hh"

#include "Xenon100RunAction.hh"

Xenon100RunAction::Xenon100RunAction(Xenon100AnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

Xenon100RunAction::~Xenon100RunAction()
{

}

void
Xenon100RunAction::BeginOfRunAction(const G4Run *pRun)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->BeginOfRun(pRun);

	struct timeval hTimeValue;
	gettimeofday(&hTimeValue, NULL);

	CLHEP::HepRandom::setTheEngine(new CLHEP::DRand48Engine);
	CLHEP::HepRandom::setTheSeed(hTimeValue.tv_usec);
}

void
Xenon100RunAction::EndOfRunAction(const G4Run *pRun)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->EndOfRun(pRun);
}

