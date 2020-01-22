#include <G4Event.hh>

#include "Xenon100EventAction.hh"

Xenon100EventAction::Xenon100EventAction(Xenon100AnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

Xenon100EventAction::~Xenon100EventAction()
{
}

void
Xenon100EventAction::BeginOfEventAction(const G4Event *pEvent)
{
	if(pEvent->GetEventID() % 1000 == 0)
	{
		G4cout << G4endl;
		G4cout << "------ Begin event " << pEvent->GetEventID()
			<< "------" << G4endl;
	}
	
	if(m_pAnalysisManager)
		m_pAnalysisManager->BeginOfEvent(pEvent);
}

void Xenon100EventAction::EndOfEventAction(const G4Event *pEvent)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->EndOfEvent(pEvent);
}


