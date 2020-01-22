#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "Xenon100PmtHit.hh"

G4Allocator<Xenon100PmtHit> Xenon100PmtHitAllocator;

Xenon100PmtHit::Xenon100PmtHit() {}

Xenon100PmtHit::~Xenon100PmtHit() {}

Xenon100PmtHit::Xenon100PmtHit(const Xenon100PmtHit &hXenon100PmtHit):G4VHit()
{
	m_hPosition = hXenon100PmtHit.m_hPosition;
	m_dTime = hXenon100PmtHit.m_dTime;
	m_iPmtNb = hXenon100PmtHit.m_iPmtNb;
}

const Xenon100PmtHit &
Xenon100PmtHit::operator=(const Xenon100PmtHit &hXenon100PmtHit)
{
	m_hPosition = hXenon100PmtHit.m_hPosition;
	m_dTime = hXenon100PmtHit.m_dTime;
	m_iPmtNb = hXenon100PmtHit.m_iPmtNb;
	
	return *this;
}

G4int
Xenon100PmtHit::operator==(const Xenon100PmtHit &hXenon100PmtHit) const
{
	return ((this == &hXenon100PmtHit) ? (1) : (0));
}

void Xenon100PmtHit::Draw()
{
//    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//    
//    if(pVVisManager)
//    {
//        G4Circle hCircle(m_hPosition);
//        G4Colour hColour(1.000, 0.973, 0.184);
//        G4VisAttributes hVisAttributes(hColour);
//        
//        hCircle.SetScreenSize(0.1);
//        hCircle.SetFillStyle(G4Circle::filled);
//        hCircle.SetVisAttributes(hVisAttributes);
//        pVVisManager->Draw(hCircle);
//    }
}

void Xenon100PmtHit::Print()
{
	G4cout << "Pmt hit ---> " 
		<< "Pmt#" << m_iPmtNb
		<< " Position: " << m_hPosition.x()/mm
		<< " " << m_hPosition.y()/mm
		<< " " << m_hPosition.z()/mm
		<< " mm"
		<< " Time: " << m_dTime/s << " s" << G4endl;
}

