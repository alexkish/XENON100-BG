#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "Xenon100LXeHit.hh"

G4Allocator<Xenon100LXeHit> Xenon100LXeHitAllocator;

Xenon100LXeHit::Xenon100LXeHit() {}

Xenon100LXeHit::~Xenon100LXeHit()
{
	if(m_pParticleType) delete m_pParticleType;
	if(m_pParentType) delete m_pParentType;
	if(m_pCreatorProcess) delete m_pCreatorProcess;
	if(m_pDepositingProcess) delete m_pDepositingProcess;
}

Xenon100LXeHit::Xenon100LXeHit(const Xenon100LXeHit &hXenon100LXeHit):G4VHit()
{
	m_iTrackId = hXenon100LXeHit.m_iTrackId;
	m_iParentId = hXenon100LXeHit.m_iParentId;
	m_pParticleType = hXenon100LXeHit.m_pParticleType;
	m_pParentType = hXenon100LXeHit.m_pParentType ;
	m_pCreatorProcess = hXenon100LXeHit.m_pCreatorProcess ;
	m_pDepositingProcess = hXenon100LXeHit.m_pDepositingProcess ;
	m_hPosition = hXenon100LXeHit.m_hPosition;
	m_dEnergyDeposited = hXenon100LXeHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXenon100LXeHit.m_dKineticEnergy ;
	m_dTime = hXenon100LXeHit.m_dTime;
}

const Xenon100LXeHit &
Xenon100LXeHit::operator=(const Xenon100LXeHit &hXenon100LXeHit)
{
	m_iTrackId = hXenon100LXeHit.m_iTrackId;
	m_iParentId = hXenon100LXeHit.m_iParentId;
	m_pParticleType = hXenon100LXeHit.m_pParticleType;
	m_pParentType = hXenon100LXeHit.m_pParentType ;
	m_pCreatorProcess = hXenon100LXeHit.m_pCreatorProcess ;
	m_pDepositingProcess = hXenon100LXeHit.m_pDepositingProcess ;
	m_hPosition = hXenon100LXeHit.m_hPosition;
	m_dEnergyDeposited = hXenon100LXeHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXenon100LXeHit.m_dKineticEnergy ;
	m_dTime = hXenon100LXeHit.m_dTime;
	
	return *this;
}

G4int
Xenon100LXeHit::operator==(const Xenon100LXeHit &hXenon100LXeHit) const
{
	return ((this == &hXenon100LXeHit) ? (1) : (0));
}

void Xenon100LXeHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	
	if(pVVisManager)
	{
		G4Circle hCircle(m_hPosition);
		G4Colour hColour(1.000, 0.973, 0.184);
		G4VisAttributes hVisAttributes(hColour);
		
		hCircle.SetScreenSize(0.1);
		hCircle.SetFillStyle(G4Circle::filled);
		hCircle.SetVisAttributes(hVisAttributes);
		pVVisManager->Draw(hCircle);
	}
}

void Xenon100LXeHit::Print()
{
	G4cout << "-------------------- LXe hit --------------------" 
		<< "Id: " << m_iTrackId
		<< " Particle: " << *m_pParticleType
		<< " ParentId: " << m_iParentId
		<< " ParentType: " << *m_pParentType << G4endl
		<< "CreatorProcess: " << *m_pCreatorProcess
		<< " DepositingProcess: " << *m_pDepositingProcess << G4endl
		<< "Position: " << m_hPosition.x()/mm
		<< " " << m_hPosition.y()/mm
		<< " " << m_hPosition.z()/mm
		<< " mm" << G4endl
		<< "EnergyDeposited: " << m_dEnergyDeposited/keV << " keV"
		<< " KineticEnergyLeft: " << m_dKineticEnergy/keV << " keV"
		<< " Time: " << m_dTime/s << " s" << G4endl;
}

