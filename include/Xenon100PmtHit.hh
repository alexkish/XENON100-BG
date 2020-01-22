#ifndef __XENON10PPMTHIT_H__
#define __XENON10PPMTHIT_H__

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>

class Xenon100PmtHit: public G4VHit
{
public:
	Xenon100PmtHit();
	~Xenon100PmtHit();
	Xenon100PmtHit(const Xenon100PmtHit &);
	const Xenon100PmtHit & operator=(const Xenon100PmtHit &);
	G4int operator==(const Xenon100PmtHit &) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	void Draw();
	void Print();

public:
	void SetPosition(G4ThreeVector hPosition) { m_hPosition = hPosition; }
	void SetTime(G4double dTime) { m_dTime = dTime; }
	void SetPmtNb(G4int iPmtNb) { m_iPmtNb = iPmtNb; }

	G4ThreeVector GetPosition() { return m_hPosition; }
	G4double GetTime() { return m_dTime; }
	G4int GetPmtNb() { return m_iPmtNb; }

private:
	G4ThreeVector m_hPosition;
	G4double m_dTime;
	G4int m_iPmtNb;
};

typedef G4THitsCollection<Xenon100PmtHit> Xenon100PmtHitsCollection;

extern G4Allocator<Xenon100PmtHit> Xenon100PmtHitAllocator;

inline void*
Xenon100PmtHit::operator new(size_t)
{
	return((void *) Xenon100PmtHitAllocator.MallocSingle());
}

inline void
Xenon100PmtHit::operator delete(void *pXenon100PmtHit)
{
	Xenon100PmtHitAllocator.FreeSingle((Xenon100PmtHit*) pXenon100PmtHit);
}

#endif // __XENON10PPMTHIT_H__

