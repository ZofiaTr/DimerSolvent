#pragma once


#include "SBCContainerBuffer.hpp"
#include "SBCContainerIndex.hpp"
#include "particle.h"
#include <vector>

/// The class NeighborSearchGrid implements a grid-based neighbor search algorithm that can be applied to particles

class NeighborSearchGrid {

public: 
	
	/// \name Constructor and destructor
	//@{
	
	NeighborSearchGrid(const std::vector<Particle*>& nodeVector, const double& cutoffDistance);															
	virtual ~NeighborSearchGrid();										

	//@}

	/// \name Neighbor lists
	//@{

	void														signalPositionUpdate(Particle* particle);

	virtual void												initializeNeighborLists();
	virtual void												updateNeighborLists();

	SBCContainerIndex<Particle*> const*							getNeighborIndex(unsigned int i);
	SBCContainerIndex<Particle*> const*							getNeighborIndex(Particle* particle);

	//@}

	/// \name Debugging
	//@{

	virtual void												print(unsigned int offset=0) const;

	//@}


protected:
	
	class GridKey {

	public:

		GridKey() { x=0;y=0;z=0; }
		GridKey(const GridKey& gridKey) { x=gridKey.x;y=gridKey.y;z=gridKey.z; }
		GridKey(int xx, int yy, int zz) { x=xx;y=yy;z=zz; }
		virtual ~GridKey() {}

		bool														operator==(const GridKey& gridKey) { return ((x==gridKey.x)&&(y==gridKey.y)&&(z==gridKey.z)); }

		int x;
		int y;
		int z;

	};

	class GridKeyFunctor {

	public:

		GridKeyFunctor() {}
		virtual ~GridKeyFunctor() {}

		static unsigned long										hashCode(const GridKey& key) { 

			unsigned long hash = (key.x << 6) + (key.y << 16) - key.z;
			return hash;

		}

	};

	class GridKeyComparator {

	public: 

		GridKeyComparator() {}
		virtual ~GridKeyComparator() {}

		static bool													equal(const GridKey& key1, const GridKey& key2) { 

			return ((key1.x==key2.x)&&(key1.y==key2.y)&&(key1.z==key2.z));

		}

	};

	class GridCell {

	public: 

		GridCell(Particle* p) { particleIndex=new SBCContainerIndex<Particle*>;particleIndex->push_back(p); }
		virtual ~GridCell() { delete particleIndex; }

		void														addParticle(Particle* p) { particleIndex->push_back(p); }
		void														removeParticle(Particle* p) { particleIndex->eraseObject(p); }

		SBCContainerIndex<Particle*>*								particleIndex;

	};

	/// \name Grid attributes
	//@{

	double*														cutoffDistanceSquared;
	double*														inverseCutoffDistance;

	//@}

	/// \name Particle attributes
	//@{

	SBIndex<Particle*>*											particleIndex;															///< The particle system the neighbor search algorithm is attached to
	SBBuffer<Vector>*											positionBuffer;															///< The position buffer that stores the particles positions


	//@}

	/// \name Neighbor lists attributes
	//@{

	SBBuffer<SBCContainerIndex<Particle*>*>*					neighborIndexBuffer;													///< The buffer containing the neighbor lists associated to each particle

	//@}

	SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>*	grid;
	SBVector<GridKey>*												gridKeyVector;

};

