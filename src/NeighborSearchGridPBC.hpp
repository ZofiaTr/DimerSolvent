#pragma once


#include "SBCContainerBuffer.hpp"
#include "SBCContainerIndex.hpp"
#include "particle.h"
#include <vector>


/// The class NeighborSearchGridPBC is the base class to describe a neighbor search algorithm that can be applied to particle systems.


class NeighborSearchGridPBC {

public: 
	
	/// \name Destructor
	//@{
	
	NeighborSearchGridPBC(const std::vector<Particle*>& nodeVector, const IAVector& PBC, const double& cutoffDistance);
	virtual ~NeighborSearchGridPBC();										

	//@}

	/// \name Neighbor lists
	//@{

	void														signalPositionUpdate(Particle* particle);

	virtual void												initializeNeighborLists();
	virtual void												updateNeighborLists();
	virtual void												flushPositionBuffer();

	void														onPositionChanged(Particle* particle);

	SBCContainerIndex<Particle*> const*				getNeighborIndex(unsigned int i);
	SBCContainerIndex<Particle*> const*             getNeighborIndex(Particle* structuralParticle);
	
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

		bool													operator==(const GridKey& gridKey) { return ((x==gridKey.x)&&(y==gridKey.y)&&(z==gridKey.z)); }
		bool													operator<(const GridKey& gridKey) { if (x>=gridKey.x) return false;if (y>=gridKey.y) return false;if (z>=gridKey.z) return false;return true;  }
		
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

	void														updateLists();
	

	int															roundForPBC(const double& d);
	void														correctDistanceWithPBC(Vector &distance);
	void														checkCell(GridCell* gridCellJ, int i, std::set<unsigned int>& neighbors);
	bool														particleInTheBox(Vector const& position);
	void														checkCellsAroundSmallBox(GridKey& gridKeyI, int i);
	void														checkCellsAroundBigBox(GridKey& gridKeyI, int i);
	void														getCellCoordinates(Vector const& position, Vector& coordinates);
	bool														cellCrossesTheBox(GridKey& gridKey);
	bool														cellInTheBox(GridKey& gridKey);
	bool														cellOutsideTheBox(GridKey& gridKey);
	void														getMirrorCellKey(GridKey& gridKey, GridKey& gridKeyMirror);
	void														getMirrorCellKey(GridKey& gridKey, std::vector<GridKey>& gridKeyMirror);
	void														checkCellsAroundCellByCell();

	void														checkCellsCross(GridCell* gridCellI, GridCell* gridCellJ);
	
	double*														cutoffDistance;

	SBIndex<Particle*>*											particleIndex;
	SBBuffer<Vector>*											positionBuffer;															///< The position buffer that stores the particles positions

	double*														cutoffDistanceSquared;
	double*														inverseCutoffDistance;

	SBHashMap<GridKey, GridCell*, GridKeyFunctor, GridKeyComparator>*	grid;
	SBVector<GridKey>*													gridKeyVector;

	//@}

	/// \name Neighbor lists attributes
	//@{

	SBBuffer<SBCContainerIndex<Particle*>*>*					neighborIndexBuffer;													///< The buffer containing the neighbor lists associated to each particle

	//@}

	IAVector*													boxPBC;
	double*														boxSizeX;
	double*														boxSizeY;
	double*														boxSizeZ;

};


