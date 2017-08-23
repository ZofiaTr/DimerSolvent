#include "NeighborSearchGrid.hpp"
#include "math.h"

NeighborSearchGrid::NeighborSearchGrid(const std::vector<Particle*>& nodeVector, const double& c) {

	// collect particles from the node vector

	particleIndex=new SBIndex<Particle*>;
	unsigned int nodeVectorSize=(unsigned int)nodeVector.size();

	for (unsigned int i=0;i<nodeVectorSize;i++) {

		Particle* currentNode=nodeVector[i];
		particleIndex->push_back(currentNode);

	}

	unsigned int nParticles=particleIndex->size(); // the number of particles

	// setup the position buffer

	positionBuffer=new SBBuffer<Vector>(nParticles);

	for (unsigned int i=0;i<nParticles;i++) {

		Particle* currentParticle=particleIndex->getObject(i);
		Vector currentPosition=currentParticle->getPosition();
		positionBuffer->setElement(i,currentPosition);

	}

	positionBuffer->flush();
	
	// setup the neighbor index buffer

	neighborIndexBuffer=new SBBuffer<SBCContainerIndex<Particle*>*>(nParticles);
	
	for (unsigned int i=0;i<nParticles;i++) {

		SBCContainerIndex<Particle*>* currentVector=new SBCContainerIndex<Particle*>;
		neighborIndexBuffer->setElement(i,currentVector);
		
	}

	// setup the grid

	inverseCutoffDistance=new (double)(1.0/c);
	cutoffDistanceSquared=new (double)(c*c);

	grid=new SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>;
	gridKeyVector=new SBVector<GridKey>(nParticles);

	for (unsigned int i=0;i<nParticles;i++) {

		Particle* currentParticle=particleIndex->getObject(i);
		Vector currentPosition=currentParticle->getPosition();

		GridKey gridKey((int)floor(currentPosition.v[0] * (*inverseCutoffDistance)), (int)floor(currentPosition.v[1] * (*inverseCutoffDistance)), (int)floor(currentPosition.v[2]*(*inverseCutoffDistance)));
		gridKeyVector->push_back(gridKey);

		// query the grid

		GridCell* gridCell=0;
		grid->getValue(gridKey,gridCell);
		
		if (gridCell!=0) gridCell->addParticle(currentParticle);
		else {

			gridCell=new GridCell(currentParticle);
			grid->insert(gridKey,gridCell);

		}

	}
	
}

NeighborSearchGrid::~NeighborSearchGrid() {

	// delete the grid

	for (SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>::iterator i=grid->begin(); i!=grid->end(); ++i) delete i->getValue();
	delete grid;
	delete gridKeyVector;

	// delete the neighbor index buffer

	delete neighborIndexBuffer;

	// delete the position buffer

	delete positionBuffer;

	// delete the pointer index

	delete particleIndex;

	// delete cutoffs

	delete inverseCutoffDistance;
	delete cutoffDistanceSquared;

}

void NeighborSearchGrid::initializeNeighborLists() {

	unsigned int nParticles=particleIndex->size();

	for (unsigned int i=0;i<nParticles;i++) {

		Vector const& positionI=(*positionBuffer)[i];
		Particle* particleI=particleIndex->getObject(i);

		GridKey& gridKeyI=(*gridKeyVector)[i];

		for (int x=gridKeyI.x-1;x<=gridKeyI.x+1;x++) {

			for (int y=gridKeyI.y-1;y<=gridKeyI.y+1;y++) {

				for (int z=gridKeyI.z-1;z<=gridKeyI.z+1;z++) {

					GridKey gridKey(x,y,z);

					// query the grid

					GridCell* gridCell=0;
					grid->getValue(gridKey,gridCell);

					if (gridCell==0) continue; // empty neighboring cell

					unsigned int nParticlesJ=gridCell->particleIndex->size();

					for (unsigned int j=0;j<nParticlesJ;j++) {

						Particle* particleJ=gridCell->particleIndex->getObject(j);
						unsigned int particleJIndex=particleIndex->getIndex(particleJ);

						if (particleJIndex<=i) continue;

						Vector const& positionJ=(*positionBuffer)[particleJIndex];

						if ((positionJ-positionI).norm2()<=*cutoffDistanceSquared) {

							(*neighborIndexBuffer)[i]->push_back(particleJ);
							(*neighborIndexBuffer)[particleJIndex]->push_back(particleI);

						}

					} // particle j loop

				} // z loop

			} // y loop

		} // x loop

	} // particle i loop

}

void NeighborSearchGrid::updateNeighborLists() {

	unsigned int nParticles=positionBuffer->getNumberOfChangedElements();

	// prune the list of neighbors of particles (because of particles which have moved apart)

	for (unsigned int i=0;i<nParticles;i++) {

		unsigned int particleIIndex=positionBuffer->getChangedElementIndex(i);
		Particle* particleI=particleIndex->getObject(particleIIndex);
		Vector const& positionI=(*positionBuffer)[particleIIndex];

		SBCContainerIndex<Particle*>* neighborIndexI=(*neighborIndexBuffer)[particleIIndex];
	
		unsigned int j=0;

		while (j<neighborIndexI->size()) {

			Particle* particleJ=neighborIndexI->getObject(j);
			unsigned int particleJIndex=particleIndex->getIndex(particleJ);
			Vector const& positionJ=(*positionBuffer)[particleJIndex];

			if ((positionJ-positionI).norm2()<=*cutoffDistanceSquared) j++; // keep this neighbor
			else {

				neighborIndexI->eraseIndex(j);
				(*neighborIndexBuffer)[particleJIndex]->eraseObject(particleI);

			}

		}

	}

	// update the grid

	for (unsigned int i=0;i<nParticles;i++) {

		unsigned int particleIIndex=positionBuffer->getChangedElementIndex(i);
		Particle* particleI=particleIndex->getObject(particleIIndex);
		Vector const& positionI=(*positionBuffer)[particleIIndex];

		GridKey& gridKeyIOld=(*gridKeyVector)[particleIIndex];
		GridKey  gridKeyINew((int)floor(positionI.v[0]*(*inverseCutoffDistance)),(int)floor(positionI.v[1]*(*inverseCutoffDistance)),(int)floor(positionI.v[2]*(*inverseCutoffDistance))); // new grid key

		if (gridKeyIOld==gridKeyINew) continue; // the particle stays in the same grid cell

		// remove the particle from the old grid cell

		GridCell* gridCellIOld=0;
		grid->getValue(gridKeyIOld,gridCellIOld);

		gridCellIOld->removeParticle(particleI);

		if (gridCellIOld->particleIndex->size()==0) { // the old cell is now empty, remove it from the grid

			delete gridCellIOld;
			grid->erase(gridKeyIOld);

		}

		// add the particle to the new grid cell

		GridCell* gridCellINew=0;
		grid->getValue(gridKeyINew,gridCellINew);

		if (gridCellINew!=0) gridCellINew->addParticle(particleI);
		else {

			gridCellINew=new GridCell(particleI);
			grid->insert(gridKeyINew,gridCellINew);

		}

		// store the new grid key

		(*gridKeyVector)[particleIIndex]=gridKeyINew;

	}

	// detect new neighbors

	for (unsigned int i=0;i<nParticles;i++) {

		unsigned int particleIIndex=positionBuffer->getChangedElementIndex(i);
		Particle* particleI=particleIndex->getObject(particleIIndex);
		Vector const& positionI=(*positionBuffer)[particleIIndex];

		GridKey  gridKeyI((int)floor(positionI.v[0] * (*inverseCutoffDistance)), (int)floor(positionI.v[1] * (*inverseCutoffDistance)), (int)floor(positionI.v[2] * (*inverseCutoffDistance))); // new grid key

		for (int x=gridKeyI.x-1;x<=gridKeyI.x+1;x++) {

			for (int y=gridKeyI.y-1;y<=gridKeyI.y+1;y++) {

				for (int z=gridKeyI.z-1;z<=gridKeyI.z+1;z++) {

					GridKey gridKey(x,y,z);

					// query the grid

					GridCell* gridCell=0;
					grid->getValue(gridKey,gridCell);

					if (gridCell==0) continue; // empty neighboring cell

					unsigned int nParticlesJ=gridCell->particleIndex->size();

					for (unsigned int j=0;j<nParticlesJ;j++) {

						Particle* particleJ=gridCell->particleIndex->getObject(j);

						if (particleJ==particleI) continue;

						unsigned int particleJIndex=particleIndex->getIndex(particleJ);

						Vector const& positionJ=(*positionBuffer)[particleJIndex];

						if ((positionJ-positionI).norm2()<=*cutoffDistanceSquared) {

							(*neighborIndexBuffer)[particleIIndex]->push_back(particleJ);
							(*neighborIndexBuffer)[particleJIndex]->push_back(particleI);

						}

					} // particle j loop

				} // z loop

			} // y loop

		} // x loop

	} // particle i loop

	// flush the position buffer, since we have used all the updates

	positionBuffer->flush();

}

void NeighborSearchGrid::signalPositionUpdate(Particle* particle) {

	unsigned int index = particleIndex->getIndex(particle);
	positionBuffer->setElement(index, particle->getPosition());

}

void NeighborSearchGrid::print(unsigned int offset) const {

	unsigned int nParticles=particleIndex->size();

	std::cout << "Neighboring particles: " << std::endl;

	unsigned int nPairs=0;
	for (unsigned int i=0;i<nParticles;i++) {

		SBCContainerIndex<Particle*>* neighborIndexI=(*neighborIndexBuffer)[i];
		
		for (unsigned int j=0;j<neighborIndexI->size();j++) {
			
			std::cout << "Neighbor pair: " << i << " " << particleIndex->getIndex(neighborIndexI->getObject(j)) << std::endl;
			nPairs++;

		}

	}

	std::cout << "Number of grid cells: " << grid->size() << std::endl;
	std::cout << "Number of pairs: " << nPairs << std::endl << std::endl;

}

SBCContainerIndex<Particle*> const* NeighborSearchGrid::getNeighborIndex(unsigned int i) { return (*neighborIndexBuffer)[i]; }

SBCContainerIndex<Particle*> const* NeighborSearchGrid::getNeighborIndex(Particle* particle) { 

	if (!particleIndex->hasIndex(particle)) return 0;
	return (*neighborIndexBuffer)[particleIndex->getIndex(particle)]; 

}
