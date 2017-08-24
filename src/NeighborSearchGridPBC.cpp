#include "NeighborSearchGridPBC.hpp"
#include "math.h"

#include <set>

NeighborSearchGridPBC::NeighborSearchGridPBC(const std::vector<Particle*>& nodeVector, const IAVector& PBC, const double& c) {

	cutoffDistance = new double(c);
	inverseCutoffDistance = new double(1.0 / c);
	cutoffDistanceSquared = new double(c*c);

	boxPBC = new IAVector(PBC);
	boxSizeX = new double(boxPBC->i[0].i[1] - boxPBC->i[0].i[0]);
	boxSizeY = new double(boxPBC->i[1].i[1] - boxPBC->i[1].i[0]);
	boxSizeZ = new double(boxPBC->i[2].i[1] - boxPBC->i[2].i[0]);

	// collect particles from the node vector

	particleIndex = new SBIndex<Particle*>;
	unsigned int nodeVectorSize = (unsigned int)nodeVector.size();

	for (unsigned int i = 0; i<nodeVectorSize; i++) {

		Particle* currentNode = nodeVector[i];
		particleIndex->push_back(currentNode);

	}

	unsigned int nParticles = particleIndex->size(); // the number of particles

	neighborIndexBuffer = new SBBuffer<SBCContainerIndex<Particle*>*>(nParticles);

	for (unsigned int i = 0; i<nParticles; i++) {

		SBCContainerIndex<Particle*>* currentVector = new SBCContainerIndex<Particle*>;
		neighborIndexBuffer->setElement(i, currentVector);

	}

	neighborIndexBuffer->flush();

	// setup the position buffer

	positionBuffer = new SBBuffer<Vector>(nParticles);

	for (unsigned int i = 0; i<nParticles; i++) {

		Particle* currentParticle = particleIndex->getObject(i);
		Vector currentPosition = currentParticle->getPosition();
		positionBuffer->setElement(i, currentPosition);

	}

	// allocate memory

	grid=new SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>;
	gridKeyVector=new SBVector<GridKey>(nParticles);

	// build grid

	for (unsigned int i=0;i<nParticles;i++) {

		Vector const& currentPosition = (*positionBuffer)[i];
		Particle* currentParticle=particleIndex->getObject(i);

		Vector cellCoordinates;
		getCellCoordinates(currentPosition, cellCoordinates);
		GridKey gridKey((int)cellCoordinates[0], (int)cellCoordinates[1], (int)cellCoordinates[2]);

		gridKeyVector->push_back(gridKey);

		// query the grid

		GridCell* gridCell=0;
		grid->getValue(gridKey,gridCell);

		if (gridCell!=0)
			gridCell->addParticle(currentParticle);
		else {

			gridCell=new GridCell(currentParticle);
			grid->insert(gridKey,gridCell);

		}

	}

}

NeighborSearchGridPBC::~NeighborSearchGridPBC() {

	// delete what's inside the map

	for (SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>::iterator i=grid->begin(); i!=grid->end(); ++i) delete i->getValue();

	// delete the map

	delete grid;

	// delete the vector

	delete gridKeyVector;

	delete cutoffDistance;

}

void NeighborSearchGridPBC::getCellCoordinates(Vector const& currentPosition, Vector& coordinates) {

	int x = (int)floor(currentPosition.v[0]*(*inverseCutoffDistance));
	int y = (int)floor(currentPosition.v[1]*(*inverseCutoffDistance));
	int z = (int)floor(currentPosition.v[2]*(*inverseCutoffDistance));

	// x
	int gridKeyMostLeftX = (int)floor((boxPBC->i[0].i[0])*(*inverseCutoffDistance));
	int gridKeyMostRightX = (int)floor((boxPBC->i[0].i[1])*(*inverseCutoffDistance));

	// most left x
	if ((x == gridKeyMostLeftX) && (fabs(gridKeyMostLeftX*(*cutoffDistance) - boxPBC->i[0].i[0]) > 0.1))
		x+=1;
	// most right x
	else if ((x == gridKeyMostRightX) && (fabs((gridKeyMostRightX+1)*(*cutoffDistance) - boxPBC->i[0].i[1]) > 0.1))
		x-=1;

	// y
	int gridKeyMostLeftY = (int)floor((boxPBC->i[1].i[0])*(*inverseCutoffDistance));
	int gridKeyMostRightY = (int)floor((boxPBC->i[1].i[1])*(*inverseCutoffDistance));

	// most left y
	if ((y == gridKeyMostLeftY) && (fabs(gridKeyMostLeftY*(*cutoffDistance) - boxPBC->i[1].i[0]) > 0.1))
		y+=1;
	// most right y
	else if ((y == gridKeyMostRightY) && (fabs((gridKeyMostRightY+1)*(*cutoffDistance) - boxPBC->i[1].i[1]) > 0.1))
		y-=1;

	// z
	int gridKeyMostLeftZ = (int)floor((boxPBC->i[2].i[0])*(*inverseCutoffDistance));
	int gridKeyMostRightZ = (int)floor((boxPBC->i[2].i[1])*(*inverseCutoffDistance));

	// most left z
	if (z == gridKeyMostLeftZ && (fabs(gridKeyMostLeftZ*(*cutoffDistance) - boxPBC->i[2].i[0]) > 0.1))
		z+=1;
	// most right z
	else if (z == gridKeyMostRightZ && (fabs((gridKeyMostRightZ+1)*(*cutoffDistance) - boxPBC->i[2].i[1]) > 0.1))
		z-=1;

	coordinates.v[0]=x;
	coordinates.v[1]=y;
	coordinates.v[2]=z;

}

// determines whether the particle with a certain position lies inside the PBC box
bool NeighborSearchGridPBC::particleInTheBox(Vector const& position){

	if (position[0] >= boxPBC->i[0].i[0] && position[0] <= boxPBC->i[0].i[1] &&
		position[1] >= boxPBC->i[1].i[0] && position[1] <= boxPBC->i[1].i[1] &&
		position[2] >= boxPBC->i[2].i[0] && position[2] <= boxPBC->i[2].i[1] )
		return true;

	return false;
}


int NeighborSearchGridPBC::roundForPBC(const double& d) {

	if (d >= 0) {

		if (d - ((int)d) < 0.5)
			return (int)d;

		return (int)d + 1;

	}
	else {

		if (d - ((int)d) > -0.5)
			return (int)d;

		return (int)d - 1;

	}

}
void NeighborSearchGridPBC::correctDistanceWithPBC(Vector& distance) {

	distance[0] -= roundForPBC(distance[0] / (*boxSizeX))*(*boxSizeX);
	distance[1] -= roundForPBC(distance[1] / (*boxSizeY))*(*boxSizeY);
	distance[2] -= roundForPBC(distance[2] / (*boxSizeZ))*(*boxSizeZ);

}

void NeighborSearchGridPBC::checkCell(GridCell* gridCellJ, int i, std::set<unsigned int>& neighbors){

	Vector const& positionI=(*positionBuffer)[i];
	Particle* particleI=particleIndex->getObject(i);

	unsigned int nParticlesJ=gridCellJ->particleIndex->size();

	for (unsigned int j=0;j<nParticlesJ;j++) {

		Particle* particleJ=gridCellJ->particleIndex->getObject(j);
		unsigned int particleJIndex=particleIndex->getIndex(particleJ);

		//if (particleJIndex<=(unsigned int)i) continue;
		if (particleJIndex == (unsigned int)i) continue;

		if (neighbors.count(particleJIndex)>0)
			// the cell was already visited
			return;

		Vector const& positionJ=(*positionBuffer)[particleJIndex];
		if (!particleInTheBox(positionJ)) continue;

		Vector distance = positionJ-positionI;

		bool replica = false;
		if (distance.norm2() > (*cutoffDistanceSquared)) {

			replica = true;
			correctDistanceWithPBC(distance);

			if (distance.norm2() > (*cutoffDistanceSquared))
				continue;

		}

		if (!(*neighborIndexBuffer)[i]->hasIndex(particleJ))
			(*neighborIndexBuffer)[i]->push_back(particleJ);

		if (!(*neighborIndexBuffer)[particleJIndex]->hasIndex(particleI))
			(*neighborIndexBuffer)[particleJIndex]->push_back(particleI);

		neighbors.insert(particleJIndex);

	}

}

void NeighborSearchGridPBC::checkCellsCross(GridCell* gridCellI, GridCell* gridCellJ) {

	unsigned int nParticlesI=gridCellI->particleIndex->size();
	unsigned int nParticlesJ=gridCellJ->particleIndex->size();

	for (unsigned int i=0; i<nParticlesI; i++) {

		Particle* particleI=gridCellI->particleIndex->getObject(i);
		unsigned int particleIIndex=particleIndex->getIndex(particleI);
		Vector const& positionI=(*positionBuffer)[particleIIndex];
		if (!particleInTheBox(positionI)) {
			std::cout << "particle out !\n";
			positionI.print();
			continue;
		}

		for (unsigned int j=0; j<nParticlesJ; j++) {

			Particle* particleJ=gridCellJ->particleIndex->getObject(j);
			unsigned int particleJIndex=particleIndex->getIndex(particleJ);

			if (gridCellI == gridCellJ && particleJIndex<=particleIIndex)
				continue;

			Vector const& positionJ=(*positionBuffer)[particleJIndex];
			if (!particleInTheBox(positionJ)) {
				std::cout << "particle out !\n";
				positionJ.print();
				continue;
			}

			Vector distance = positionJ-positionI;

			bool replica = false;
			if (distance.norm2() > (*cutoffDistanceSquared)) {

				replica = true;
				correctDistanceWithPBC(distance);

				if (distance.norm2() > (*cutoffDistanceSquared))
					continue;

			}

			(*neighborIndexBuffer)[particleIIndex]->push_back(particleJ);
			(*neighborIndexBuffer)[particleJIndex]->push_back(particleI);

		}

	}

}

void NeighborSearchGridPBC::initializeNeighborLists() {

	updateLists();

	positionBuffer->flush();

}

void NeighborSearchGridPBC::updateLists() {

	unsigned int nChangedParticles = positionBuffer->getNumberOfChangedElements();

	if (*boxSizeX >= *cutoffDistance * 3 && *boxSizeY >= *cutoffDistance * 3 && *boxSizeZ >= *cutoffDistance * 3) {

		for (unsigned int i = 0; i<nChangedParticles; i++) {

			// check whether particle belongs to the box

			unsigned int particleIIndex = positionBuffer->getChangedElementIndex(i);

			if (!particleInTheBox((*positionBuffer)[particleIIndex])) continue;

			checkCellsAroundBigBox((*gridKeyVector)[particleIIndex], particleIIndex);

		}

	}
	else  {

		for (unsigned int i = 0; i<nChangedParticles; i++) {

			// check whether particle belongs to the box

			unsigned int particleIIndex = positionBuffer->getChangedElementIndex(i);

			if (!particleInTheBox((*positionBuffer)[particleIIndex])) continue;

			checkCellsAroundSmallBox((*gridKeyVector)[particleIIndex], particleIIndex);

		}

	}

	//unsigned int nChangedParticles = positionBuffer->getNumberOfChangedElements();

	//if (*boxSizeX >= *cutoffDistance*3 && *boxSizeY >= *cutoffDistance*3 && *boxSizeZ >= *cutoffDistance*3) {

	//	checkCellsAroundCellByCell();

	//}
	//else  {

	//	for (unsigned int i = 0; i<nChangedParticles; i++) {

	//		// check whether particle belongs to the box

	//		if (!particleInTheBox((*positionBuffer)[i])) continue;

	//		checkCellsAroundSmallBox((*gridKeyVector)[i], i);

	//	}

	//}

}

void NeighborSearchGridPBC::updateNeighborLists() {

	unsigned int nChangedParticles=positionBuffer->getNumberOfChangedElements();
		//std::cout << nChangedParticles << std::endl;

	// prune the list of neighbors of particles (because of particles which have moved apart)

	for (unsigned int i = 0; i<nChangedParticles; i++) {

		unsigned int particleIIndex=positionBuffer->getChangedElementIndex(i);
		Particle* particleI=particleIndex->getObject(particleIIndex);
		Vector const& positionI=(*positionBuffer)[particleIIndex];

		SBCContainerIndex<Particle*>* neighborIndexI=(*neighborIndexBuffer)[particleIIndex];

		unsigned int j=0;

		while (j<neighborIndexI->size()) {

			Particle* particleJ=neighborIndexI->getObject(j);
			unsigned int particleJIndex=particleIndex->getIndex(particleJ);
			Vector const& positionJ = (*positionBuffer)[particleJIndex];

			if ((positionJ-positionI).norm2()<=(*cutoffDistanceSquared))
				j++; // keep this neighbor
			else {

				neighborIndexI->eraseIndex(j);
				(*neighborIndexBuffer)[particleJIndex]->eraseObject(particleI);

			}

		}

	}

	// update the grid

	for (unsigned int i=0; i<nChangedParticles; i++) {

		unsigned int particleIIndex=positionBuffer->getChangedElementIndex(i);
		Particle* particleI=particleIndex->getObject(particleIIndex);
		Vector const& positionI=(*positionBuffer)[particleIIndex];

		GridKey& gridKeyIOld=(*gridKeyVector)[particleIIndex];
		Vector cellCoordinates;
		getCellCoordinates(positionI, cellCoordinates);
		GridKey gridKeyINew((int)cellCoordinates[0], (int)cellCoordinates[1], (int)cellCoordinates[2]);

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

	// update neighbor lists

	updateLists();

	positionBuffer->flush();

}

void NeighborSearchGridPBC::print(unsigned int offset) const {

	unsigned int nParticles=particleIndex->size();

	std::cout << "Neighboring particles: " << std::endl;

	unsigned int nPairs=0;
	for (unsigned int i=0;i<nParticles;i++) {

		SBCContainerIndex<Particle*>* neighborIndexI=(*neighborIndexBuffer)[i];

		if (neighborIndexI->size() > 0)
		//std::cout << "Neighbor pair: " << i << std::endl;

		for (unsigned int j=0;j<neighborIndexI->size();j++) {

			std::cout << "Neighbor pair: " << i << " " << particleIndex->getIndex(neighborIndexI->getObject(j)) << std::endl;
			//std::cout <<   particleIndex->getIndex(neighborIndexI->getObject(j)) << ", ";
			nPairs++;

		}
		//std::cout << std::endl;
	}

	std::cout << "Number of grid cells: " << grid->size() << std::endl;
	std::cout << "Number of pairs: " << nPairs << std::endl << std::endl;

	std::cout << "Number of changed: " << positionBuffer->getNumberOfChangedElements() << std::endl << std::endl;


}



void NeighborSearchGridPBC::flushPositionBuffer() {

	positionBuffer->flush();
}


bool NeighborSearchGridPBC::cellInTheBox(GridKey& gridKey) {

	return
		((gridKey.x+1)*(*cutoffDistance) <= boxPBC->i[0].i[1]) && (gridKey.x*(*cutoffDistance) >= boxPBC->i[0].i[0])
		&&
		((gridKey.y+1)*(*cutoffDistance) <= boxPBC->i[1].i[1]) && (gridKey.y*(*cutoffDistance) >= boxPBC->i[1].i[0])
		&&
		((gridKey.z+1)*(*cutoffDistance) <= boxPBC->i[2].i[1]) && (gridKey.z*(*cutoffDistance) >= boxPBC->i[2].i[0]);

}

bool NeighborSearchGridPBC::cellCrossesTheBox(GridKey& gridKey) {

	return
		((gridKey.x + 1)*(*cutoffDistance) >= boxPBC->i[0].i[1]) && (gridKey.x*(*cutoffDistance) <= boxPBC->i[0].i[1])
		||
		((gridKey.x + 1)*(*cutoffDistance) >= boxPBC->i[0].i[0]) && (gridKey.x*(*cutoffDistance) <= boxPBC->i[0].i[0])
		||
		((gridKey.y + 1)*(*cutoffDistance) >= boxPBC->i[1].i[1]) && (gridKey.y*(*cutoffDistance) <= boxPBC->i[1].i[1])
		||
		((gridKey.y + 1)*(*cutoffDistance) >= boxPBC->i[1].i[0]) && (gridKey.y*(*cutoffDistance) <= boxPBC->i[1].i[0])
		||
		((gridKey.z + 1)*(*cutoffDistance) >= boxPBC->i[2].i[1]) && (gridKey.z*(*cutoffDistance) <= boxPBC->i[2].i[1])
		||
		((gridKey.z + 1)*(*cutoffDistance) >= boxPBC->i[3].i[0]) && (gridKey.z*(*cutoffDistance) <= boxPBC->i[2].i[0]);

}

bool NeighborSearchGridPBC::cellOutsideTheBox(GridKey& gridKey) {

	return

		((gridKey.x + 1)*(*cutoffDistance) <= boxPBC->i[0].i[0]) || (gridKey.x*(*cutoffDistance) >= boxPBC->i[0].i[1])
		||
		((gridKey.y + 1)*(*cutoffDistance) <= boxPBC->i[1].i[0]) || (gridKey.y*(*cutoffDistance) >= boxPBC->i[1].i[1])
		||
		((gridKey.z + 1)*(*cutoffDistance) <= boxPBC->i[2].i[0]) || (gridKey.z*(*cutoffDistance) >= boxPBC->i[2].i[1]);

}

void NeighborSearchGridPBC::getMirrorCellKey(GridKey& gridKey, GridKey& gridKeyMirror) {

	gridKeyMirror = gridKey;

	if (boxPBC->i[0].i[1] < (gridKey.x + 1)*(*cutoffDistance)) {
		// most right x
		gridKeyMirror.x = (int)floor((boxPBC->i[0].i[0])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.x)*(*cutoffDistance) - boxPBC->i[0].i[0]) > 0.1)
			gridKeyMirror.x++;
	}
	else if (gridKey.x*(*cutoffDistance) < boxPBC->i[0].i[0]) {
		// most left x
		gridKeyMirror.x = (int)floor((boxPBC->i[0].i[1])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.x+1)*(*cutoffDistance) - boxPBC->i[0].i[1]) > 0.1)
			gridKeyMirror.x--;
	}

	// y
	if (boxPBC->i[1].i[1] < (gridKey.y + 1)*(*cutoffDistance)) {
		// most right y
		gridKeyMirror.y = (int)floor((boxPBC->i[1].i[0])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.y)*(*cutoffDistance) - boxPBC->i[1].i[0]) > 0.1)
			gridKeyMirror.y++;
	}
	else if (gridKey.y*(*cutoffDistance) < boxPBC->i[1].i[0]) {
		// most left y
		gridKeyMirror.y = (int)floor((boxPBC->i[1].i[1])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.y + 1)*(*cutoffDistance) - boxPBC->i[1].i[1]) > 0.1)
			gridKeyMirror.y--;
	}

	// z
	if (boxPBC->i[2].i[1] < (gridKey.z + 1)*(*cutoffDistance)) {
		// most right z
		gridKeyMirror.z = (int)floor((boxPBC->i[2].i[0])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.z)*(*cutoffDistance) - boxPBC->i[2].i[0]) > 0.1)
			gridKeyMirror.z++;
	}
	else if (gridKey.z*(*cutoffDistance) < boxPBC->i[2].i[0]) {
		// most left z
		gridKeyMirror.z = (int)floor((boxPBC->i[2].i[1])*(*inverseCutoffDistance));
		if (fabs((gridKeyMirror.z + 1)*(*cutoffDistance) - boxPBC->i[2].i[1]) > 0.1)
			gridKeyMirror.z--;
	}

	// x
	//if (boxPBC->i[0].i[1] <= (gridKey.x+1)*(*cutoffDistance)) {
	//	// most right x
	//	gridKeyMirror.x = (int)floor((boxPBC->i[0].i[0])*(*inverseCutoffDistance));
	//	if (fabs((gridKeyMirror.x+1)*(*cutoffDistance) - boxPBC->i[0].i[1]) > 0.1)
	//		gridKeyMirror.x++;
	//}
	//else if (gridKey.x*(*cutoffDistance) <= boxPBC->i[0].i[0]) {
	//	// most left x
	//	gridKeyMirror.x = (int)floor((boxPBC->i[0].i[1])*(*inverseCutoffDistance));
	//	if (fabs(gridKeyMirror.x*(*cutoffDistance) - boxPBC->i[0].i[0]) > 0.1)
	//		gridKeyMirror.x--;
	//}

	//// y
	//if (boxPBC->i[1].i[1] <= (gridKey.y+1)*(*cutoffDistance)) {
	//	// most right y
	//	gridKeyMirror.y = (int)floor((boxPBC->i[1].i[0])*(*inverseCutoffDistance));
	//	if (fabs((gridKeyMirror.y+1)*(*cutoffDistance) - boxPBC->i[1].i[1]) > 0.1)
	//		gridKeyMirror.y++;
	//}
	//else if (gridKey.y*(*cutoffDistance) <= boxPBC->i[1].i[0]) {
	//	// most left y
	//	gridKeyMirror.y = (int)floor((boxPBC->i[1].i[1])*(*inverseCutoffDistance));
	//	if (fabs(gridKeyMirror.y*(*cutoffDistance) - boxPBC->i[1].i[0]) > 0.1)
	//		gridKeyMirror.y--;
	//}

	//// z
	//if (boxPBC->i[2].i[1] <= (gridKey.z+1)*(*cutoffDistance)) {
	//	// most right z
	//	gridKeyMirror.z = (int)floor((boxPBC->i[2].i[0])*(*inverseCutoffDistance));
	//	if (fabs((gridKeyMirror.z+1)*(*cutoffDistance) - boxPBC->i[2].i[1]) > 0.1)
	//		gridKeyMirror.z++;
	//}
	//else if (gridKey.z*(*cutoffDistance) <= boxPBC->i[2].i[0]) {
	//	// most left z
	//	gridKeyMirror.z = (int)floor((boxPBC->i[2].i[1])*(*inverseCutoffDistance));
	//	if (fabs(gridKeyMirror.z*(*cutoffDistance) - boxPBC->i[2].i[0]) > 0.1)
	//		gridKeyMirror.z--;
	//}

}

//void NeighborSearchGridPBC::getMirrorCellKey(GridKey& gridKey, std::vector<GridKey>& gridKeyMirror) {
//
//	std::vector<int> possibleX;
//	std::vector<int> possibleY;
//	std::vector<int> possibleZ;
//
//	// x
//	if (boxPBC->i[0].i[1] <= (gridKey.x + 1)*(*cutoffDistance)) {
//		// most right x
//		possibleX.push_back((int)floor((boxPBC->i[0].i[0])*(*inverseCutoffDistance)));
//		if (fabs((gridKeyMirror.x + 1)*(*cutoffDistance) - boxPBC->i[0].i[1]) > 0.1)
//			gridKeyMirror.x++;
//	}
//	else if (gridKey.x*(*cutoffDistance) <= boxPBC->i[0].i[0]) {
//		// most left x
//		possibleX.push_back((int)floor((boxPBC->i[0].i[1])*(*inverseCutoffDistance)));
//		if (fabs(gridKeyMirror.x*(*cutoffDistance) - boxPBC->i[0].i[0]) > 0.1)
//			gridKeyMirror.x--;
//	}
//
//	// y
//	if (boxPBC->i[1].i[1] <= (gridKey.y + 1)*(*cutoffDistance)) {
//		// most right y
//		gridKeyMirror.y = (int)floor((boxPBC->i[1].i[0])*(*inverseCutoffDistance));
//		if (fabs((gridKeyMirror.y + 1)*(*cutoffDistance) - boxPBC->i[1].i[1]) > 0.1)
//			gridKeyMirror.y++;
//	}
//	else if (gridKey.y*(*cutoffDistance) <= boxPBC->i[1].i[0]) {
//		// most left y
//		gridKeyMirror.y = (int)floor((boxPBC->i[1].i[1])*(*inverseCutoffDistance));
//		if (fabs(gridKeyMirror.y*(*cutoffDistance) - boxPBC->i[1].i[0]) > 0.1)
//			gridKeyMirror.y--;
//	}
//
//	// z
//	if (boxPBC->i[2].i[1] <= (gridKey.z + 1)*(*cutoffDistance)) {
//		// most right z
//		gridKeyMirror.z = (int)floor((boxPBC->i[2].i[0])*(*inverseCutoffDistance));
//		if (fabs((gridKeyMirror.z + 1)*(*cutoffDistance) - boxPBC->i[2].i[1]) > 0.1)
//			gridKeyMirror.z++;
//	}
//	else if (gridKey.z*(*cutoffDistance) <= boxPBC->i[2].i[0]) {
//		// most left z
//		gridKeyMirror.z = (int)floor((boxPBC->i[2].i[1])*(*inverseCutoffDistance));
//		if (fabs(gridKeyMirror.z*(*cutoffDistance) - boxPBC->i[2].i[0]) > 0.1)
//			gridKeyMirror.z--;
//	}
//
//
//	gridKeyMirror.push_back(gridKey);
//
//}

void NeighborSearchGridPBC::checkCellsAroundSmallBox(GridKey& gridKeyI, int i) {

	std::set<unsigned int> neighbors;

	// loop through all neighbouring cells
	for (int x=gridKeyI.x-1;x<=gridKeyI.x+1;x++) {

		for (int y=gridKeyI.y-1;y<=gridKeyI.y+1;y++) {

			for (int z=gridKeyI.z-1;z<=gridKeyI.z+1;z++) {

				// get the neighbouring cell J
				GridKey gridKeyJ(x,y,z);

				// real cell inside the box
				if (cellInTheBox(gridKeyJ)) {

					GridCell* gridCellJ=0;
					// if cell exists
					if (grid->getValue(gridKeyJ,gridCellJ)) {
						checkCell(gridCellJ, i, neighbors);
					}

				}

				// mirror cell
				else {

					GridKey gridKeyJMirror;
					getMirrorCellKey(gridKeyJ, gridKeyJMirror);

					// check mirror cell
					GridCell* gridCellJMirror=0;
					if (grid->getValue(gridKeyJMirror,gridCellJMirror)) {
						checkCell(gridCellJMirror, i, neighbors);
					}

				}

			} // z loop

		} // y loop

	} // x loop

}

void NeighborSearchGridPBC::checkCellsAroundBigBox(GridKey& gridKeyI, int i) {

	std::set<unsigned int> neighbors;

	// loop through all neighbouring cells
	for (int x=gridKeyI.x-1;x<=gridKeyI.x+1;x++) {

		for (int y=gridKeyI.y-1;y<=gridKeyI.y+1;y++) {

			for (int z=gridKeyI.z-1;z<=gridKeyI.z+1;z++) {

				// get the neighbouring cell J
				GridKey gridKeyJ(x,y,z);

				// real cell inside the box
				if (cellInTheBox(gridKeyJ)) {

					GridCell* gridCellJ=0;
					// if cell exists
					if (grid->getValue(gridKeyJ,gridCellJ))
						checkCell(gridCellJ, i, neighbors);

				}

				// mirror cell
				else {

					GridKey gridKeyJMirror;
					getMirrorCellKey(gridKeyJ, gridKeyJMirror);

					if ((abs(gridKeyJMirror.x-gridKeyI.x) > 1 || abs(gridKeyJMirror.y-gridKeyI.y) > 1) || abs(gridKeyJMirror.z-gridKeyI.z) > 1) {
						// check mirror cell
						GridCell* gridCellJMirror=0;
						if (grid->getValue(gridKeyJMirror,gridCellJMirror))
							checkCell(gridCellJMirror, i, neighbors);

					}

				}

			} // z loop

		} // y loop

	} // x loop

}

void NeighborSearchGridPBC::checkCellsAroundCellByCell() {

	for (SBHashMap<GridKey,GridCell*,GridKeyFunctor,GridKeyComparator>::iterator it = grid->begin(); it!=grid->end(); ++it) {

		GridKey gridKeyI = (*it).getKey();
		GridCell* gridCellI=0;
		grid->getValue(gridKeyI,gridCellI);

		if (!gridCellI)
			continue;

		int inside = 0;

		// loop through 13 neighbouring cells and the cell itself
		for (int x=gridKeyI.x-1;x<=gridKeyI.x+1;x++) {

			for (int y=gridKeyI.y-1;y<=gridKeyI.y+1;y++) {

				for (int z=gridKeyI.z-1;z<=gridKeyI.z+1;z++) {

					inside++;
					if (inside > 14)
						continue;

					// get the neighbouring cell J
					GridKey gridKeyJ(x,y,z);

					// real cell inside the box
					if (cellInTheBox(gridKeyJ)) {

						GridCell* gridCellJ=0;
						// if cell exists
						if (grid->getValue(gridKeyJ,gridCellJ))
							checkCellsCross(gridCellI, gridCellJ);

					}

					// mirror cell
					else {

						GridKey gridKeyJMirror;
						getMirrorCellKey(gridKeyJ, gridKeyJMirror);

						if ((abs(gridKeyJMirror.x-gridKeyI.x) > 1 || abs(gridKeyJMirror.y-gridKeyI.y) > 1) || abs(gridKeyJMirror.z-gridKeyI.z) > 1) {
							// check mirror cell
							GridCell* gridCellJMirror=0;
							if (grid->getValue(gridKeyJMirror,gridCellJMirror))
								checkCellsCross(gridCellI, gridCellJMirror);

						}

					}

				} // z loop

			} // y loop

		} // x loop

	}

}

void NeighborSearchGridPBC::signalPositionUpdate(Particle* particle) {

	unsigned int index = particleIndex->getIndex(particle);
	positionBuffer->setElement(index, particle->getPosition());

}

SBCContainerIndex<Particle*> const* NeighborSearchGridPBC::getNeighborIndex(unsigned int i) {
	return (*neighborIndexBuffer)[i];
}

SBCContainerIndex<Particle*> const* NeighborSearchGridPBC::getNeighborIndex(Particle* structuralParticle) {

	if (!particleIndex->hasIndex(structuralParticle)) return 0;
	return
		(*neighborIndexBuffer)[particleIndex->getIndex(structuralParticle)];

}

void NeighborSearchGridPBC::onPositionChanged(Particle* particle) {

	unsigned int index = particleIndex->getIndex(particle);
	positionBuffer->setElement(index, particle->getPosition());

}
