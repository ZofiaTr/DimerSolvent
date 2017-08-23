/*************************************************************************\

	SAMSON - Software for Adaptive Modeling and Simulation Of Nanosystems

	Copyright INRIA - NANO-D
	All Rights Reserved.

\**************************************************************************/


#pragma once


#include "SBCContainerVector.hpp" 
#include "SBCContainerHashMap.hpp" 

#include <set>
#include <vector>


/// SBCContainerIndex is a template class used to efficiently index objects. 

/// Fast indexing of objects is required in SAMSON because objects may operate on lists of other objects (e.g. a collision
/// detection solver operates on a list of structural model components). An indexer provides a way to map the set of n objects to a set of ordered
/// indices, from 0 to n-1, which can then be used for fast access to information associated to the mapped objects. 
///
/// The set of indices remains contiguous even when an object is removed from the indexer. To achieve this, when an object with index
/// \a objectIndex is removed, the last object in the set of indexed objects is reindexed and is given the index \a objectIndex (unless
/// the last object was the one removed).
///
/// Internally, SBCContainerIndex uses a <em>hash map</em> (by default SBCContainerHashMap) to assign indices and a \e vector (by defaultSBCContainerVector) to store
/// objects. As a result, retrieving indices takes an almost constant time, while retrieving objects takes constant time.
/// 
/// \shortname{SBIndex}

template <class ObjectType, class HashMap=SBCContainerHashMap<ObjectType,unsigned int>, class Vector=SBCContainerVector<ObjectType> > class SBCContainerIndex {

public:

	/// \name Constructors and destructors
	//@{

	SBCContainerIndex();																												///< Creates an indexer.
	SBCContainerIndex(unsigned int initialSize);																						///< Creates an indexer with a default pre-allocated size.
	SBCContainerIndex(std::set<ObjectType> &objectset);																					///< Indexes a set of objects.
	SBCContainerIndex(std::vector<ObjectType> &objectVector);																			///< Indexes a vector of objects.
	virtual ~SBCContainerIndex();

	//@}

	/// \name Indexing
	//@{

	void														clear();																///< Clears the index

	unsigned int												push_back(const ObjectType& object);									///< Adds an object in the indexer and returns the index of the object.
	unsigned int												pop_back();																///< Adds an object in the indexer and returns the index of the object.

	unsigned int												eraseObject(const ObjectType& object);									///< Erases the \a object from the indexer.
	unsigned int												eraseIndex(unsigned int objectIndex);									///< Erases object \a objectIndex from the indexer.

	unsigned int												size() const;															///< Returns the number of indexed objects.

	bool														hasIndex(const ObjectType& object) const;								///< Returns true if the \a object has an index.
	unsigned int												getIndex(const ObjectType& object) const;								///< Returns the index associated to the \a object.
	bool														getIndex(const ObjectType& object, unsigned int& index) const;			///< Returns the index associated to the \a object.

	ObjectType													getObject(unsigned int index) const;									///< Returns the object associated to the \a index.

	//@}

	/// \name Debugging
	//@{

	void														print() const;

	//@}

protected:

	/// \name The hash table
	//@{

	HashMap*													indexMap;
	Vector*														objectVector;

	//@}

};

template <class ObjectType, class HashMap, class Vector> SBCContainerIndex<ObjectType,HashMap,Vector>::SBCContainerIndex() {

	indexMap=new HashMap();
	objectVector=new Vector();

}

template <class ObjectType, class HashMap, class Vector> SBCContainerIndex<ObjectType,HashMap,Vector>::SBCContainerIndex(unsigned int initialSize) {

	indexMap=new HashMap(initialSize);
	objectVector=new Vector(initialSize);

}

template <class ObjectType, class HashMap, class Vector> SBCContainerIndex<ObjectType,HashMap,Vector>::SBCContainerIndex(std::vector<ObjectType> &objects) {

	unsigned int nObjects=(unsigned int)objects.size();
	indexMap=new HashMap(nObjects);
	objectVector=new Vector(nObjects);

	// insert entries into the hash table

	unsigned int currentIndex=0;

	for (typename std::vector<ObjectType>::iterator i=objects.begin();i!=objects.end();i++) {

		indexMap->insert(*i,currentIndex);
		objectVector->push_back(*i);
		currentIndex++;

	}

}

template <class ObjectType, class HashMap, class Vector> SBCContainerIndex<ObjectType,HashMap,Vector>::SBCContainerIndex(std::set<ObjectType> &objects) {

	unsigned int nObjects=(unsigned int)objects.size();
	indexMap=new HashMap(nObjects);
	objectVector=new Vector(nObjects);

	// insert entries into the hash table

	unsigned int currentIndex=0;

	for (typename std::set<ObjectType>::iterator i=objects.begin();i!=objects.end();i++) {

		indexMap->insert(*i,currentIndex);
		objectVector->push_back(*i);
		currentIndex++;

	}

}

template <class ObjectType, class HashMap, class Vector> SBCContainerIndex<ObjectType,HashMap,Vector>::~SBCContainerIndex() {

	delete indexMap;
	delete objectVector;

}

template <class ObjectType, class HashMap, class Vector> void SBCContainerIndex<ObjectType,HashMap,Vector>::clear() {

	indexMap->clear();
	objectVector->clear();

}

template <class ObjectType, class HashMap, class Vector> bool SBCContainerIndex<ObjectType,HashMap,Vector>::hasIndex(const ObjectType& object) const {

	return indexMap->hasValue(object);

}

template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::getIndex(const ObjectType& object) const {

	unsigned int index=-1;
	indexMap->getValue(object,index);
	return index;

}

template <class ObjectType, class HashMap, class Vector> bool SBCContainerIndex<ObjectType,HashMap,Vector>::getIndex(const ObjectType& object, unsigned int& index) const {

	return indexMap->getValue(object,index);

}

template <class ObjectType, class HashMap, class Vector> ObjectType SBCContainerIndex<ObjectType,HashMap,Vector>::getObject(unsigned int index) const {

	return (*objectVector)[index];

}

template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::push_back(const ObjectType& object) {

	unsigned int objectIndex=-1;
	if (indexMap->getValue(object,objectIndex)) 
		return objectIndex;

	unsigned int newIndex=objectVector->size();
	objectVector->push_back(object);
	indexMap->insert(object,newIndex);
	return newIndex;

}

template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::pop_back() {

	if (objectVector->size()==0) return -1;
	unsigned int removedIndex=objectVector->size()-1;

	eraseIndex(removedIndex);

	return removedIndex;

}

/// The function returns an unsigned int with the following meaning:
///    - if the \a object was not present in the index, the returned value is the number of indexed objects.
///    - else, the returned value is the index of the object that was removed, in [0,size-1].
/// Note that these two returned values coincide when the client (unsuccessfully) attempts to erase object with index \e size.


template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::eraseObject(const ObjectType& object) {

	unsigned int objectIndex=-1;
	if (!indexMap->getValue(object,objectIndex)) return objectVector->size();
	return eraseIndex(objectIndex);

}

/// The function returns an unsigned int with the following meaning:
///    - if the \a object was not present in the index, the returned value is the number of indexed objects.
///    - else, the returned value is the index of the object that was removed, in [0,size-1].
/// Note that these two returned values coincide when the client (unsuccessfully) attempts to erase object with index \e size.

template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::eraseIndex(unsigned int objectIndex) {

	if (objectIndex>=objectVector->size()) return objectVector->size();

	// remove the object from the index map

	indexMap->erase((*objectVector)[objectIndex]);
	
	// remove the object from the object vector. Since we want a contiguous array of indices, 
	// the last object is reindexed: it is given the newly available index (unless the last object is the one being erased).

	if (objectIndex!=objectVector->size()-1) {

		(*objectVector)[objectIndex]=(*objectVector)[objectVector->size()-1];
		indexMap->insert((*objectVector)[objectIndex],objectIndex);

	}

	objectVector->pop_back();

	return objectIndex;

}

template <class ObjectType, class HashMap, class Vector> unsigned int SBCContainerIndex<ObjectType,HashMap,Vector>::size() const {

	return objectVector->size();

}

template <class ObjectType, class HashMap, class Vector> void SBCContainerIndex<ObjectType,HashMap,Vector>::print() const {

	std::cout << "Container index: " << this << std::endl;

	if (objectVector->size()==0) std::cout << "Empty.";

	std::cout << std::endl;

	for (unsigned int i=0;i<objectVector->size();i++) {

		std::cout << i << "\t" << (*objectVector)[i] << "\n";

	}

	std::cout << std::endl;

}


#define SBIndex													SBCContainerIndex
