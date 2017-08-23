/*************************************************************************\

	SAMSON - Software for Adaptive Modeling and Simulation Of Nanosystems

	Copyright INRIA - NANO-D
	All Rights Reserved.

\**************************************************************************/


#pragma once

#include "SBCContainerVector.hpp"


#include <set>
#include <string>
#include <vector>
#include <assert.h>


template <class KeyType> class SBCContainerHashMapFunctor {

public: 

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const KeyType& key) { 

		unsigned long hash=reinterpret_cast<unsigned long>(key);
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<long> {

public:

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const long& key) {

		unsigned long hash = (unsigned long)key;
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<unsigned long> {

public:

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const unsigned long& key) {

		unsigned long hash = key;
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<int> {

public:

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const int& key) {

		unsigned long hash = (unsigned long)key;
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<unsigned int> {

public:

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const unsigned int& key) {

		unsigned long hash = key;
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<std::string> {

public: 

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const std::string& key) { 

		unsigned long hash=0;
		unsigned long length=(unsigned long)key.size();

		for (unsigned long i=0;i<length;i++) {

			hash+=key[i];
			hash = (hash << 6) + (hash << 16) - hash;

		}

		return hash;

	}

};

template <> class SBCContainerHashMapFunctor<char> {

public: 

	SBCContainerHashMapFunctor() {}
	virtual ~SBCContainerHashMapFunctor() {}

	static unsigned long										hashCode(const char& key) { 

		unsigned long hash=key;
		hash = (hash << 6) + (hash << 16) - hash;

		return hash;

	}

};

template <class KeyType> class SBCContainerHashMapKeyComparator {

public: 

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const KeyType& key1, const KeyType& key2) { 

		return (reinterpret_cast<unsigned long>(key1)==reinterpret_cast<unsigned long>(key2));

	}

};

template <> class SBCContainerHashMapKeyComparator<long> {

public:

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const long& key1, const long& key2) {

		return (key1 == key2);

	}

};

template <> class SBCContainerHashMapKeyComparator<unsigned long> {

public:

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const unsigned long& key1, const unsigned long& key2) {

		return (key1 == key2);

	}

};

template <> class SBCContainerHashMapKeyComparator<int> {

public:

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const int& key1, const int& key2) {

		return (key1 == key2);

	}

};

template <> class SBCContainerHashMapKeyComparator<unsigned int> {

public:

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const unsigned int& key1, const unsigned int& key2) {

		return (key1 == key2);

	}

};

template <> class SBCContainerHashMapKeyComparator<std::string> {

public: 

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const std::string& key1, const std::string& key2) { 

		return (key1.compare(key2)==0);

	}

};

template <> class SBCContainerHashMapKeyComparator<char> {

public: 

	SBCContainerHashMapKeyComparator() {}
	virtual ~SBCContainerHashMapKeyComparator() {}

	static bool													equal(const char& key1, const char& key2) { 

		return (key1==key2);

	}

};

template <class KeyType, class ValueType> class SBCContainerHashMapEntry {

public:

	SBCContainerHashMapEntry(const KeyType& k, const ValueType& v, SBCContainerHashMapEntry* n) { key = k; value = v; next = n; }
	virtual ~SBCContainerHashMapEntry() { if (next != NULL) delete next; }

	KeyType													getKey() const { return key; }
	ValueType												getValue() const { return value; }
	SBCContainerHashMapEntry*								getNext() const { return next; }

	void													setValue(const ValueType& v) { value = v; }
	void													setNext(SBCContainerHashMapEntry* n) { next = n; }

	unsigned int											getMemoryFootprint() const { return sizeof(KeyType)+sizeof(ValueType)+sizeof(SBCContainerHashMapEntry*); }

private:

	KeyType													key;
	ValueType												value;
	SBCContainerHashMapEntry*								next;

};


/// \brief SBCContainerHashMap is a template class which implements a hashtable-based map.
///
/// \param KeyType The type of keys
/// \param ValueType The type of values
/// \param HashMapFunctor The hash map functor, used to hash keys (by default, this is a SBCContainerHashMapFunctor<KeyType>)
/// \param KeyComparator The key comparator, used to compare keys (by default, this is a SBCContainerHashMapKeyComparator<KeyType>)
///
/// SBCContainerHashMap is a template class which implements a hashtable-based map. 
/// This container is used for example to index objects (see SBCContainerIndex).

template <class KeyType, class ValueType, class HashMapFunctor=SBCContainerHashMapFunctor<KeyType>, class KeyComparator=SBCContainerHashMapKeyComparator<KeyType> > class SBCContainerHashMap {

public:

	class iterator {

		friend class SBCContainerHashMap;

	public:

		iterator() { entry=0;hashMapIndex=0;hashMap=0; }
		iterator(SBCContainerHashMapEntry<KeyType,ValueType>* e, unsigned int i, SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>* h) { entry=e;hashMapIndex=i;hashMap=h; }
		iterator(const iterator& i) { entry=i.entry;hashMapIndex=i.hashMapIndex;hashMap=i.hashMap; }

		virtual ~iterator() {}

		iterator&													operator++() { 
			
			if (entry==0) return *this;	// the iterator is not referencing any entry

			entry=entry->getNext();
			if (entry!=0) return *this;	// the entry was not the last one in the list attached to the vector element

			if (entry==0) { // the entry was the last one in the list, find the next list if it exists

				hashMapIndex++;
				unsigned int hashMapSize=hashMap->size();

				while (hashMapIndex<hashMapSize) {

					SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=hashMap->getElement(hashMapIndex);

					if (currentEntry!=0) { // we found the next entry
						
						entry=currentEntry;
						return *this;

					}

					hashMapIndex++;

				}

			}

			// no next entry was found, the iterator does not reference any entry

			return *this; 
		
		}

		iterator													operator++(int) { 
			
			iterator i(entry);
			++(*this);
			return i; 
		
		}

		iterator&													operator--() { entry--;return this; }
		iterator													operator--(int) { entry--;return this; }

		iterator&													operator=(const iterator& i) { entry=i.entry;hashMapIndex=i.hashMapIndex;hashMap=i.hashMap;return (*this); }

		bool														operator==(iterator& i) { return entry==i.entry; }
		bool														operator!=(iterator& i) { return entry!=i.entry; }
		bool														operator==(const iterator& i) { return entry==i.entry; }
		bool														operator!=(const iterator& i) { return entry!=i.entry; }

		SBCContainerHashMapEntry<KeyType,ValueType>*				operator->() const { return entry; }
		SBCContainerHashMapEntry<KeyType,ValueType>&				operator*() const { return *entry; }
		
	private:

		SBCContainerHashMapEntry<KeyType,ValueType>*						entry;														///< The poitner to the hash map entry
		unsigned int														hashMapIndex;												///< The index of the vector element in which the map entry is stored
		SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>*	hashMap;													///< The pointer to the hash map vector
	
	};

	/// \name Constructors and destructors
	//@{

	SBCContainerHashMap(unsigned int nElements=0);																						///< Creates a new hash map.
	virtual ~SBCContainerHashMap();																										///< Deletes the hash map.

	//@}

	/// \name Insertion and look-up
	//@{

	void														clear();																///< Erases all entries in the hash map.
	void														insert(const KeyType& key, const ValueType& value);						///< Inserts an element in the hash map.
	bool														erase(const KeyType& key);												///< Attempts to erase an entry with a specific \a key and returns true if successful.

	bool														hasValue(const KeyType& key) const;										///< Returns true if the \a key has an associated value.
	bool														getValue(const KeyType& key, ValueType& value) const;					///< Returns the value associated to the key.
	unsigned int												size() const;															///< The number of entries in the hash map.
	bool														empty() const;															///< Returns true when the hash map contains no entry.

	//@}

	/// \name Iterators
	//@{

	iterator													begin() { 
	
		iterator it(0,0,hashMap);
		if (numberOfEntries==0) return it;

		// find the first entry in the map

		unsigned int hashMapSize=hashMap->size();

		for (unsigned int i=0;i<hashMapSize;i++) {

			SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=hashMap->getElement(i);
			if (currentEntry!=0) { 
				
				it.entry=currentEntry;
				it.hashMapIndex=i;
				return it;

			}

		}

		return it; 
	
	}

	iterator													end() { return iterator(); }

	//@}

	/// \name Debugging
	//@{

	unsigned int												getMemoryFootprint() const;

	//@}

protected:

	/// \name The hash table
	//@{

	unsigned int												numberOfEntries;			///< The number of entries in the hash map.

	SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>*		hashMap;					///< The hash table

	//@}

};

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::SBCContainerHashMap(unsigned int nElements) {

	numberOfEntries=0;

	// allocate memory for the hash map

	unsigned int hashMapSize=nElements*3+11;
	hashMap=new SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>(hashMapSize);
	for (unsigned int i=0;i<hashMapSize;i++) hashMap->push_back(0);

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::~SBCContainerHashMap() {

	clear();
	delete hashMap;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> void SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::insert(const KeyType& key, const ValueType& value) {

	if (numberOfEntries*3+11>hashMap->size()) {

		// we would like the vector to contain at least numberOfEntries*3+11 elements.
		// If that's not the case, we create a new one and rehash the keys.

		unsigned int oldHashMapSize=hashMap->size();
		unsigned int newHashMapSize=(unsigned int)(2*oldHashMapSize);
		SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>* oldHashMap=hashMap;
		hashMap=new SBCContainerVector<SBCContainerHashMapEntry<KeyType,ValueType>*>(newHashMapSize);
		for (unsigned int i=0;i<newHashMapSize;i++) hashMap->push_back(0);

		// rehash

		for (unsigned int i=0;i<oldHashMapSize;i++) {

			while ((*oldHashMap)[i]!=0) {

				// pop the old hash map entry

				SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=(*oldHashMap)[i];
				(*oldHashMap)[i]=(*oldHashMap)[i]->getNext();

				// insert the entry into the new hash map

				unsigned int hashIndex=HashMapFunctor::hashCode(currentEntry->getKey())%hashMap->size();
				currentEntry->setNext((*hashMap)[hashIndex]);
				(*hashMap)[hashIndex]=currentEntry;

			}

		}

		delete oldHashMap;

	}

	unsigned int hashIndex=HashMapFunctor::hashCode(key)%hashMap->size();
	SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=(*hashMap)[hashIndex];

	while (currentEntry!=NULL) {

		if (KeyComparator::equal(currentEntry->getKey(),key)) {

			currentEntry->setValue(value);
			break;

		}

		currentEntry=currentEntry->getNext();

	}

	if (currentEntry==0) {

		(*hashMap)[hashIndex]=new SBCContainerHashMapEntry<KeyType,ValueType>(key,value,(*hashMap)[hashIndex]);
		numberOfEntries++;

	}

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> bool SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::hasValue(const KeyType& key) const {

	unsigned long hashIndex=HashMapFunctor::hashCode(key)%hashMap->size();

	SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=(*hashMap)[hashIndex];

	while (currentEntry!=NULL) {

		if (KeyComparator::equal(currentEntry->getKey(),key)) return true;
		currentEntry=currentEntry->getNext();

	}

	return false;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> bool SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::getValue(const KeyType& key, ValueType& value) const {

	unsigned long hashIndex=HashMapFunctor::hashCode(key)%hashMap->size();

	SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=(*hashMap)[hashIndex];

	while (currentEntry!=NULL) {

		if (KeyComparator::equal(currentEntry->getKey(),key)) {

			value=currentEntry->getValue();
			return true;

		}

		currentEntry=currentEntry->getNext();

	}

	return false;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> unsigned int SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::size() const {

	return numberOfEntries;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> bool SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::empty() const {

	return numberOfEntries==0;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> void SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::clear() {

	unsigned int vectorSize=hashMap->size();

	for (unsigned int i=0;i<vectorSize;i++) 
		
		if ((*hashMap)[i]!=0) {

			delete (*hashMap)[i];
			(*hashMap)[i]=0;

		}

	numberOfEntries=0;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> bool SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::erase(const KeyType& key) {

	unsigned int hashIndex=HashMapFunctor::hashCode(key)%hashMap->size();

	SBCContainerHashMapEntry<KeyType,ValueType>* currentEntry=(*hashMap)[hashIndex];
	SBCContainerHashMapEntry<KeyType,ValueType>* previousEntry=(*hashMap)[hashIndex];

	while (currentEntry!=NULL) {

		if (KeyComparator::equal(currentEntry->getKey(),key)) {

			// remove the entry from the list

			if (currentEntry==(*hashMap)[hashIndex]) (*hashMap)[hashIndex]=(*hashMap)[hashIndex]->getNext(); // the entry was at the beginning of the list
			else previousEntry->setNext(currentEntry->getNext());

			currentEntry->setNext(0);
			delete currentEntry;

			numberOfEntries--;

			return true;

		}

		if (currentEntry!=(*hashMap)[hashIndex]) previousEntry=previousEntry->getNext();
		currentEntry=currentEntry->getNext();

	}

	return false;

}

template <class KeyType, class ValueType, class HashMapFunctor, class KeyComparator> unsigned int SBCContainerHashMap<KeyType,ValueType,HashMapFunctor,KeyComparator>::getMemoryFootprint() const {

	return sizeof(unsigned int)+hashMap->getMemoryFootPrint();

}


#define SBHashMap												SBCContainerHashMap

