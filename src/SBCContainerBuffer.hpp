#pragma once


#include "SBCTimeClock.hpp"
#include "SBCContainerVector.hpp"


/// The class SBCContainerBuffer is a template defining an extensible array of elements that may be incrementally changed. It is 
/// a key element used in incremental computations.

template <class T> class SBCContainerBuffer {


public:

	/// \name Constructors and destructors
	//@{

	SBCContainerBuffer() { 
	
		typedVector=new SBVector<T>();
		indexVector=new SBVector<unsigned int>();
		xedniVector=new SBVector<unsigned int>();
		stampVector=new SBVector<SBTime>();
	
		flush();

	}
	
	SBCContainerBuffer(unsigned int i) { 
	
		typedVector=new SBVector<T>(i);
		indexVector=new SBVector<unsigned int>(i);
		xedniVector=new SBVector<unsigned int>(i);
		stampVector=new SBVector<SBTime>(i);
	
		for (unsigned int j=0;j<i;j++) {

			typedVector->push_back(T());
			indexVector->push_back(0);
			xedniVector->push_back(0);
			stampVector->push_back(0);

		}

		flush();

	}
	
	virtual ~SBCContainerBuffer() { 
		
		delete typedVector;
		delete indexVector;
		delete xedniVector;
		delete stampVector;

	}

	//@}

	/// \name Basic operations
	//@{
	
	T const&													operator[](unsigned int i) const { return (*typedVector)[i]; }

	unsigned int												size() const { return typedVector->size(); }
	bool														empty() const { return (typedVector->empty()); }

	//@}

	/// \name Buffer management
	//@{
	
	void														addElement(const T& t) {

		typedVector->push_back(t);
		indexVector->push_back((unsigned int)(typedVector->size()-1));
		xedniVector->push_back((unsigned int)(indexVector->size()-1));
		stampVector->push_back(SBCTimeClock());

	}

	void														removeElement(unsigned int i) {

		if (i>=typedVector->size()) return;

		if ((*stampVector)[i]>=lastFlushTime) {

			// element i was in the list of changed elements, we need to remove it

			if ((*xedniVector)[i]<(unsigned int)(indexVector->size()-1)) {
				
				(*indexVector)[(*xedniVector)[i]]=(*indexVector)[(unsigned int)(indexVector->size()-1)];
				(*xedniVector)[(*indexVector)[(*xedniVector)[i]]]=(*xedniVector)[i];
			
			}
			
			indexVector->pop_back();
			
		}

		if (i<(unsigned int)(typedVector->size()-1)) {

			// element typedVector->size()-1 becomes element i

			if ((*stampVector)[(unsigned int)(typedVector->size()-1)]>=lastFlushTime) {

				// element typedVector->size()-1 is in the list of changed elements, we need to rename it

				(*indexVector)[(*xedniVector)[(unsigned int)(typedVector->size()-1)]]=i;

			}

			(*typedVector)[i]=(*typedVector)[(unsigned int)(typedVector->size()-1)];
			(*xedniVector)[i]=(*xedniVector)[(unsigned int)(xedniVector->size()-1)];
			(*stampVector)[i]=(*stampVector)[(unsigned int)(stampVector->size()-1)];

		}

		typedVector->pop_back();
		xedniVector->pop_back();
		stampVector->pop_back();

	}

	SBTime														getLastFlushTime() const { return lastFlushTime; }

	void														flush() { 

		indexVector->clear();
		lastFlushTime=SBCTimeClock();

	}

	//@}

	/// \name Change operations
	//@{
	
	unsigned int												getNumberOfChangedElements() { return indexVector->size(); }
	unsigned int												getChangedElementIndex(unsigned int i) { return (*indexVector)[i]; }
	bool														getChangedElementFlag(unsigned int i) { return (*stampVector)[i]>=lastFlushTime; }

	void														setElement(unsigned int i, const T& t) {

		if (i>=typedVector->size()) return;

		if ((*typedVector)[i] == t) 
			return;

		(*typedVector)[i]=t;

		if ((*stampVector)[i]<lastFlushTime) {

			indexVector->push_back(i);
			(*xedniVector)[i]=indexVector->size()-1;
			(*stampVector)[i]=SBCTimeClock();

		}

	}

	//@}

	/// \name Debugging
	//@{

	void														print(unsigned int offset=0) const {

		for (unsigned int i=0;i<offset;i++) std::cout << "\t";
		std::cout << "Buffer (last flush time = " << lastFlushTime << "):\n";

		for (unsigned int i=0;i<offset+1;i++) std::cout << "\t";
		std::cout << "Typed vector: ";
		typedVector->print();

		for (unsigned int i=0;i<offset+1;i++) std::cout << "\t";
		std::cout << "Index vector: ";
		indexVector->print();

		for (unsigned int i=0;i<offset+1;i++) std::cout << "\t";
		std::cout << "Xedni vector: ";
		xedniVector->print();

		for (unsigned int i=0;i<offset+1;i++) std::cout << "\t";
		std::cout << "Stamp vector: ";
		stampVector->print();

		// check consistency of stamps

		for (unsigned int i=0;i<stampVector->size();i++) {
			
			bool stampOk=true;

			if ((*stampVector)[i]<lastFlushTime) {
				
				for (unsigned int j=0;j<indexVector->size();j++) if ((*indexVector)[j]==i) stampOk=false;

			}
			else {
				
				stampOk=false;
				for (unsigned int j=0;j<indexVector->size();j++) if ((*indexVector)[j]==i) stampOk=true;

			}

			if (!stampOk) std::cout << "Stamp incorrect for element " << i << "\n";

		}

		// check consistency of indices

		for (unsigned int i=0;i<indexVector->size();i++) {
			
			if ((*xedniVector)[(*indexVector)[i]]!=i) std::cout << "Index vector incorrect at position " << i << "\n";

		}


		std::cout << std::endl;

	}

	unsigned int												getMemoryFootPrint() const  { return typedVector->getMemoryFootPrint()+indexVector->getMemoryFootPrint()+xedniVector->getMemoryFootPrint()+stampVector->getMemoryFootPrint()+sizeof(SBTime); }

	//@}


private:

	SBCContainerVector<T>*										typedVector;															///< The vector of elements
	SBCContainerVector<unsigned int>*							indexVector;															///< The vector of indices of changed elements between two buffer flushes
	SBCContainerVector<unsigned int>*							xedniVector;															///< The vector of reverse indices : indexVector[xedniVector[i]]=i and xedniVector[indexVector[j]]=j
	SBCContainerVector<SBTime>*									stampVector;															///< The vector of time stamps, to determine whether an element has already been marked as having been changed

	SBTime														lastFlushTime;															///< Last flush time

};


#define SBBuffer												SBCContainerBuffer

