#pragma once


#include <iostream>


/// The class SBCContainerVector is a template defining an extensible array.

template <class T> class SBCContainerVector {

public:

	class iterator {

	public:

		iterator() { t=0; }
		iterator(T* tt) { t=tt; }
		iterator(const iterator& i) { t=i.t; }

		virtual ~iterator() {}

		iterator*													operator++() { t++;return this; }
		iterator*													operator++(int) { t++;return this; }
		iterator*													operator--() { t--;return this; }
		iterator*													operator--(int) { t--;return this; }

		iterator&													operator=(T* tt) { t=tt;return (*this); }
		iterator&													operator=(const iterator& i) { t=i.t;return (*this); }

		bool														operator==(iterator& i) const { return t==i.t; }
		bool														operator!=(iterator& i) const { return t!=i.t; }
		bool														operator==(const iterator& i) const { return t==i.t; }
		bool														operator!=(const iterator& i) const { return t!=i.t; }

		T															operator*() const { return *t; }
		
		T*															getPointer() const { return t; }

	private:

		T*															t;

	};

	class reverse_iterator {

	public:

		reverse_iterator() { t=0; }
		reverse_iterator(T* tt) { t=tt; }
		reverse_iterator(const reverse_iterator& i) { t=i.t; }

		virtual ~reverse_iterator() {}

		reverse_iterator*											operator++() { t--;return this; }
		reverse_iterator*											operator++(int) { t--;return this; }
		reverse_iterator*											operator--() { t++;return this; }
		reverse_iterator*											operator--(int) { t++;return this; }
		
		reverse_iterator&											operator=(T* tt) { t=tt;return (*this); }
		reverse_iterator&											operator=(const reverse_iterator& i) { t=i.t;return (*this); }

		bool														operator==(reverse_iterator& i) const { return t==i.t; }
		bool														operator!=(reverse_iterator& i) const { return t!=i.t; }
		bool														operator==(const reverse_iterator& i) const { return t==i.t; }
		bool														operator!=(const reverse_iterator& i) const { return t!=i.t; }

		T															operator*() const { return *t; }
		
		T*															getPointer() const { return t; }

	private:

		T*															t;

	};

public:

	SBCContainerVector() { firstUnusedIndex=0;currentSize=0;array=0; }
	SBCContainerVector(unsigned int i) { firstUnusedIndex=0;currentSize=i;array=new T[i]; }

	virtual ~SBCContainerVector() { 
		
		if (array) delete [] array; 
	
	}

	T&															getElement(unsigned int i) const { return array[i]; }
	T&															operator[](unsigned int i) const { return array[i]; }
	T&															back() const { return array[firstUnusedIndex-1]; }
	T*															getPointer(unsigned int i) const { return &array[i]; }

	T&															getElement(unsigned int i) { return array[i]; }
	T&															operator[](unsigned int i) { return array[i]; }
	T&															back() { return array[firstUnusedIndex-1]; }
	T*															getPointer(unsigned int i) { return &array[i]; }

	void														push_back(const T& t) {

		if (firstUnusedIndex>=currentSize) {

			unsigned int newSize=(unsigned int)(1.1*currentSize+1);
			T* newArray=new T[newSize];
			
			for (unsigned int i=0;i<currentSize;i++) newArray[i]=array[i];

			delete [] array;
			
			array=newArray;
			currentSize=newSize;

		}

		array[firstUnusedIndex]=t;
		firstUnusedIndex++;

	}

	void														pop_back() {

		if (firstUnusedIndex>0) firstUnusedIndex--;

	}

	bool														empty() const { return (firstUnusedIndex==0); }
	void														clear() { firstUnusedIndex=0; }

	unsigned int												size() const { return firstUnusedIndex; }

	iterator													begin() const { return iterator(array); }
	iterator													end() const { return iterator(array+firstUnusedIndex); }

	reverse_iterator											rbegin() const { return reverse_iterator(array + firstUnusedIndex - 1); }
	reverse_iterator											rend() const { return reverse_iterator(array - 1); }

	void														print(unsigned int offset=0) const {

		for (unsigned int i=0;i<offset;i++) std::cout << "\t";
		for (unsigned int i=0;i<size();i++) std::cout << array[i] << "\t";
		std::cout << "\n";

	}

	unsigned int												getMemoryFootPrint() const { return sizeof(T*)+2*sizeof(unsigned int)+currentSize*sizeof(T); }


private:

	T*															array;
	unsigned int												currentSize;
	unsigned int												firstUnusedIndex;

};


#define SBVector												SBCContainerVector
