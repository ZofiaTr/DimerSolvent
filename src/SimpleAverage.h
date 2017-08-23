//#pragma once
#include<iostream>

//Stephane old average class - block averaging with constant blocksize 

template <class T> class SimpleAverage {

public:

	SimpleAverage() {


		blockSize=1000;
		currentBlockIndex=0;				
		currentIndexInCurrentBlock=0;		
		
		currentBlockSum=T(0);				
		sumOfBlockAverages=T(0);				
		
		}


	SimpleAverage(unsigned int blockSize) {

		this->blockSize=blockSize;
		currentBlockIndex=0;				
		currentIndexInCurrentBlock=0;		
		
		currentBlockSum=T(0);				
		sumOfBlockAverages=T(0);				
		
		}

	~SimpleAverage() {
		

	}

	void														addSample(const T& t) {

		// add the sample to the current block

		currentBlockSum+=t;
		currentIndexInCurrentBlock++;

		if (currentIndexInCurrentBlock==blockSize) { // the block is full
			

			// add the block average to the sum of block averages
			currentBlockAverage= currentBlockSum/blockSize;
			sumOfBlockAverages+=currentBlockSum/blockSize;
			
			// begin a new block

			currentBlockSum=T(0);
			currentIndexInCurrentBlock=0;
			currentBlockIndex++;

		}

	}

	T															getAverage() {

		// the average of completed blocks is sumOfBlockAverages/currentBlockIndex (which represents blockSize*currentBlockIndex samples)
		// the average on the current block is currentBlockSum/currentIndexInCurrentBlock (which represents currentIndexInCurrentBlock samples)
		// the current average is thus : 
		// (sumOfBlockAverages/currentBlockIndex*blockSize*currentBlockIndex+currentBlockSum/currentIndexInCurrentBlock*currentIndexInCurrentBlock) / (blockSize*currentBlockIndex+currentIndexInCurrentBlock)

		unsigned int numberOfSamples=blockSize*currentBlockIndex+currentIndexInCurrentBlock;
		if (numberOfSamples==0) return T(0);

		T currentAverage = currentBlockSum/numberOfSamples;

		if (currentBlockIndex>0) {

			currentAverage	+= sumOfBlockAverages/currentBlockIndex;
			currentAverage	-= (currentIndexInCurrentBlock*sumOfBlockAverages)/numberOfSamples/currentBlockIndex;

		}

		return currentAverage;

	}

	T															getCurrentBlockAverage(){
		return currentBlockAverage;
	
	}

	unsigned int												getCurrentIndex(){
	
	
		return currentIndexInCurrentBlock;
	
	}

	void														clean(){

		currentBlockIndex=0;				
		currentIndexInCurrentBlock=0;		
		
		currentBlockSum=T(0);				
		sumOfBlockAverages=T(0);
	
	
	}

	void restartAverage(T startingValue, double numberOfTimeSteps){
	
		sumOfBlockAverages = startingValue;
		currentBlockIndex = numberOfTimeSteps / blockSize;
	}

	
	

private:

	unsigned int												blockSize;																///< The size of a period on which to perform a partial average
	unsigned int												currentBlockIndex;														///< The index of the current averaging period
	unsigned int												currentIndexInCurrentBlock;												///< The current position index in the current averaging period

	T															currentBlockSum;														///< The current block sum
	T															sumOfBlockAverages;														///< The sum of block averages
	T															currentBlockAverage;

			
};





//#include<iostream>
//
//
//template <class T> class SimpleAverage {
//
//public:
//
//	SimpleAverage() {
//
//
//		
//		currentIndex=0;		
//		currentSum=T(0);				
//		
//		
//		}
//	
//	
//	~SimpleAverage() {
//		
//
//	}
//
//	void														addSample(const T& t) {
//
//		// add the sample
//
//		currentSum+=t;
//		currentIndex++;
//			
//	}
//
//	
//
//	T															getAverage(){
//
//		if(currentIndex != 0)		return currentSum/currentIndex;
//		else return T(0);
//	
//	}
//
//	void														clean(){
//
//		currentIndex=0;		
//		currentSum=T(0);
//
//	
//	
//	}
//
//	
//	
//
//private:
//
//
//
//	unsigned int												currentIndex;		
//	T															currentSum;
//
//			
//};
//
//
//
//
