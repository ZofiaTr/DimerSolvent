#pragma once

#include<iostream>
#include"math.h"
#include <cstdlib> 
#include<fstream>




class Vector {

public:

	Vector() { v[0] = v[1] = v[2] = 0.0; }
	Vector(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
	~Vector() {}

	double& operator[](int index) { return v[index]; };
	const double& operator[](int index) const { return v[index]; };

	Vector operator-(const Vector& u) const { return Vector(v[0] - u.v[0], v[1] - u.v[1], v[2] - u.v[2]); }
	bool operator == (const Vector& u) const { return (v[0] == u.v[0]) && (v[1] == u.v[1]) && (v[2] == u.v[2]); }
	double norm2() const { return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
	void print() const { std::cout << v[0] << " " << v[1] << " " << v[2] << "\n"; }
	double v[3];

};

class Interval {

public:

	Interval() { i[0] = i[1] = 0.0; }
	Interval(double x) { i[0] = x; i[1] = x; }
	Interval(double x, double y) { i[0] = x; i[1] = y; }
	~Interval() {}

	double& operator[](int index) { return i[index]; };
	const double& operator[](int index) const { return i[index]; };

	Interval operator-(const Interval& u) const { return Interval(i[0] - u.i[0], i[1] - u.i[1]); }

	double i[2];

};


class IAVector {

public:

	IAVector() { i[0] = i[1] = i[2] = Interval(0.0, 0.0); }
	IAVector(double x, double y, double z) { i[0] = Interval(x); i[1] = Interval(y); i[2] = Interval(z); }
	IAVector(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) { i[0] = Interval(xmin, xmax); i[1] = Interval(ymin, ymax); i[2] = Interval(zmin, zmax); }
	IAVector(Interval x, Interval y, Interval z) { i[0] = x; i[1] = y; i[2] = z; }
	~IAVector() {}

	Interval& operator[](int index) { return i[index]; };
	const Interval& operator[](int index) const { return i[index]; };

	IAVector operator-(const IAVector& u) const { return IAVector(i[0] - u.i[0], i[1] - u.i[1], i[2] - u.i[2]); }

	Interval i[3];

};

class Particle{

	// will contain the x coordinates of positions and momenta
public:
	// constructors
	Particle(){

		itsPositionX=0;
		itsMomentumX=0;
		itsPositionY=0;
		itsMomentumY=0;
		itsPositionZ = 0;
		itsMomentumZ = 0;


		itsPositionXOld = 0.0;
		itsPositionYOld = 0.0;
		itsPositionZOld = 0.0;
		itsMomentumXOld = 0.0;
		itsMomentumYOld = 0.0;
		itsMomentumZOld = 0.0;
		
		itsRX = 0;
		itsRY = 0;
		itsRZ = 0;

		itsMass=1.0;

		//std::cout << "\n\n\n\n Its mass = " << itsMass << std::endl;
		itsFx=0;
		itsFy=0;
		itsFz = 0;

	//	itsFcomingFromInteractions=0;
		active=1;

		itsepsr=0.0;
		itsepsf=0.0;

		type=0;
		index = 0;

		itsGOldx = 0;
		itsGx = 0;
		itsGOldy = 0;
		itsGy = 0;
		itsGOldz = 0;
		itsGz = 0;


		itsFxOld = 0;
		itsFyOld = 0;
		itsFzOld = 0;
	}

	Particle(double px, double py, double pz, double qx, double qy, double qz, double m, double er, double ef){
		
		itsMomentumX = px, itsPositionX = qx, itsMomentumY = py, itsPositionY = qy, itsMomentumZ = pz, itsPositionZ = qz,
	
		itsMass=m;
		itsFx=0;
		itsFy=0;
		itsFz = 0;

		itsFxOld = 0;
		itsFyOld = 0;
		itsFzOld = 0;
		//itsFcomingFromInteractions=0.0;
	
		itsepsr=er;
		itsepsf=ef;

		type=0;

		itsRX = 0;
		itsRY = 0;
		itsRZ = 0; 

		itsPositionXOld = 0.0;
		itsPositionYOld = 0.0;
		itsPositionZOld = 0.0;
		itsMomentumXOld = 0.0;
		itsMomentumYOld = 0.0;
		itsMomentumZOld = 0.0;
		index = 0;

		 itsGOldx=0;
		 itsGx=0;
		 itsGOldy=0;
		 itsGy=0;
		 itsGOldz=0;
		 itsGz=0;


	}

	~Particle(){};

	// accessors
	Vector getMomentum() const { return Vector(itsMomentumX, itsMomentumY, itsMomentumZ); }
	double getMomentumX() const { return itsMomentumX; }
	double getMomentumY() const { return itsMomentumY; }	
	double getMomentumZ() const { return itsMomentumZ; }

	Vector getPosition() const { return Vector(itsPositionX, itsPositionY, itsPositionZ); }
	double getPositionX() const { return itsPositionX; }
	double getPositionY() const { return itsPositionY; }	
	double getPositionZ() const { return itsPositionZ; }

	Vector getPositionOld() const { return Vector(itsPositionXOld, itsPositionYOld, itsPositionZOld); }
	double getPositionXOld() const { return itsPositionXOld; }
	double getPositionYOld() const { return itsPositionYOld; }
	double getPositionZOld() const { return itsPositionZOld; }

	Vector getMomentumOld() const { return Vector(itsMomentumXOld, itsMomentumYOld, itsMomentumZOld); }
	double getMomentumXOld() const { return itsMomentumXOld; }
	double getMomentumYOld() const { return itsMomentumYOld; }
	double getMomentumZOld() const { return itsMomentumZOld; }


	double getRX() const { return itsRX; }
	double getRY() const { return itsRY; }
	double getRZ() const { return itsRZ; }

	double getMass() const{return itsMass;}
	
	void setMomentum(double px, double py, double pz){
	
		
		itsMomentumX=px;
		itsMomentumY=py;
		itsMomentumZ = pz;
	
	
	};

	void setR(double px, double py, double pz){


		itsRX = px;
		itsRY = py;
		itsRZ = pz;


	};
	
	void setPosition(double qx, double qy, double qz){
	
				
		itsPositionX=qx; 
		itsPositionY=qy; 
		itsPositionZ = qz;
		
	

	};


	void setPositionOld(){
		
		itsPositionXOld = itsPositionX;
		itsPositionYOld = itsPositionY;
		itsPositionZOld = itsPositionZ;
		
		itsFxOld = itsFx;
		itsFyOld = itsFy;
		itsFzOld = itsFz;
	};

	void setMomentumOld(){

		itsMomentumXOld = itsMomentumX;
		itsMomentumYOld = itsMomentumY;
		itsMomentumZOld = itsMomentumZ;


	};


	void setFOld(){

		itsFx = itsFxOld;
		itsFy = itsFyOld;
		itsFz = itsFzOld;


	};


	double getFX() const {return itsFx;};	
	double getFY() const {return itsFy;};	
	double getFZ() const {return itsFz;};	
	void setF(double fx, double fy, double fz){itsFx=fx;itsFy=fy;itsFz=fz;}



	double getEpsr() const {return itsepsr;};
	double getEpsf() const {return itsepsf;};

	void setParticleType(unsigned int t){
		type=t;
	
	/*	if (type == 0)
			itsMass = 1.0;
		else
			itsMass = 1.0;*/
		
	};
	unsigned int getParticleType() const{return type;};




	void decideActiveNonActive(){
	
	

		if (fabs(itsMomentumX) / itsMass < itsepsr && fabs(itsMomentumY) / itsMass < itsepsr && fabs(itsMomentumZ) / itsMass< itsepsr)
		{
			howManyRestrained++;
			active = 0;
		}

		else{
			active = 1;
			howManyActive++;
		}




		
	
			
	}
	
	double decideActiveNonActiveOneCoordinate(double px){

		double fpx = fabs(px)/ itsMass;
			
		if (fpx  >= itsepsf)
			{
				return 1;
				
				
			}
		if (fpx  <= itsepsr)
			{
				return 0;
			}

			else
				return 2;	


	}


	void updateActiveStatus(){


		

		if (fabs(itsMomentumX) / itsMass<itsepsr && fabs(itsMomentumY) / itsMass<itsepsr && fabs(itsMomentumZ) / itsMass<itsepsr)

			active = 0;
		
		
		else 
			active = 1;
		

	}



	unsigned int getActiveStatus() const {return active;}
	double getEkin() const{ return (itsMomentumX*itsMomentumX + itsMomentumY*itsMomentumY + itsMomentumZ*itsMomentumZ) / (2.0*itsMass); }
	double getEkinOld() const{ return (itsMomentumXOld*itsMomentumXOld + itsMomentumYOld*itsMomentumYOld + itsMomentumZOld*itsMomentumZOld) / (2.0*itsMass); }

	static int howManyRestrained;
	static int howManyActive;

	void setEps(double er, double ef){ itsepsr=er; itsepsf=ef;};

	//void cleanCounterOfRestrainedParticles(){howManyRestrained =0;};

//	std::map<Particle*, double> forceMap;

	void setParticleIndex(unsigned int i){ index = i; };
	unsigned int getParticleIndex(){return index; };



	void setG(double gx, double gy, double gz){ itsGx = gx; itsGy = gy; itsGz = gz; };

	void setGOld(){
		itsGOldx = itsGx;
		itsGOldy = itsGy;
		itsGOldz = itsGz;
	};
	double GetGx(){ return itsGx; };
	double GetGy(){ return itsGy; };
	double GetGz(){ return itsGz; };

	double GetGOldx(){ return itsGOldx; };
	double GetGOldy(){ return itsGOldy; };
	double GetGOldz(){ return itsGOldz; };

private:

	double itsMomentumX;
	double itsPositionX;

	double itsMomentumY;
	double itsPositionY;

	double itsMomentumZ;
	double itsPositionZ;


	double itsPositionXOld;
	double itsPositionYOld;
	double itsPositionZOld;
	
	double itsMomentumXOld;
	double itsMomentumYOld;
	double itsMomentumZOld;

	double itsRX;
	double itsRY;
	double itsRZ;


	double itsMass;
	double itsFx;
	double itsFy;
	double itsFz;

	double itsFxOld;
	double itsFyOld;
	double itsFzOld;
	//double itsFcomingFromInteractions;

	unsigned int active; // 1 activ, 0 non activ (restrained), 2 transition phase
	
	double itsepsr;
	double itsepsf;

	unsigned int type; //0 solvent, 1 dimer

	unsigned int index;

	double itsGOldx;
	double itsGx;
	double itsGOldy;
	double itsGy;
	double itsGOldz;
	double itsGz;

};

