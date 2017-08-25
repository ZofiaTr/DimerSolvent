#include "particle.h"

#include "SimpleAverage.h"
#include "SBCTimeClock.hpp"
#include "SBDTypeRandom.hpp"
#include <iostream>

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>

class NeighborSearchGridPBC;

//************** DIFFERENT FORCES UPDATE ALGORITHMS
// newton -> use f_ij=-f_ji -> i=0:N j=i+1:N
// without newton-> i=0:N and j=0:N


//#define SPLINE_FIRST_ORDER
//#define SPLINE_SECOND_ORDER
//#define SPLINE_THIRD_ORDER
//#define SPLINE_FOURTH_ORDER
//#define EXP_PERTUBATION
//#define QUARTIC_POLYNOMIAL
#define POLYNOMIAL_POWERS

//#define GAMMA_ANALYSIS //goes with metropolis !!!
//#define WRITE_INSTANT_REJECTION_RATE
#define RATE_HAM 0.002
#define RATE_FD  0.002


// choose standard (without arps overhead in the kinetic energy functions) or adaptive (allows to take general kinetic energy function)
//#define STANDARD_UPDATE_FORCES //usual algo for standard case- no active status
#define ADAPTIV_UPDATE_FORCES

#define CONDITION_ADAPTIVE_FORCES_UPDATE 0.7071067812 // condition from complexity analysis for the forces update - Redon Trstanova 2016

//#define CHANGE_PARAMETERS_ON_THE_FLY //for time measurements : allows to change aprs parameters and to restart the trajectory
#define MEASURE_TIME_PER_TIMESTEP
#define MEASURE_TIME_PER_FORCE_UPDATE //measures only the ADD function (see arps algo), time is then multiplied by 2 to express also SUBTRACT

//#define OPTIMAL_PARAMETERS // in main.cpp: 3 points, goes with VARIANCE_BY_AUTOCORRELATION_FUNCTION
//#define RUN_REFERENCE_SIMULATION_WITH_STD

//#define ERROR_ALERT_PBC //check blow up by checking whether the position is correct after applying of PBC

//#define NO_FORCES

//**************************AVERAGES AND TIME STEP ANALYSIS ***********

#define SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS

#define NUMBER_OF_AVERAGES 11 //for file AVERAGES.m - how many values written
#define HOW_MANY_PHYSICAL_TIMES 1000 //for time step analysis: divide physical time interval into HOW_MANY_PHYSICAL_TIMES subintervals to be written into file AVERAGES.m

#define ADD_SAMPLES

//*******************************************
//#define NVE

//#define ONE_PARTICLE_ONLY //in update forces: it allows to fix the particle to a "wall" (otherwise noone to interact with)

//******* WRITING INTO A FILE *****************************************

//#define WRITE_INST_ENERGY
//#define WRITE_POSITIONS_MATLAB
//#define WRITE_MOMENTA_MATLAB
//#define WRITE_MOMENTA_MATLAB_HISTOGRAM
//#define WRITE_DIMER_DISTANCE_MATLAB
//#define WRITE_ERROR_LANGEVIN
//#define RDF
//#define DISTRIBUTION_RESTR_PART
//#define WRITE_BA_RESULT // write the block averaging into matlab file in order to verify the "plateau" that estimated the stdev
//#define WRITE_P_SQUARE_2D

#define WRITE_AVER_TIME_PER_TIME_STEP //
//************************************************
//#define BLOCK_AVERAGING

//in which domain of the kinetic energy function i am at most - simpleaverages full/restr/spline
//#define KINETIC_ENERGY_DOMAIN_COUNTERS



//************** VARIANCE_BY_AUTOCORRELATION_FUNCTION**************************
// measuring variance: uncomment VARIANCE_BY_AUTOCORRELATION_FUNCTION

#define VARIANCE_BY_AUTOCORRELATION_FUNCTION

//#define ACF_SIMPSON
#define ACF_TRAPEZ

//do not comment out! (constants)
// variance 1 -> total potential and variance 2 dimerpotential
//PHYSICAL TIME !!!
#define TIME_AUTOCORRELATION_1 2//0 // solvent solvent potential
#define TIME_AUTOCORRELATION_2 15//0 // dimer position//potential
#define TIME_AUTOCORRELATION_3 5//2 // temperature//solvent dimer
//**********************************************************

#define VARIANCE_IR
#define LENGTH_N_IN_REPLICAS 100 // physical time !!!

#define MORE_VARIANCE_IR // with different length in replicas
#define LENGTH_N_IN_REPLICAS2 10 // physical time !!!
//**********************************************************

//#define BACKUP_VARIANCE_EVOLUTION // write variance over time step

#define MAX_LENGTH_OF_VARIANCE_ARRAY 100000
#define MAX_LENGTH_OF_ACF_ARRAY 500000
#define WRITING_AUTOCORRELATION_TO_FILE 100000//1000000

//*****************************************************************************
// plot histograms of max nabla U -> save also correspondant p

//#define HIST_MAX_NABLA_U // to be rewritten for 3d!!!
#define pMAX 6
#define NABLA_U_MAX 6

#define pSTEP 0.05

//#define HISTOGRAM // to be rewritten for 3d!!!
#define MAX_VALUE_HISTOGRAM 2

//*****************************************************************************
// DISCRETIZATION SCHEMES

#define CHECK_SPLITTING_SCHEME //write on screen the names of functions called from PerformTimeStep in order to check the integration scheme (first/second order (Midpoint, CN, strang))

#define METROPOLIZATION //can stand (uncomment) alone
#define METROPOLIS_FDR
//#define METROPOLIS_FDR_ADAPTIVE

//#define METROPOLIS_LANGEVIN

//#define FIRST_ORDER_SPLITTING //can stand (uncomment) alone, oposite of second order spliiting which need either second_order_fd_cn or ..midpont
//#define SECOND_ORDER_SPLITTING
//#define SECOND_ORDER_FD_CN
//#define SECOND_ORDER_STRANG //MIDPOINT VERLET MIDPOINT

//#define SECOND_ORDER_FD_MIDPOINT
//#define ONLY_FLUCTUATION_DISSIPATION_PART //optional for second order splitting
//#define OTHER_KINETIC_ENERGY //take U(p)=p^4/4 instead of the standard kinetic energy or the arps kinetic energy


#define BAOAB_SCHEMES
//#define Z_BAOAB
//#define BAOAB
//#define ABOBA
//#define OBABO
//#define BABM
//#define MBABM
//#define BABE
//#define BAB_MOMENTA_LANGEVIN
//#define FD_EXACT
//#define FD_ONLY
//#define OBABO_noNBL

//-----------------------------------

//#define READ_INITIAL_CONDITION_FROM_FILE
#define SAME_SEED 0
#define SEED_VALUE 0

#define PI 3.14159265359
#define MAX_ARRAY_LENGTH 10000000
#define BA_STYLE 2 // 1: take last MAX_ARRAY_LENGTH time steps
					// 2: take every WritingPeriodBlockaveraging'th time step and then unscale

#define SOLVENT_LJ
//#define SOLVENT_COS //we should adjust the force to be zero at the cut_off..

//#define COMPUTE_MOMENTA_SQUARE_AVERAGE


class Simulation {

public:

	Simulation(int NrParticles);
	~Simulation();

	void	ActiveStatus();
	void	AddSamples();		// add observables into average objects
	void	AddSamplesChangingParameters(); // add averages for time mesumerent for restrainedPartPercentageTime and numberOfInteractionsPerRestrPart
	double*	BlockAveraging();														// block averaging for variance estimation
	void	CleanSamples();															// clean averages objects
												// convert seconds into hours/minutes/seconds and write it to screen

	double*	ComputeBlockAveraging(double *posBA, int NM);							// block averaging function
	double	ComputeEnergy();
	double	ComputeKineticEnergy();
	void	ComputeRadialDistributionFunction();
	void	ComputeDistributionOfRestrainedParticles();
	float	ComputeBoxMuller(float m, float s);										// normal random variate generator, mean m, standard deviation s
	double ComputeDimerDistance();

	double	ComputePotentialEnergy();
	double	ComputePressure();
	double	ComputeTemperature();
	double ComputeConfigurationalTemperature();
	double	ComputeTimeLeftToTheEndOfEXE(double s);
	void	ConvertTimeUnits(double s);


	double	EstimateVarianceFromBA(double *v, int lengthv);							// block averaging function - search for plateau - take first local max as variance value
	void	DisplayTimeLeft();														// estimate time left from current time step and write it to the screen
	void	DisplayTimeLeft(int n, int writingperiod, int N, int measure_t);


	// kinetic energy function

	double HkineticEnergyFunction(double px, double py, double pz, double mass, double epsrp, double epsfp);

	double	Hp(double px, double py, double pz, double mass, double epsrp, double epsfp);					// derivative of the Hamiltonian H_arps in p direction
	double	Hpx(double px, double mass, double epsrp, double epsfp);					// derivative of the Hamiltonian H_arps in p direction
	double	dHdp(double px, double py, double pz, double mass, double epsrp, double epsfp);					// derivative of the Hamiltonian H_arps in p direction
	double	dHdp3(double px, double py, double pz, double mass, double epsrp, double epsfp);					// derivative of the Hamiltonian H_arps in p direction
	double	dHdp2xx(double px, double py, double pz, double mass, double epsrp, double epsfp);					//second derivative of the Hamiltonian H_arps in p direction
	double	dHdp2xy(double px, double py, double pz, double mass, double epsrp, double epsfp);

	double	dHdpNablaLaplacexx(double px, double py, double pz, double mass, double epsrp, double epsfp);					// third derivative of the Hamiltonian H_arps in p direction
	double	dHdpNablaLaplacexy(double px, double py, double pz, double mass, double epsrp, double epsfp);					// third derivative of the Hamiltonian H_arps in p direction

	// ARPS FUNCTIONS
	double	Rho_p(double px, double py, double pz, double epsrp, double epsfp, double mass);				// arps function rho (zeta)
	double	dRho_p(double px, double py, double pz, double epsrp, double epsfp, double mass);				// derivative of arps function rho (zeta) in p direction
	double	dRho_p2xx(double px, double py, double pz, double epsrp, double epsfp, double mass);
	double	dRho_p2xy(double px, double py, double pz, double epsrp, double epsfp, double mass);
	double	dRho_p3xx(double px, double py, double pz, double epsrp, double epsfp, double mass);
	double	dRho_p3xy(double px, double py, double pz, double epsrp, double epsfp, double mass);
	double	dRho_p3(double px, double py, double pz,double epsrp, double epsfp, double mass);

	// nabla U (kinetic energy)
	double dUdpx(double px, double py, double pz, double epsri, double epsfi, double mass);
	double dUdpy(double px, double py, double pz, double epsri, double epsfi, double mass);
	double dUdpz(double px, double py, double pz, double epsri, double epsfi, double mass);


	double computeDeterminant2(double** A);
	double computeDeterminant3(double** A);

	double	getCurrentTimeInSeconds() const;										// return currentTimeInSeconds
	double* getArrayTimePerTimeStep() const;
	double GetPercentageOfRestrained(){ return RestrainedParticles.getAverage(); };

	double*	getSamples();	// return vector of samples
	double	getVar_mu() const; // first entry from BA to compute correlation length N_corr
	double	getStandardDeviation();											// return standard deviation
	double	getValues(int i) const;
	double	getObservableForVariance(int varianceOneOrVarianceTwo);

//	double	getVariance() ;// return internal constants - check function definition

	void	InitialCondition();
	double	MidPoint(double x,double y);





//	 metropolization
//#ifdef METROPOLIZATION

	double MetropolisHamiltonianARPS(double E_old);
	double TproposalARPS(double pOldx, double pOldy, double pOldz, double pNewx, double pNewy, double pNewz, double m, double epsrn, double epsfn, double t);
	void MetropolisStochARPS(double Told, double t);
	void MetropolisEulerMaruyama(double t);
	double TproposalEulerMaruyama(double pxNew, double pyNew, double pzNew, double pxOld, double pyOld, double pzOld, double m, double epsrn, double epsfn, double t);
	void MetropolisStochLangevin(double t);
	void MetropolisStochLangevinSeparableSpace(double t);

	void MetropolisFDR(double* oldKinEn, double h);
	void SampleG();
	double* ComputeEnergyFDR();


//#endif
	double Observable1();
	double Observable2();
	double Observable3();

	void TproposalAlternative(double Told, double t);

		// integration of the dynamics
	void	PerformTimeStep();

//	void	PerformTimeStepMetropolisARPS();


	//void	PerformLearningTrajectory();

	void	PerformTrajectory();
	void	PerformTrajectoryChangingParameters(); //used for time estimation per time step !!!

//	void	PerformTrajectoryIndependentReplicas();

	//void	PerformTrajectoryNVE(double er, double ef);
//	void	PerformTrajectoryMetropolis();
//	void	PerformTrajectorySecondOrder();

	void	PeriodicBoundaryCondition();
	int		roundMy(double d1);
	int		roundUp(double x);
	double** readMatrixFromFile(int nrRows, int nrColumns);
	void	setCurrentTimeInSeconds(double sec);									// set currentTimeInSeconds
	void	setNumberStepsFORloop(double physTime, double dtMinimal, int NumberOfDifferentDt, double deltaDT);
	void	setSimulationParameters(double NrOfTimeSteps, double TimeStepSize,  double er, double ef, int nrDimer, int nrRepl, double initialDensity, double nrEquil); // set parameters: nrOftimesteps, dt, epsilon_r, epsilon_f

	void	SampleVarianceIR();
	void	SampleVarianceIR2();
	void	Tick();																	// time measuring from inside - change t0 variable
	void	Tock();																	// time measuting from inside - change measureTime variable- estimate average time per time step
	//double	TproposalARPS(double pOld, double pNew, double m, double epsrn, double epsfn);		// metropolis-hastings function

	// discretization
	void	UpdateMomentaParticles(double t);
	void	UpdateMomentaStochOrder1( double t);
	void	UpdateMomentaStochOrder2( double t);
	void	UpdateMomentaStochOrder2_CN(double t);
	void	UpdateMomentaStochExact(double t);
	void	UpdateMomentaStochModifiedEqCN(double t);
	void	UpdateMomentaStochSecondOrder(double t);
	void	UpdateMomentaStochLangevin(double t);

	double AnalyticalSolutionFD(double t, double px, double mass, double epsrp, double epsfp);

	void	UpdateMomentaStochExactARPS(double t);

	void	UpdateMomentaStochExtendedMidpoint(double t);
	void	ComputeAdditionalTermU3(double t);
	void	UpdateMomentaARPSperturbation(double t);
	void	UpdateMomentaStochModifiedEq(double t);
	void	UpdateMomentaVerlet_FDR(double h);
	void	UpdateMomentaVerlet_FDR_MetropolizePerParticle(double h);
	void	UpdateMomentaFD_BDV(double t); // Bou-Rabee Donev Vanden-Eijnden 2013


	void	UpdatePositionsParticles(double t);
	void	UpdatePositionsParticlesStandard(double t);
	void	UpdateForcesParticles();
	void	UpdateForcesParticlesAdaptively();
	void	UpdateForcesParticlesAdaptivelyARPS_ADD();
	void	UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
	void	UpdateForcesParticlesAdaptivelyNewtonARPS_ADD(); //newton between active neighbors
	void	UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT(); //newton between active neighbors
	void	UpdateForcesParticlesNewton();
	void	UpdateForcesParticlesAdaptivelyNewton();

	/*void	ForcesCalculationActiveParticles(int i, int j);
	void	ForcesCalculationActiveParticles(Particle* particleI, Particle* particleJ);
	void	ForcesCalculationActiveParticlesNewton(int i, int j);
	void	ForcesCalculationRestrainedParticles(int i, int j);
	void	ForcesCalculationRestrainedParticlesNewton(int i, int j);
	void	ForcesCalculation(int i, int j);
	void	ForcesCalculation(Particle* particleI, Particle* particleJ);*/
	double	InteractionForce(double r, unsigned int typeParticle1, unsigned int typeParticle2);
	double	InteractionPotential(double r, unsigned int typeParticle_i, unsigned int typeParticle_j);

	double WCA(double x);
		double fWCA(double x)  ;
		double DW(double x);
		double fDW(double x)  ;

		void SaveCurrentPositionsAsOld();
		void SaveCurrentMomentaAsOld();

	void	WritePositionForBlockAveraging();										// part of AddSamples()
	void	WriteOnScreen();
	void	WritePositionsInFile();
	void	SaveCurrentSimulation();

#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS
	void	SaveCurrentAverages();
	void SaveCurrentAveragesOnlyInstantOnes();
#endif

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION
	void	WriteAutocorrelationAveragedInReplicas();
#endif

protected:

	int NumberOfParticles;
	int NumberOfDimerParticles;

	double NumberOfTimeSteps;
	double NumberOfEquilibrationSteps;
	int NumberOfReplicas;

	double dt;

	unsigned int PBC;

	NeighborSearchGridPBC* grid;


	unsigned long current_n;
	int numberOfCurrent_n;
	double MaxLengthOfCurrentN;
	double currentTimeInSeconds;


	double percentageRestrainedParticles;
	double countPercRestr;

	double epsCoupling;

	double sigma_LJ ,epsilon_LJ,cut_off,d_s,wDimer,	hDimer;

	double* epsr;
	double* epsf;

	int i_rep;

	double omega;
	double fix_point_tolerance;

	double T;
	double k;

	double kT;
	double beta;
	//double alpha;

	double d;
//	double m;


	double gamma;
	double sigma;


	double BoxLength;
	double L;

	double intialEnergyNVE;

	SimpleAverage<double> Temperature;
	SimpleAverage<double> ConfigTemperature;
	SimpleAverage<double> Pressure;
	SimpleAverage<double> Position;
	SimpleAverage<double> Potential;
	SimpleAverage<double> RestrainedParticles;
	SimpleAverage<double> MeasureTime;
	SimpleAverage<double> ForceAppliedOnDimer;
	SimpleAverage<double>	MomentaSquare;

	SimpleAverage<double> VarianceOne;
	SimpleAverage<double> VarianceTwo;
	SimpleAverage<double> VarianceThree;


	SimpleAverage<double> ReplicasAverage;
	SimpleAverage<double> ReplicasAverageThroughPathN;
	SimpleAverage<double> VarianceIR;


	SimpleAverage<double> ReplicasAverage2;
	SimpleAverage<double> ReplicasAverageThroughPathN2;
	SimpleAverage<double> VarianceIR2;

	SimpleAverage<double> AutocorrelationAverage;
	// this could be variance of dimer force or dimer center... -> SimpleAverage<double> VarianceThree;
	SimpleAverage<double> NumberOfInteractions;
	SimpleAverage<double> AverageNumberOfNeighbors;

	double VarianceIndependentReplicas;
	double VarianceIndependentReplicas2;

#ifdef KINETIC_ENERGY_DOMAIN_COUNTERS
	//counters for domain identification - in which part of the kinetic energy I am the most of the time
	SimpleAverage<double> FullDynamicsDomain;
	SimpleAverage<double> RestrainedDynamicsDomain;
	SimpleAverage<double> SplineDomain;
#endif


	//SimpleAverage<double> MeanSquareDisplacement;

	SimpleAverage<double> EnergyErrorNVE;

#ifdef HIST_MAX_NABLA_U

	SimpleAverage<double>* HistNablaU; // average ||\nabla U||Linfty over corresponding pmax

	double* NablaUArray; // ||\nabla U||Linfty over corresponding pmax
	SimpleAverage<double>* NablaUCounterArray; // counting for histogram of NablaUArray

	int numberOfBinsHistpMax;
	int numberOfBinsHistNablaU;


	double** ProbaMatrixNablaU; // matrix for histogram pmax times || NablaU ||

	double probaMatrixNablaUCounterTotal;

#endif

#ifdef HISTOGRAM

	double* Histogram;
	int numerOfBinsInHistogram;
	double histogramCounter = 0;
	double stepHist;

#endif


	Particle** pParticles;


	int WritingPeriodOnScreen;
	int WritingPeriodFile;
	int WritingPeriodBlockaveraging;


	double  t0;
	double  t1;

	double measureTime;

	// if doing for loop with different time step size
	double physicalTime;
	double dtMinimalLoop; //dtMin
	int nrOfDTloop;
	double deltaTimeStepSize;

	double SumOfnrOfStepsInTimeStepSizeLoop; // number of steps already passed

	SBRandom* randomGenerator;
	SBRandom* rMetropolis;

	double* positionBA;
//	double* positionAverage;
	double* positionBA_scaled;
	double* positionBA_scaled2;

	double ValueObservable_IC;

	int BA_count;

	double transitionTime ;
	SimpleAverage<double> TransitionTime;
	int numberOfStepsFromPreviousState;
	int dimerState;

	double standardDeviation;
	double standardDeviation2;

	double potentialEnergy;
	double potentialEnergyOld;

	double kineticEnergy;
	double temperatureMomenta;
	double temperatureForSamples;

	double DimerDistance;
	double DimerPotential;
	double DimerSolventPotential;
	double SolventSolventPotential;

	double oldEnergy;
	bool errorIndicator;

	double interactionsCounter;
	double Var_mu;



	// adaptive


	Particle** activeParticleArray;
	int* indexArrayActiveParticles;
	int* indexArrayRestrainedParticles;

	SimpleAverage<double>* restrainedPartPercentageTime;
	SimpleAverage<double>* restrainedPartPercentageTimePerTimeStep;
	SimpleAverage<double>* numberOfInteractionsPerRestrPart;
	SimpleAverage<double>  timePerTimeStep;
	SimpleAverage<double>  timePerForceUpdate;


	//SimpleAverage<double>	 AveragePercentage;

	SBTime tPerTimeStep0;
	SBTime tPerTimeStep1;
	//other two for add/substract functions in arps..
	SBTime tPerTimeStep3;
	SBTime tPerTimeStep4;


	SimpleAverage<double> RejectionRate;
	SimpleAverage<double> RejectionRateFD;

#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS
	double** AverageMatrix;
	int AverageCounter; //how many times an average was written into file upto NUMBER_OF_AVERAGES
#endif

#ifdef RDF
	//--------------Radial distribution function---------

	double* gRDF;
	double numberOfBins;
	double drRDF;

#endif


#ifdef DISTRIBUTION_RESTR_PART
	//--------------Distribution of restrained particles relatively to center of dimer ---------

	/*double* gDRP;
	double numberOfBins_DRP;
	double drDRP;
*/


	SimpleAverage<double>* gDRP;
	SimpleAverage<double>* gDRPRestrained;

	int numberOfBins_DRP;
	double drDRP;

#endif

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

	SimpleAverage<double>* AutocorrelationArrayReplica1;
	SimpleAverage<double>* AutocorrelationArrayReplica2;
	SimpleAverage<double>* AutocorrelationArrayReplica3;

	SimpleAverage<double>* AveragesArrayReplica1;
	SimpleAverage<double>* AveragesArrayReplica2;
	SimpleAverage<double>* AveragesArrayReplica3;

	double NumberOfStepsInTimeInterval1;
	double NumberOfStepsInTimeInterval2;
	double NumberOfStepsInTimeInterval3;

	double InitialValueAutocorrelation1 ;
	double InitialValueAutocorrelation2 ;
	double InitialValueAutocorrelation3 ;

	double countReplicasAutocorrelation1;
	double countReplicasAutocorrelation2;
	double countReplicasAutocorrelation3;

	double* ConvergenceEstimatedVarianceInNumberOfReplicas1;
	double* ConvergenceEstimatedVarianceInNumberOfReplicas2;
	double* ConvergenceEstimatedVarianceInNumberOfReplicas3;


	//Variance computed by ACF
	double IntegralAutocorrelation1;
	double IntegralAutocorrelation2;
	double IntegralAutocorrelation3;

#endif

#ifdef VARIANCE_IR
	double LengthOfPartialTrajectoryN;
#endif

#ifdef MORE_VARIANCE_IR
	double LengthOfPartialTrajectoryN2;
#endif


#ifdef GAMMA_ANALYSIS


	double fixedRejectionRateHam;
	double fixedRejectionRateFD;
	double toleranceRejectionRate;



#ifdef WRITE_INSTANT_REJECTION_RATE

	double instantRejectionHam;
	double instantRejectionFD;

#endif

#endif


	unsigned int doAddSubst ;

};

#define S Simulation
