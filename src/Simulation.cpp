#include "Simulation.h"
#include "NeighborSearchGridPBC.hpp"


Simulation::Simulation(int NrParticles){

	dt = 0.01;

	NumberOfParticles = NrParticles;
	NumberOfDimerParticles = 0;
	NumberOfTimeSteps = 0.0;
	NumberOfEquilibrationSteps = 0;//


	WritingPeriodOnScreen = 10000;//for simulations on cluster use 100000000;
	WritingPeriodFile = 100000;
	WritingPeriodBlockaveraging = 10;//10;

	BA_count = 0;
	current_n = 0; // time step n in perform trajectory
	numberOfCurrent_n = 0;
	currentTimeInSeconds = 0.0; // time measured from outside of the class
	MaxLengthOfCurrentN = pow(10.0, 9);


	PBC = 1;

	intialEnergyNVE = 0.0;

	epsCoupling = 0.0;

	percentageRestrainedParticles = 0.0;
	countPercRestr = 0;

	//arps thesholds
	epsr = new double[NumberOfParticles];
	epsf = new double[NumberOfParticles];

	for (int n = 0; n < NumberOfParticles; n++){

		epsr[n] = 0.0;
		epsf[n] = 0.0;

	}


	omega = 0.1;
	fix_point_tolerance = pow(10, -10);

	sigma_LJ = 1.0;
	epsilon_LJ = 1.0;

	cut_off = pow(2., 1. / 6)*sigma_LJ;
	//cut_off=2.5*sigma_LJ;
	d_s = (cut_off - 0.1);

	wDimer = 1;
	hDimer = 1;

	T = 1.0;
	k = 1.0;
	std::cout << "set Temperature is " << T << std::endl;

	kT = k*T;
	beta = 1 / kT;


	//alpha = 1;//0.25;

	d = 3.0;
	//m = 1.0;


	gamma = 1;
	sigma = sqrt(2.0*kT*gamma);

	BoxLength = cut_off * 4;//2*PI;
	L = BoxLength;//sqrt(NumberOfParticles*k*T/alpha);


	measureTime = 0.0;
	t0 = 0.0;
	t1 = 0.0;


	physicalTime = 0.0;
	dtMinimalLoop = 0.0;
	nrOfDTloop = 0;
	deltaTimeStepSize = 0.0;

	SumOfnrOfStepsInTimeStepSizeLoop = 0.0;

	i_rep = 1;


#if SAME_SEED == 1
	unsigned int seed = SEED_VALUE;
#endif
#if SAME_SEED == 0
	unsigned int seed = (unsigned int)SBCTimeClock();
#endif

	randomGenerator = new SBRandom(seed);

	//unsigned int seed = 0;// (unsigned int)SBCTimeClock();
	//randomGenerator= new SBRandom(seed);
	rMetropolis = new SBRandom(SEED_VALUE + 1);

	ValueObservable_IC = 0; // value of observable A(q_0,p_0) for autocorrelation

	positionBA = new double[MAX_ARRAY_LENGTH];
	//positionAverage=new double[MAX_ARRAY_LENGTH];
	positionBA_scaled = new double[MAX_ARRAY_LENGTH];
	positionBA_scaled2 = new double[MAX_ARRAY_LENGTH];

	for (int n = 0; n < MAX_ARRAY_LENGTH; n++){

		positionBA_scaled[n] = 0.0;
		positionBA_scaled2[n] = 0.0;
		positionBA[n] = 0.0;
		//	positionAverage[n]=0.0;


	}

	potentialEnergy = 0.0;
	potentialEnergyOld = 0.0;

	temperatureMomenta = 0.0;
	kineticEnergy = 0.0;
	temperatureForSamples = 0.0;

	DimerDistance = 0.0;
	DimerPotential = 0.0;

	standardDeviation = 0.0;
	standardDeviation2 = 0.0;

	VarianceIndependentReplicas = 0.0;
	VarianceIndependentReplicas2 = 0.0;

	transitionTime = 0.0;
	numberOfStepsFromPreviousState=0;
	dimerState=1;

	/*
		KineticEnergyParticle= new SimpleAverage<double>[NumberOfParticles];
		VarianceReplicas = new SimpleAverage<double>[NumberOfParticles];
		ObservableReplicas = new SimpleAverage<double>[NumberOfParticles];
		*/

	errorIndicator = false;

	interactionsCounter = 0.0;
	Var_mu = 0.0;

	//adaptive

	tPerTimeStep0 = SBCTimeClock();
	tPerTimeStep1 = SBCTimeClock();
	tPerTimeStep3 = SBCTimeClock();
	tPerTimeStep4 = SBCTimeClock();

	restrainedPartPercentageTime = new SimpleAverage<double>[100];
	restrainedPartPercentageTimePerTimeStep = new SimpleAverage<double>[100];
	numberOfInteractionsPerRestrPart = new SimpleAverage<double>[100];
#ifdef RDF
	//--------------Radial distribution function---------
	numberOfBins = 500;

	gRDF = new double[(int)numberOfBins];

	for (int i = 0; i<numberOfBins; i++)
	{
		gRDF[i] = 0;
	}

	drRDF = 0;

#endif



#ifdef DISTRIBUTION_RESTR_PART
	//--------------Radial distribution function---------
	numberOfBins_DRP = 50;//500;

	gDRP = new SimpleAverage<double>[(int)numberOfBins_DRP];
	gDRPRestrained = new SimpleAverage<double>[(int)numberOfBins_DRP];
	/*
		for (int i = 0; i<numberOfBins_DRP; i++)
		{
		gDRP[i] = 0;
		}
		*/
	drDRP = 0;

#endif


#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION


	AutocorrelationArrayReplica1 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];
	AutocorrelationArrayReplica2 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];
	AutocorrelationArrayReplica3 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];

	AveragesArrayReplica1 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];
	AveragesArrayReplica2 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];
	AveragesArrayReplica3 = new SimpleAverage<double>[MAX_LENGTH_OF_ACF_ARRAY];

	ConvergenceEstimatedVarianceInNumberOfReplicas1 = new double[MAX_LENGTH_OF_VARIANCE_ARRAY];
	ConvergenceEstimatedVarianceInNumberOfReplicas2 = new double[MAX_LENGTH_OF_VARIANCE_ARRAY];
	ConvergenceEstimatedVarianceInNumberOfReplicas3 = new double[MAX_LENGTH_OF_VARIANCE_ARRAY];

	for (int i = 0; i < MAX_LENGTH_OF_VARIANCE_ARRAY; i++) ConvergenceEstimatedVarianceInNumberOfReplicas1[i] = 0.0;
	for (int i = 0; i < MAX_LENGTH_OF_VARIANCE_ARRAY; i++) ConvergenceEstimatedVarianceInNumberOfReplicas2[i] = 0.0;
	for (int i = 0; i < MAX_LENGTH_OF_VARIANCE_ARRAY; i++) ConvergenceEstimatedVarianceInNumberOfReplicas3[i] = 0.0;

	IntegralAutocorrelation1 = 0.0;
	IntegralAutocorrelation2 = 0.0;
	IntegralAutocorrelation3 = 0.0;


	InitialValueAutocorrelation1 = 0.0;
	InitialValueAutocorrelation2 = 0.0;
	InitialValueAutocorrelation3 = 0.0;

	countReplicasAutocorrelation1 = 0.0;
	countReplicasAutocorrelation2 = 0.0;
	countReplicasAutocorrelation3 = 0.0;

#endif


#ifdef VARIANCE_IR

	LengthOfPartialTrajectoryN = 0;// (int)(LENGTH_N_IN_REPLICAS / dt);
#endif

#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS

	AverageMatrix = new double*[NUMBER_OF_AVERAGES];
	for (int i = 0; i < NUMBER_OF_AVERAGES; i++){
		AverageMatrix[i] = new double[MAX_ARRAY_LENGTH];

		for (int j = 0; j < MAX_ARRAY_LENGTH; j++)AverageMatrix[i][j] = 0.0;
	};

	AverageCounter = 0;
#endif

#ifdef HIST_MAX_NABLA_U

	double pmax = pMAX;
	double maximalNablaU = NABLA_U_MAX;
	double pstep = pSTEP;
	numberOfBinsHistpMax = (int)(2.0 * pmax / pstep); //  discr -pmax.......0........pmax
	numberOfBinsHistNablaU = (int)(maximalNablaU / pstep);

	HistNablaU = new SimpleAverage<double>[numberOfBinsHistpMax];


	NablaUArray= new double[numberOfBinsHistpMax];
	for (int i = 0; i < numberOfBinsHistpMax; i++) NablaUArray[i] = 0.0;

	NablaUCounterArray=new SimpleAverage<double>[numberOfBinsHistpMax];




	probaMatrixNablaUCounterTotal = 0.0;


	ProbaMatrixNablaU = new double*[numberOfBinsHistNablaU];
	for (int i = 0; i < numberOfBinsHistNablaU; i++){

		ProbaMatrixNablaU[i] = new double[numberOfBinsHistpMax];
		for (int j = 0; j < numberOfBinsHistpMax; j++){

			ProbaMatrixNablaU[i][j] = 0.0;
		}

	}


#endif

#ifdef HISTOGRAM

	double maxValueHistogram = MAX_VALUE_HISTOGRAM;
	stepHist = pSTEP;
	numerOfBinsInHistogram = (int)(maxValueHistogram / stepHist);


	Histogram = new double[numerOfBinsInHistogram];
	for (int i = 0; i < numerOfBinsInHistogram; i++) Histogram[i] = 0.0;
	histogramCounter = 0;

#endif

	doAddSubst = 0;



#ifdef GAMMA_ANALYSIS

	fixedRejectionRateHam = RATE_HAM;

	//fixedRejectionRateHam = fixedRejectionRateHam / 100;

	fixedRejectionRateFD = RATE_FD;
	//fixedRejectionRateFD = fixedRejectionRateFD /100;

	toleranceRejectionRate = 0.1;



#ifdef WRITE_INSTANT_REJECTION_RATE

	 instantRejectionHam=0.0;
	 instantRejectionFD=0.0;

#endif


#endif

};

Simulation::~Simulation(){

	delete epsr;
	delete epsf;

	delete[] positionBA;
	delete[] positionBA_scaled;
	delete[] positionBA_scaled2;

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

	delete[] ConvergenceEstimatedVarianceInNumberOfReplicas1;
	delete[] ConvergenceEstimatedVarianceInNumberOfReplicas2;
	delete[] ConvergenceEstimatedVarianceInNumberOfReplicas3;
#endif
	//for (int i = 0; i<100; i++) 	restrainedPartPercentageTime[i].SimpleAverage<double>::~SimpleAverage();

#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS
	for (int i = 0; i < NUMBER_OF_AVERAGES; i++) delete[] AverageMatrix[i];
	delete[] AverageMatrix;
#endif

	std::cout << "i am in simulation destructor" << std::endl;

};


void Simulation::ActiveStatus(){

	int count = 0;

	Particle::howManyRestrained = 0;
	Particle::howManyActive = 0;

	int	numberOfActiveParticles = 0;
	int numberOfRestrainedParticles = 0;

	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->decideActiveNonActive();

		if (pParticles[i]->getActiveStatus() == 1) {

			//activeParticleArray[numberOfActiveParticles] = pParticles[i];
			//indexArrayActiveParticles[numberOfActiveParticles] = i;
			numberOfActiveParticles++;


		}
		else {
			//indexArrayRestrainedParticles[numberOfRestrainedParticles] = i;
			numberOfRestrainedParticles++;
		}

	}

	Particle::howManyActive = numberOfActiveParticles;
	Particle::howManyRestrained = numberOfRestrainedParticles;

	/*if ((Particle::howManyActive + Particle::howManyRestrained) != NumberOfParticles) {
		std::cout << "NrActiv+NrRestr != NrPart\n" << std::endl;

		while (1);
	}*/

	for (int i = 0; i < NumberOfDimerParticles; i++){

		if (pParticles[i]->getActiveStatus() != 1){

			std::cout << "Error in ActiveStatus()- dimer particle defined as not active!\n";
			while (1);
		}
	}


}

double Simulation::Observable1(){


	return ComputePotentialEnergy();// DimerPotential;//potentialEnergy;//SolventSolventPotential+DimerPotential+DimerSolventPotential;

};



double Simulation::Observable2(){



	if (NumberOfDimerParticles > 2) std::cout << "Observable 2 dimer distance only for 2 particles!" << std::endl;

	return ComputeDimerDistance();//DimerPotential;

};


double Simulation::Observable3(){


	return ComputeTemperature();//temperatureForSamples;//DimerPotential;

};

void Simulation::AddSamples(){


	temperatureForSamples =  Observable3();
	//ComputePotentialEnergy();

	potentialEnergy = Observable1();
	DimerDistance = Observable2();

	//ConfigTemperature.addSample(ComputeConfigurationalTemperature());

	//AutocorrelationAverage.addSample(ValueObservable_IC*(DimerDistance - VarianceOne.getAverage()));

	Temperature.addSample(temperatureForSamples);
	//Pressure.addSample(ComputePressure());

	Position.addSample(DimerDistance);

	Potential.addSample(potentialEnergy);//(ComputePotentialEnergy());

	RestrainedParticles.addSample(100 * (double)Particle::howManyRestrained / NumberOfParticles);// 100 * percentageRestrainedParticles);

	VarianceOne.addSample(Observable1());
	VarianceTwo.addSample(Observable2());
	VarianceThree.addSample(Observable3());

	// transition time between the two states of dimer
	if(NumberOfDimerParticles == 2){
		int stateFlag =0;

		double tol = 0.01;
			if (fabs(DimerDistance - cut_off)/fabs(cut_off) < tol ){

				stateFlag = 1;


			}
			else if (fabs(DimerDistance - (cut_off + 2.0*wDimer))/ fabs(cut_off) < tol){

				stateFlag = 2;
			}

			//std::cout<<DimerDistance<<std::endl;

		numberOfStepsFromPreviousState++;
		if(dimerState == stateFlag){
			//std::cout<<"transition occured\n"<<std::endl;
			transitionTime = numberOfStepsFromPreviousState * dt;
			numberOfStepsFromPreviousState =0;
			TransitionTime.addSample(transitionTime);
			dimerState = stateFlag;

		}
	}


	int nrNeighbors = 0;
	for (int i = 0; i < NumberOfParticles; i++){

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);
		nrNeighbors=nrNeighbors+neighborIndexI->size();
		//std::cout<<"neighborIndexI->size()"<<neighborIndexI->size()<<std::endl;
	}
	AverageNumberOfNeighbors.addSample(nrNeighbors/NumberOfParticles);

#ifdef COMPUTE_MOMENTA_SQUARE_AVERAGE

	for (int i = 0; i < NumberOfParticles; i++){

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		MomentaSquare.addSample(px*px + py*py+pz*pz);

	}
#endif


#ifdef NVE

	double eCurrent = ComputeEnergy();// ComputeKineticEnergy();
	double errorEnergyPercentage = 100 * abs((eCurrent - intialEnergyNVE)) / abs(intialEnergyNVE);

	EnergyErrorNVE.addSample(errorEnergyPercentage);
#endif

#ifdef BLOCK_AVERAGING

	WritePositionForBlockAveraging();

#endif

	//temperatureForSamples=0.0; //I want to be sure it will be recomputed in pressure if needed

	//std::cout << Particle::howManyRestrained/(double)NumberOfParticles << std::endl;
	//std::cout << (int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles)) << std::endl;

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
	if (doAddSubst == 1){
		restrainedPartPercentageTime[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample(2*((tPerTimeStep1 - tPerTimeStep0) / C_PROC));
		timePerForceUpdate.addSample(2 * ((tPerTimeStep1 - tPerTimeStep0) / C_PROC));
	}
	else	{
		restrainedPartPercentageTime[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample((tPerTimeStep1 - tPerTimeStep0) / C_PROC);
		timePerForceUpdate.addSample( ((tPerTimeStep1 - tPerTimeStep0) / C_PROC));
	}
#endif
#ifdef MEASURE_TIME_PER_TIMESTEP
	restrainedPartPercentageTimePerTimeStep[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample((tPerTimeStep4 - tPerTimeStep3) / C_PROC);
	timePerTimeStep.addSample((tPerTimeStep4 - tPerTimeStep3) / C_PROC);
#endif

	//std::cout << restrainedPartPercentageTime[50].getAverage() << std::endl;
	numberOfInteractionsPerRestrPart[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample(interactionsCounter);

	// norm of force applied on dimer
	double fDimer = 0;
	/*for(int i=0; i<NumberOfDimerParticles; i++){

		fDimer = fDimer+(pParticles[i]->getFX())*(pParticles[i]->getFX()) + (pParticles[i]->getFY())*(pParticles[i]->getFY());

		}
		*/

	// %%%%%%%%%%%%%%%%%%%%%%%%%%% F dimer - || nabla U ||_Linfty ... F dimer only a name.. it is taken over all particles

#ifdef HIST_MAX_NABLA_U

	double maxU = 0.0;
	double pOfmaxU = 0.0;
	double pDelta = pSTEP;
	double pMaxHist = pMAX;
	double maxnablaUHist = NABLA_U_MAX;

	for (int i = 0; i < NumberOfParticles; i++) {

		double m = pParticles[i]->getMass();
		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double epsfi = pParticles[i]->getEpsf();
		double epsri = pParticles[i]->getEpsr();

		double tmpx = 0;
		double tmpy = 0;
		double tmpN = 0;
		double pOfmaxUtmp = 0;

		tmpx = abs(px*(1 - Rho_p(px, py, epsri, epsfi, m)) / m - (px*px + py*py)*((dRho_p(px, py, epsri, epsfi, m)*px)) / (2.0*m));
		tmpy = abs(py*(1 - Rho_p(px, py, epsri, epsfi, m)) / m - (px*px + py*py)*((dRho_p(py, px, epsri, epsfi, m)*py)) / (2.0*m));



		if ((abs(tmpx) < abs(maxnablaUHist)) && (abs(px) < abs(pMaxHist)))	{


			int indexpMax = (int)((pMaxHist + px) / pDelta);
			int indexNablaU = (int)(tmpx / pDelta);

			probaMatrixNablaUCounterTotal = probaMatrixNablaUCounterTotal + 1.0;
			ProbaMatrixNablaU[indexNablaU][indexpMax] = ProbaMatrixNablaU[indexNablaU][indexpMax] + 1.0;


		}



		if (tmpx > tmpy){

			tmpN = tmpx;
			pOfmaxUtmp = px;
		}
		else {

			tmpN = tmpy;
			pOfmaxUtmp = py;

		}

		if (tmpN > maxU){

			maxU = tmpN;
			pOfmaxU = pOfmaxUtmp;

		}
	}

	ForceAppliedOnDimer.addSample(maxU);


	// shift to right- also negative values into the histogram which is made for maximal value pmax

	if (abs(pOfmaxU) < abs(pMaxHist))	{

		HistNablaU[(int)((pMaxHist + pOfmaxU) / pDelta)].addSample(maxU);


		NablaUArray[(int)((pMaxHist + pOfmaxU) / pDelta)] = maxU;

		int pMaxUIndex = (int)((pMaxHist + pOfmaxU) / pDelta);

		for (int i = 0; i < numberOfBinsHistpMax; i++){

			if (i == pMaxUIndex)			NablaUCounterArray[i].addSample(1);
			else							NablaUCounterArray[i].addSample(0);
		}
	}



#endif


#ifdef HISTOGRAM

	double maxValueHistogram = MAX_VALUE_HISTOGRAM;

	for (int i = 0; i < NumberOfParticles; i++){

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();

		//double p2 = px*px + py*py;

		//if (p2 < maxValueHistogram) {

		//	Histogram[(int)(p2 / stepHist)]++;
		//	histogramCounter++;

		//}


		if (px < 0.5* maxValueHistogram) {

			Histogram[(int)((px - 0.5* maxValueHistogram) / stepHist)]++;
			histogramCounter++;

		}
		if (px < 0.5* maxValueHistogram) {

			Histogram[(int)((px - 0.5* maxValueHistogram) / stepHist)]++;
			histogramCounter++;

		}

	}





#endif

#ifdef RDF

	ComputeRadialDistributionFunction();

#endif

#ifdef DISTRIBUTION_RESTR_PART

	ComputeDistributionOfRestrainedParticles();

#endif
};


void Simulation::AddSamplesChangingParameters(){

	RestrainedParticles.addSample(100 * (double)Particle::howManyRestrained / NumberOfParticles);// 100 * percentageRestrainedParticles);

	int nrNeighbors = 0;
	for (int i = 0; i < NumberOfParticles; i++){

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);
		nrNeighbors = nrNeighbors + neighborIndexI->size();
		//std::cout<<"neighborIndexI->size()"<<neighborIndexI->size()<<std::endl;
	}
	AverageNumberOfNeighbors.addSample(nrNeighbors / NumberOfParticles);


	//temperatureForSamples=0.0; //I want to be sure it will be recomputed in pressure if needed

	//std::cout << Particle::howManyRestrained/(double)NumberOfParticles << std::endl;
	//std::cout << (int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles)) << std::endl;
	//restrainedPartPercentageTime[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample(((tPerTimeStep1 - tPerTimeStep0) + (tPerTimeStep4 - tPerTimeStep3)) / C_PROC);
	//std::cout << restrainedPartPercentageTime[50].getAverage() << std::endl;



#ifdef MEASURE_TIME_PER_FORCE_UPDATE
	if (doAddSubst == 1){
		restrainedPartPercentageTime[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample(2 * ((tPerTimeStep1 - tPerTimeStep0) / C_PROC));
	}
	else restrainedPartPercentageTime[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample((tPerTimeStep1 - tPerTimeStep0) / C_PROC);

#endif
#ifdef MEASURE_TIME_PER_TIMESTEP
	restrainedPartPercentageTimePerTimeStep[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample((tPerTimeStep4 - tPerTimeStep3) / C_PROC);
#endif

	numberOfInteractionsPerRestrPart[(int)((100 * (double)Particle::howManyRestrained / (double)NumberOfParticles))].addSample(interactionsCounter);

};

double* Simulation::BlockAveraging(){


	std::cout << "\nBlock averaging...\n" << std::endl;

	int maxLength = BA_count;
	//	std::cout<<"BA_count is "<<BA_count<<std::endl;


	if (BA_count > MAX_ARRAY_LENGTH) std::cout << "Error: BA_count> Max length!\n";

	int pPower = (int)(log(maxLength) / log(2));

	double* varianceBA = new double[pPower]; // standard deviation from block averaging
	double* varianceBA2 = new double[pPower]; // standard deviation from block averaging

#if BA_STYLE ==2

	varianceBA = ComputeBlockAveraging(positionBA_scaled, maxLength);
	varianceBA2 = ComputeBlockAveraging(positionBA_scaled2, maxLength);

#endif

#if BA_STYLE == 1
	varianceBA=ComputeBlockAveraging(positionBA);
#endif

	for (int i = 0; i < pPower; i++) varianceBA[i] = varianceBA[i] * sqrt(WritingPeriodBlockaveraging);
	for (int i = 0; i < pPower; i++) varianceBA2[i] = varianceBA2[i] * sqrt(WritingPeriodBlockaveraging);

	for (int n = 0; n < MAX_ARRAY_LENGTH; n++) {

		//	positionBA[n]=0.0;
		positionBA_scaled[n] = 0.0;
		positionBA_scaled2[n] = 0.0;

	}

	Var_mu = varianceBA[0];


	//std::cout << "From first block averaging value for DimerDistance: Var_mu(A)=" << Var_mu<<std::endl;
	//std::cout << "From first block averaging value for DimerPotential: Var_mu(A)=" << varianceBA2[0] << std::endl;

#ifdef WRITE_BA_RESULT
	// FILE init.
	std::ofstream fileBA("data/BlockAveraging");
	for (int i = 0; i < pPower; i++) fileBA << varianceBA[i] << "\t";
	fileBA <<  "\n";
	for (int i = 0; i < pPower; i++) fileBA << varianceBA2[i] << "\t";
	fileBA << "\n";
#endif

	double sigmaBA = EstimateVarianceFromBA(varianceBA, pPower);
	double sigmaBA2 = EstimateVarianceFromBA(varianceBA2, pPower);

	double* sigmaBAvector = new double[2];
	sigmaBAvector[0] = sigmaBA;
	sigmaBAvector[1] = sigmaBA2;

	return sigmaBAvector;



};





void Simulation::CleanSamples(){

	Temperature.clean();
	Pressure.clean();
	Position.clean();
	Potential.clean();
	RestrainedParticles.clean();
	countPercRestr = 0;
	percentageRestrainedParticles = 0;
	MeasureTime.clean();

	VarianceOne.clean();
	VarianceTwo.clean();
	VarianceThree.clean();

	ForceAppliedOnDimer.clean();
	AutocorrelationAverage.clean();

	NumberOfInteractions.clean();

	for (int i = 0; i < 100; i++) restrainedPartPercentageTime[i].clean();
	for (int i = 0; i < 100; i++) restrainedPartPercentageTimePerTimeStep[i].clean();
	for (int i = 0; i < 100; i++) numberOfInteractionsPerRestrPart[i].clean();


#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

	for (int i = 0; i < MAX_LENGTH_OF_ACF_ARRAY; i++){

		AutocorrelationArrayReplica1[i].clean();
		AutocorrelationArrayReplica2[i].clean();
		AutocorrelationArrayReplica3[i].clean();

		AveragesArrayReplica1[i].clean();
		AveragesArrayReplica2[i].clean();
		AveragesArrayReplica3[i].clean();
	}
#endif

#ifdef VARIANCE_IR

	ReplicasAverageThroughPathN.clean();
	ReplicasAverage.clean();
	VarianceIR.clean();


#endif

#ifdef METROPOLIZATION

	RejectionRate.clean();
	RejectionRateFD.clean();
#endif

	/*
	SimpleAverage<double> Temperature;
	SimpleAverage<double> Pressure;
	SimpleAverage<double> Position;
	SimpleAverage<double> Potential;
	SimpleAverage<double> RestrainedParticles;
	SimpleAverage<double> MeasureTime;
	SimpleAverage<double> ForceAppliedOnDimer;

	SimpleAverage<double> VarianceOne;
	SimpleAverage<double> VarianceTwo;
	SimpleAverage<double> VarianceThree
	SimpleAverage<double> AutocorrelationAverage;
	// this could be variance of dimer force or dimer center... -> SimpleAverage<double> VarianceThree;
	SimpleAverage<double> NumberOfInteractions;

	*/

};


double* Simulation::ComputeBlockAveraging(double *posBA, int NM){


	int pPower;

	pPower = (int)(log(NM) / log(2));

	NM = pow(2, pPower);

	//	std::cout<<"NM="<<NM <<" and SIMULATION_LENGTH "<<SIMULATION_LENGTH<<std::endl;

	double* varianceBlockSize = new double[pPower];

	for (int i = 0; i < pPower; i++) {
		varianceBlockSize[i] = 0.0;

	}


	double* currentBlockAverage = new double[NM];
	for (int i = 0; i < NM; i++) currentBlockAverage[i] = 0.0;

	int k = pPower - 1;

	int M = 1;
	int N = NM / M;
	double A_NM;


	while (M < NM){

		double currentBlockSum = 0.0;
		int indexInBlock_N = 0;
		double sumOfBlockAverages = 0;
		int blockIndex_M = 0;
		double sumVarianceInBlock = 0;


		for (int n = 0; n < NM; n++){

			currentBlockSum = currentBlockSum + posBA[n];
			indexInBlock_N++;

			if (indexInBlock_N == N){

				currentBlockAverage[blockIndex_M] = currentBlockSum;

				currentBlockSum = 0;
				indexInBlock_N = 0;
				blockIndex_M++;

			}
		}

		if (M == 1)    A_NM = currentBlockAverage[0] / N;
		else    {

			for (int m = 0; m < M; m++) sumVarianceInBlock = sumVarianceInBlock + (currentBlockAverage[m] / N - A_NM)*(currentBlockAverage[m] / N - A_NM);
			varianceBlockSize[k] = sqrt(N*sumVarianceInBlock / M);
		}

		k--;
		N = N / 2;
		M = NM / N;

	}

	//	std::cout<<"block averaging done.."<<std::endl;

	//variance block size-> is standard deviation sigma
	return varianceBlockSize;



}



float Simulation::ComputeBoxMuller(float m, float s){

	float invRandMax = 1.0f / RAND_MAX;

	float x1, x2, w, y1;
	static float y2;

	do {


		x1 = 2.0f*((float)randomGenerator->randDouble1()) - 1.0f;
		x2 = 2.0f*((float)randomGenerator->randDouble1()) - 1.0f;

		w = x1*x1 + x2*x2;
	} while (w >= 1.0);

	w = sqrt((-2.0*log(w)) / w);
	y1 = x1*w;
	y2 = x2*w;


	return(m + y1 * s);
};

double Simulation::ComputeEnergy(){

	//return (ComputeKineticEnergy());
	return (ComputePotentialEnergy() + ComputeKineticEnergy());


};

double Simulation::ComputeKineticEnergy(){


	// compute the energy from the particles (N*d degrees of freedom)

	double particlesEnergy = 0;

	for (int i = 0; i < NumberOfParticles; i++) {

		double mass = pParticles[i]->getMass();
		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		particlesEnergy = particlesEnergy + Hp(px, py, pz, mass, epsri, epsfi);//0.5*(px*px + py*py + pz*pz)*(1 - Rho_p(px, py, pz, epsri, epsfi, mass)) / mass;

	}


	return particlesEnergy;


};

double* Simulation::ComputeEnergyFDR(){

	double *kinEnergyPointer = new double[NumberOfParticles];

	for (int i = 0; i < NumberOfParticles; i++){

		double mass = pParticles[i]->getMass();
		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double Gx = pParticles[i]->GetGx();
		double Gy = pParticles[i]->GetGy();
		double Gz = pParticles[i]->GetGz();

		kinEnergyPointer[i] =Hp(px,py,pz,mass,epsri,epsfi) + 0.5*sqrt(kT)* (Gx*Gx + Gy*Gy + Gz*Gz);

	}

	return kinEnergyPointer;

};

double Simulation::ComputeDimerDistance(){

double r =0.0;
for(int i =0; i< NumberOfDimerParticles; i++){
	for(int j =0; j< NumberOfDimerParticles; j++){

		double rx= pParticles[i]->getPositionX()  - pParticles[j]->getPositionX();
		double ry= pParticles[i]->getPositionY()  - pParticles[j]->getPositionY();
		double rz= pParticles[i]->getPositionZ()  - pParticles[j]->getPositionZ();

		r += pow(rx*rx + ry*ry + rz*rz,0.5);
	}
}

return r;


};

double Simulation::ComputePotentialEnergy(){

	double potentialEnergy_tmp = 0.0;
	SolventSolventPotential = 0.0;
	potentialEnergy = 0.0;
	DimerSolventPotential=0.0;
	DimerPotential=0.0;
	DimerDistance=0.0;

#ifdef NO_FORCES
	return potentialEnergy_tmp;

#endif
	for (int i = 0; i < NumberOfParticles; i++){


#ifdef ONE_PARTICLE_ONLY


		double qijx = pParticles[i]->getPositionX();
		double qijy = pParticles[i]->getPositionY();
		double qijz = pParticles[i]->getPositionZ();

		if (PBC == 1){
			qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
			qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
			qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
		}

		double r = sqrt(qijx*qijx + qijy*qijy++ qijz*qijz);

		potentialEnergy_tmp = potentialEnergy_tmp + 0.5*r*r;

#endif
#ifndef ONE_PARTICLE_ONLY
		// solvent type 0, dimer type 1
		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (int j = i + 1; j < NumberOfParticles; j++){



			unsigned int typeParticle_j = pParticles[j]->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double interactionPot = InteractionPotential(r, typeParticle_i, typeParticle_j);
			/*if (tmp != 0)
				int ttt = 123;*/

			if (typeParticle_i + typeParticle_j == 0) 			 {

				SolventSolventPotential = SolventSolventPotential + interactionPot;
			}
			if (typeParticle_i + typeParticle_j == 1) 			 {

				DimerSolventPotential = DimerSolventPotential + interactionPot;
			}
			if (typeParticle_i + typeParticle_j == 2) 			 {

				DimerPotential = DimerPotential + interactionPot;
				DimerDistance = DimerDistance+r;
			}



			potentialEnergy_tmp = potentialEnergy_tmp + interactionPot;


		}

#endif
	}

	potentialEnergy = potentialEnergy_tmp;

	return potentialEnergy_tmp;

};

#ifdef RDF

void Simulation::ComputeRadialDistributionFunction(){

	std::cout<< "RDF does not work in 3d "<<std::endl;

	double rx, ry, r, rxn, ryn, tmp;
	double Q = pow(BoxLength, d);
	double Q1d = BoxLength;

	double cut_off_rdf = pow(Q, 1.0 / d);
	drRDF = cut_off_rdf / numberOfBins;


	int minFitInBox = roundUp(cut_off_rdf / Q1d);
	double rho_box = NumberOfParticles / Q;
	double nr_of_periods = 0;

	if (minFitInBox % 2 == 0){ nr_of_periods = minFitInBox + 1; } //nr of repetitions in one dimension direction
	else { nr_of_periods = minFitInBox; }

	int K1 = (int)(nr_of_periods);
	int K2 = K1;
	int nr_of_boxes = (int)(nr_of_periods*nr_of_periods);

	//std::cout<<"nr of boxes="<<nr_of_boxes<<std::endl;

	for (int i = 0; i<NumberOfParticles; i++) {

		for (int j = 0; j<NumberOfParticles; j++) {

			if (j == i) continue;

			rx = pParticles[j]->getPositionX() - pParticles[i]->getPositionX();
			ry = pParticles[j]->getPositionY() - pParticles[i]->getPositionY();

			for (int k1 = -K1; k1 <= K1; k1++) {

				for (int k2 = -K2; k2 <= K2; k2++) {

					rxn = rx + k1*Q1d;
					ryn = ry + k2*Q1d;

					r = sqrt(rxn*rxn + ryn*ryn);

					if (r<cut_off_rdf) gRDF[(int)(r / drRDF)] += 1.0 / (pow(r, d - 1)*NumberOfParticles*(NumberOfParticles - 1) / Q);

				}

			}

		}

	}



};

#endif


#ifdef DISTRIBUTION_RESTR_PART

void Simulation::ComputeDistributionOfRestrainedParticles(){

	std::cout<< "ComputeDistributionOfRestrainedParticles does not work in 3d "<<std::endl;

	if (NumberOfDimerParticles != 2) {

		std::cout << "ERROR: COMPUTE distribution of restrained particles only for dimer 2!\n";
		while (1);
	}
	else{



		//double Q = pow(BoxLength, d);
		double Q1d = BoxLength;

		double cut_off_drp =  BoxLength/2.0;
		drDRP = cut_off_drp / numberOfBins_DRP;

		double nr_of_periods =  (int)(cut_off_drp / BoxLength)+1;

		//	int minFitInBox = roundUp(cut_off_drp / Q1d);
		//	double rho_box = NumberOfParticles / Q;


		//	if (minFitInBox % 2 == 0){ nr_of_periods = minFitInBox + 1; } //nr of repetitions in one dimension direction
		//	else { nr_of_periods = minFitInBox; }

		int K1 = (int)(nr_of_periods);
		int K2 = K1;
		//int nr_of_boxes = 1;// (int)(nr_of_periods*nr_of_periods);


		double*gDRPhistVector = new double[numberOfBins_DRP];
		double*gDRPhistVectorRestrained = new double[numberOfBins_DRP];

		for (int i = 0; i < numberOfBins_DRP; i++){

			gDRPhistVector[i] = 0;
			gDRPhistVectorRestrained[i] = 0;

		}

		// compute dimer - distance: centered
		double centerDimerX = 0.5*(pParticles[0]->getPositionX() + pParticles[1]->getPositionX());
		double centerDimerY = 0.5*(pParticles[0]->getPositionY() + pParticles[1]->getPositionY());

		// histogram of restrained particles around dimer center

		for (int i = 2; i < NumberOfParticles; i++) {

			double rx = centerDimerX - pParticles[i]->getPositionX();
			double ry = centerDimerY - pParticles[i]->getPositionY();

			//for (int k1 = -K1; k1 <= K1; k1++) {

			//	for (int k2 = -K2; k2 <= K2; k2++) {


			rx = rx - BoxLength*roundMy(rx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
			ry = ry - BoxLength*roundMy(ry / BoxLength);

			//double rxn = rx + k1*Q1d;
			//double ryn = ry + k2*Q1d;

			//	double r = sqrt(rxn*rxn + ryn*ryn);
			double r = sqrt(rx*rx + ry*ry);

			if (r < cut_off_drp) {

				//gDRPhistVector[(int)(r / drDRP)] = +1 / (2 * PI*r*drDRP - PI*drDRP*drDRP)  / (BoxLength*BoxLength);
				gDRPhistVector[(int)(r / drDRP)] = gDRPhistVector[(int)(r / drDRP)] + 1.0 / (pow(r, d - 1)/ pow(BoxLength, 2.0)) / (2.0*3.14159265359*drDRP);
				//if (pParticles[i]->getActiveStatus() == 0) gDRPhistVectorRestrained[(int)(r / drDRP)] = +1 / (2 * PI*r*drDRP - PI*drDRP*drDRP)/(BoxLength*BoxLength);
				if (pParticles[i]->getActiveStatus() == 0) gDRPhistVectorRestrained[(int)(r / drDRP)] = gDRPhistVectorRestrained[(int)(r / drDRP)] + 1.0 / (pow(r, d - 1) / pow(BoxLength, 2.0)) / (2.0*3.14159265359*drDRP);
			}

			//gDRP[(int)(r / drDRP)].addSample(1);  //1.0 / (pow(r, d - 1)*(NumberOfParticles - 1) / Q);
			//	if (r < cut_off_drp) gDRP[(int)(r / drDRP)] += 1.0;// / (pow(r, d - 1)*NumberOfParticles*(NumberOfParticles - 1) / Q);


			//	}

			//	}

		}

		for (int i = 0; i < numberOfBins_DRP; i++) {

			gDRP[i].addSample(gDRPhistVector[i] );
			gDRPRestrained[i].addSample(gDRPhistVectorRestrained[i]);

		}


		delete[] gDRPhistVector;
		delete[] gDRPhistVectorRestrained;
	}


	//double rx, ry, r, rxn, ryn, tmp;
	//double Q = pow(BoxLength, d);
	//double Q1d = BoxLength;

	//double cut_off_drp = pow(Q, 1.0 / d);
	//drDRP = cut_off_drp / numberOfBins_DRP;


	//int minFitInBox = roundUp(cut_off_drp / Q1d);
	//double rho_box = NumberOfParticles / Q;
	//double nr_of_periods = 0;

	//if (minFitInBox % 2 == 0){ nr_of_periods = minFitInBox + 1; } //nr of repetitions in one dimension direction
	//else { nr_of_periods = minFitInBox; }

	//int K1 = (int)(nr_of_periods);
	//int K2 = K1;
	//int nr_of_boxes = 1;// (int)(nr_of_periods*nr_of_periods);

	////std::cout<<"nr of boxes="<<nr_of_boxes<<std::endl;

	//// compute dimer - distance: centered
	//double centerDimerX = 0.5*(pParticles[0]->getPositionX() + pParticles[1]->getPositionX());
	//double centerDimerY = 0.5*(pParticles[0]->getPositionY() + pParticles[1]->getPositionY());

	//// histogram of restrained particles around dimer center

	//for (int i = 2; i < NumberOfParticles; i++) {

	//	if (pParticles[i]->getActiveStatus() == 1 || pParticles[i]->getActiveStatus() == 2) continue;
	//
	//	rx = centerDimerX - pParticles[i]->getPositionX();
	//	ry = centerDimerY - pParticles[i]->getPositionY();

	//		for (int k1 = -K1; k1 <= K1; k1++) {

	//			for (int k2 = -K2; k2 <= K2; k2++) {

	//				rxn = rx + k1*Q1d;
	//				ryn = ry + k2*Q1d;

	//				r = sqrt(rxn*rxn + ryn*ryn);

	//
	//				//if (r < cut_off_drp) gDRP[(int)(r / drDRP)] += 1.0 / (pow(r, d - 1)*NumberOfParticles*(NumberOfParticles - 1) / Q);
	//				if (r < cut_off_drp) gDRP[(int)(r / drDRP)] += 1.0 / (pow(r, d - 1)*(NumberOfParticles - 1) / Q);
	//		//	if (r < cut_off_drp) gDRP[(int)(r / drDRP)] += 1.0;// / (pow(r, d - 1)*NumberOfParticles*(NumberOfParticles - 1) / Q);
	//

	//		}

	//	}

	//}



};

#endif


double Simulation::ComputePressure(){

	double virial = 0;

	for (int i = 0; i < NumberOfParticles; i++){

		virial = +pParticles[i]->getPositionX() * pParticles[i]->getFX() + pParticles[i]->getPositionY() * pParticles[i]->getFY() + pParticles[i]->getPositionZ() * pParticles[i]->getFZ();

	}

	if (temperatureForSamples != 0)	return (temperatureForSamples - virial) / (d*pow(L, d));
	else return (ComputeTemperature() - virial) / (d*pow(L, d));

};

double Simulation::ComputeTemperature(){


	double particlesEnergy = 0;


	// configurational temperature
	/*
	for (int i = 0; i<NumberOfParticles; i++) {

	//double m = pParticles[i]->getMass();
	double qx = pParticles[i]->getPositionX();
	double qy = pParticles[i]->getPositionY();
	double fx = pParticles[i]->getFX();
	double fy = pParticles[i]->getFY();

	particlesEnergy = particlesEnergy + qx*fx +  qy*fy;

	}



	return  particlesEnergy / (NumberOfParticles*d*k);
	*/

	for (int i = 0; i < NumberOfParticles; i++) {

		double mass = pParticles[i]->getMass();

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsfi = pParticles[i]->getEpsf();
		double epsri = pParticles[i]->getEpsr();

		particlesEnergy += px*dHdp(px, py, pz, mass, epsri, epsfi) + py*dHdp(py, px, pz, mass, epsri, epsfi) + pz*dHdp(pz, px, py, mass, epsri, epsfi);
		//(px*px+py*py)*(1-Rho_p(px,py,epsri,epsfi,m))/m- (px*px+py*py)*((dRho_p(px,py,epsri,epsfi,m)*px+(dRho_p(py,px,epsri,epsfi,m)*py)))/(2.0*m);
	}


	return  particlesEnergy / (NumberOfParticles*d*k);



};


double Simulation::ComputeConfigurationalTemperature(){


	double confTemp = 0;


	 //configurational temperature

	for (int i = 0; i<NumberOfParticles; i++) {


	double qx = pParticles[i]->getPositionX();
	double qy = pParticles[i]->getPositionY();
	double qz = pParticles[i]->getPositionZ();

	double fx = pParticles[i]->getFX();
	double fy = pParticles[i]->getFY();
	double fz = pParticles[i]->getFZ();

	confTemp = confTemp - ( qx*fx + qy*fy + qz*fz);

	}

	//std::cout << confTemp / (NumberOfParticles*d*k) << std::endl;

	return  confTemp / (NumberOfParticles*d*k);




};

double Simulation::ComputeTimeLeftToTheEndOfEXE(double s){

	//input average time per time step in seconds
	//return time in hours!!!

	double timeLeft = 0.0;

	double dt_iter = dtMinimalLoop;
	double sumTime = 0.0;
	int currentLoopIndex = 0;

	for (int i = 0; i < nrOfDTloop; i++){

		if (dt == dt_iter) continue;

		currentLoopIndex++;

		dt_iter = dt_iter + deltaTimeStepSize;
	}

	for (int i = currentLoopIndex; i < nrOfDTloop; i++){

		double tmpSeconds = s*(double)((physicalTime) / dt_iter);
		double tmpHours = tmpSeconds / 360;
		sumTime = sumTime + (int)tmpHours;
		dt_iter = dt_iter + deltaTimeStepSize;

		//	std::cout<<"sumTime"<<(double)((physicalTime)/dt_iter)<<std::endl;
	}

	//if(NumberOfReplicas==1)	timeLeft=sumTime;

	//else
	timeLeft = sumTime;



	return timeLeft;
};




double Simulation::EstimateVarianceFromBA(double *v, int lengthv){

	double max = 0;
	int i = 0;
	while (max < v[i] && i < lengthv){

		max = v[i];
		i++;
	}

	return max;

}
void Simulation::DisplayTimeLeft(){

	//if (NumberOfTimeSteps < pow(10.0, 10)){

	if (current_n == 10000 || (current_n != 0 && ((current_n % WritingPeriodOnScreen) == 0))) {

		std::cout << "Current_n:" << current_n << std::endl;
		std::cout << "Current Time Step:" << (double)((double)current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN) << std::endl;

		double currentTime = (NumberOfTimeSteps - current_n)*MeasureTime.getAverage();//(measureTime/current_n);




		long numberOfHours = (long)(currentTime / (60 * 60));
		long numberOfMinutes = (long)((currentTime / 60) - numberOfHours * 60);
		long numberOfSeconds = (long)((currentTime - numberOfHours * 60 * 60 - numberOfMinutes * 60));


		std::cout << "Time left to the end of trajectory: \t" << numberOfHours << "h\t" << numberOfMinutes << "min\t" << numberOfSeconds << "s " << std::endl;// << (double)current_n/NumberOfTimeSteps*100 << "% complete) ";

		double PerformTimeStepTime = MeasureTime.getAverage();// measureTime/current_n;


		//	 double timeToTheEndInHours = ComputeTimeLeftToTheEndOfEXE(PerformTimeStepTime);


		//ConvertTimeUnits(timeToTheEndInSeconds);
		// std::cout << "Time left to the end of exe: \t \t" << timeToTheEndInHours*NumberOfReplicas <<"\t h"<< std::endl;
	}

	//}


};

void Simulation::DisplayTimeLeft(int n, int writingperiod, int N, int measure_t){

	//if (current_n < pow(10.0, 9)){
	if (current_n != 0 && (current_n% writingperiod == 0)) {

		double currentTime = (N - n)*(measure_t / n);
		long numberOfHours = (long)(n / (60 * 60));
		long numberOfMinutes = (long)((n / 60) - numberOfHours * 60);
		long numberOfSeconds = (long)((n - numberOfHours * 60 * 60 - numberOfMinutes * 60));

		std::cout << "Time for the subblock to end: " << numberOfHours << "h\t" << numberOfMinutes << "min\t" << numberOfSeconds << "s " << std::endl;// << (double)current_n/NumberOfTimeSteps*100 << "% complete) ";

		ConvertTimeUnits(currentTimeInSeconds);
	}


	//}

};



void Simulation::ConvertTimeUnits(double seconds){

	seconds = seconds*NumberOfReplicas;
	//std::cout<<NumberOfReplicas<<std::endl;


	long numberOfHours = (long)(seconds / (60 * 60));
	long numberOfMinutes = (long)((seconds / 60) - numberOfHours * 60);
	long numberOfSeconds = (long)((seconds - numberOfHours * 60 * 60 - numberOfMinutes * 60));

	std::cout << "Time left to the end of exe: \t \t" << numberOfHours << "h\t" << numberOfMinutes << "min\t" << numberOfSeconds << "s  " << std::endl;// << (double)current_n/NumberOfTimeSteps*100 << "% complete) ";




};


double Simulation::Hp(double px, double py, double pz, double mass, double epsrp, double epsfp){

#ifdef QUARTIC_POLYNOMIAL

	return pow((px*px + py*py + pz*pz),2.0) *0.25 /mass;
#endif

	double Ux = Hpx(px, mass, epsrp, epsfp);
	double Uy = Hpx(py, mass, epsrp, epsfp);
	double Uz = Hpx(pz, mass, epsrp, epsfp);

	return Ux + Uy + Uz;


};


double Simulation::Hpx(double px, double mass, double epsrp, double epsfp){

#ifdef EXP_PERTUBATION
	double ekin =  0.5*px*px / mass ;
	return ekin + exp(-epsfp * ekin);
#endif



#ifdef SPLINE_FIRST_ORDER

	if (fabs(px) /mass>= epsfp)
	{
		return 0.5*px*px/mass;
	}
	if (fabs(px)/mass <= epsrp)
	{
		double b = epsfp;
		double a = epsrp;

		return 0.5*a*b / mass;
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return -0.5*b*(a*b + 2 * a*x + x *x) / (mass*(a - b));
	}
	else
	{
		double b = epsfp;
		double a = epsrp;
		double x = px;

		return -0.5*b*(a*b - 2 * a*x + x *x) / (mass*(a - b));
	}


#endif

#ifdef SPLINE_SECOND_ORDER


	if (fabs(px)/mass >= epsfp)
	{
		return 0.5*px*px / mass;
	}

	double b = epsfp;
	double a = epsrp;

	double a2 = a*a;
	double a3 = a2*a;
	double b2 = b*b;
	double b3 = b2*b;

	if (fabs(px)/mass <= epsrp)
	{

		return (1 / 6)*(a3 * mass - 3 * a2 * b*mass + 3 * a*b2 * mass - b3 * mass + 4 * a2 * b + a*b2 + b3) / ((a + b)*mass);
	}

	double b4 = b3*b;

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;


		return (1 / 6)*(3 * a2 * b2 * mass + 6 * a2 * b*mass*x + 3 * a2 * mass*x2 - 4 * a*b3 * mass - 6 * a*b2 * mass*x + 2 * a*mass*x3 + b4 * mass - 3 * b2 * mass*x2 - 2 * b*mass*x3 - 3 * a2 * b2 - 6 * a2 * b*x - b4 + 2 * b*x3) / (mass*(a2 - b2));
	}
	else
	{
		//	double b = epsfp;
		//	double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;



		return (1 / 6)*(3 * a2 * b2 * mass - 6 * a2 * b*mass*x + 3 * a2 * mass*x2 - 4 * a*b3 * mass + 6 * a*b2 * mass*x - 2 * a*mass*x3 + b4 * mass - 3 * b2 * mass*x2 + 2 * b*mass*x3 - 3 * a2 * b2 + 6 * a2 * b*x - b4 - 2 * b*x3) / (mass*(a2 - b2));
	}




#endif



#ifdef SPLINE_THIRD_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return 0.5*px*px / mass;
	}
	if (fabs(px)/mass <= epsrp)
	{
		double b = epsfp;
		double a = epsrp;

		return (a*a + 4 * a*b + b*b) / (12 * mass);
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 12)*x*(6 * a3 * x + 24 * a2 * b2 + 6 * a2 * b*x + 8 * a2 * x2 + 24 * a*b2 * x + 8 * a*b*x2 + 3 * a*x3 + 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*b3 * (8 * a2 - a*b - b2) / (mass*(a - b)*(a2 - 2 * a*b + b2));
	}
	else
	{
		double b = epsfp;
		double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 12)*x*(6 * a3 * x - 24 * a2 * b2 + 6 * a2 * b*x - 8 * a2 * x2 + 24 * a*b2 * x - 8 * a*b*x2 + 3 * a*x3 - 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*b3 * (8 * a2 - a*b - b2) / (mass*(a - b)*(a2 - 2 * a*b + b2));
	}


#endif

#ifdef SPLINE_FOURTH_ORDER


	double b = epsfp;
	double a = epsrp;

	double a2 = a*a;
	double b2 = b*b;

	double shift = (a2 + 3 * a*b + b2) / (mass * 10);

	if (fabs(px) / mass >= epsfp)
	{
		return 0.5*px*px / mass;// -shift;
	}
	if (fabs(px) / mass <= epsrp)
	{
		return shift;
	}

	double x = px;
	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x3*x;
	double x5 = x4*x;
	double x6 = x5*x;

	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;


	double b3 = b2*b;
	double b4 = b3*b;
	double b5 = b4*b;
	double b6 = b5*b;
	double b7 = b6*b;

	if (px < 0){


		double spline = (5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 - 60 * a3 * b3 * x - 40 * a3 * b2 * x2 - 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 - 80 * a2 * b2 * x3 - 60 * a2 * b*x4 - 16 * a2 * x5 + 2 * a*b6 - 60 * a*b3 * x3 - 60 * a*b2 * x4 - 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 - 16 * b2 * x5 - 5 * b*x6) / (10 * ((a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4)*(a - b)*mass));

			// (1 / 10)*(5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 - 60 * a3 * b3 * x - 40 * a3 * b2 * x2 - 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 - 80 * a2 * b2 * x3 - 60 * a2 * b*x4 - 16 * a2 * x5 + 2 * a*b6 - 60 * a*b3 * x3 - 60 * a*b2 * x4 - 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 - 16 * b2 * x5 - 5 * b*x6) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

			// (5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 - 60 * a3 * b3 * x - 40 * a3 * b2 * x2 - 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 - 80 * a2 * b2 * x3 - 60 * a2 * b*x4 - 16 * a2 * x5 + 2 * a*b6 - 60 * a*b3 * x3 - 60 * a*b2 * x4 - 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 - 16 * b2 * x5 - 5 * b*x6) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

		//(1 / 12)*x*(6 * a3 * x + 24 * a2 * b2 + 6 * a2 * b*x + 8 * a2 * x2 + 24 * a*b2 * x + 8 * a*b*x2 + 3 * a*x3 + 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*b3 * (8 * a2 - a*b - b2) / (mass*(a - b)*(a2 - 2 * a*b + b2));

		return  spline;// -shift;
		//(1 / 12)*x*(6 * a3 * x + 24 * a2 * b2 + 6 * a2 * b*x + 8 * a2 * x2 + 24 * a*b2 * x + 8 * a*b*x2 + 3 * a*x3 + 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*b3 * (8 * a2 - a*b - b2) / (mass*(a - b)*(a2 - 2 * a*b + b2));
	}
	else
	{


		//double spline = (1 / 10)*(5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 + 60 * a3 * b3 * x - 40 * a3 * b2 * x2 + 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 + 80 * a2 * b2 * x3 - 60 * a2 * b*x4 + 16 * a2 * x5 + 2 * a*b6 + 60 * a*b3 * x3 - 60 * a*b2 * x4 + 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 + 16 * b2 * x5 - 5 * b*x6) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
		double spline = (5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 + 60 * a3 * b3 * x - 40 * a3 * b2 * x2 + 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 + 80 * a2 * b2 * x3 - 60 * a2 * b*x4 + 16 * a2 * x5 + 2 * a*b6 + 60 * a*b3 * x3 - 60 * a*b2 * x4 + 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 + 16 * b2 * x5 - 5 * b*x6) / (10 * ((a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4)*(a - b)*mass));
			//(5 * a5 * x2 - 25 * a4 * b*x2 - 15 * a3 * b4 + 60 * a3 * b3 * x - 40 * a3 * b2 * x2 + 60 * a3 * b*x3 - 15 * a3 * x4 + 4 * a2 * b5 - 90 * a2 * b3 * x2 + 80 * a2 * b2 * x3 - 60 * a2 * b*x4 + 16 * a2 * x5 + 2 * a*b6 + 60 * a*b3 * x3 - 60 * a*b2 * x4 + 28 * a*b*x5 - 5 * a*x6 - b7 - 15 * b3 * x4 + 16 * b2 * x5 - 5 * b*x6) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

		return spline;// -shift;

	}


#endif

};

double Simulation::dHdp(double px, double py, double pz, double mass, double epsrp, double epsfp){

#ifdef QUARTIC_POLYNOMIAL

 return (px*px + py*py + pz*pz) * px /mass;

#endif
	#ifdef EXP_PERTUBATION

			return (px / mass) * (1.0 - epsfp * exp( - 0.5 * px*px / mass))   ;

	#endif

#ifdef SPLINE_FIRST_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return px / mass;
	}
	if (fabs(px)/mass <= epsrp)
	{
		return 0.0;
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return -0.5*b*(2 * a + 2 * x) / (mass*(-b + a));
	}
	else{

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return -0.5*b*(-2 * a + 2 * x) / (mass*(-b + a));


	}
#endif

#ifdef SPLINE_SECOND_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return px / mass;
	}
	if (fabs(px)/mass <= epsrp)
	{
		return 0.0;
	}

	double b = epsfp;
	double a = epsrp;
	double x = px;

	double a2 = a*a;
	double a3 = a2*a;
	double b2 = b*b;
	double b3 = b2*b;

	double x2 = x*x;
	double x3 = x2*x;


	if (px < 0){


		double m = mass;

		return (1 / 6)*(6 * a2 * b*m + 6 * a2 * m*x - 6 * a*b2 * m + 6 * a*m*x2 - 6 * b2 * m*x - 6 * b*m*x2 - 6 * a2 * b + 6 * b*x2) / (m*(a2 - b2));
	}
	else{

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return  (1 / 6)*(-6 * a2 * b*mass + 6 * a2 * mass*x + 6 * a*b2 * mass - 6 * a*mass*x2 - 6 * b2 * mass*x + 6 * b*mass*x2 + 6 * a2 * b - 6 * b*x2) / (mass*(a2 - b2));


	}


#endif


#ifdef SPLINE_THIRD_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return px / mass;
	}
	if (fabs(px)/mass <= epsrp)
	{
		return 0.0;
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 12)*(6 * a3 * x + 24 * a2 * b2 + 6 * a2 * b*x + 8 * a2 * x2 + 24 * a*b2 * x + 8 * a*b*x2 + 3 * a*x3 + 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(6 * a3 + 6 * a2 * b + 16 * a2 * x + 24 * a*b2 + 16 * a*b*x + 9 * a*x2 + 16 * b2 * x + 9 * b*x2) / (mass*(a - b)*(a2 - 2 * a*b + b2));
	}
	else{

		double b = epsfp;
		double a = epsrp;
		double x = px;
		double x2 = x*x;
		double x3 = x2*x;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 12)*(6 * a3 * x - 24 * a2 * b2 + 6 * a2 * b*x - 8 * a2 * x2 + 24 * a*b2 * x - 8 * a*b*x2 + 3 * a*x3 - 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(6 * a3 + 6 * a2 * b - 16 * a2 * x + 24 * a*b2 - 16 * a*b*x + 9 * a*x2 - 16 * b2 * x + 9 * b*x2) / (mass*(a - b)*(a2 - 2 * a*b + b2));


	}
#endif

#ifdef SPLINE_FOURTH_ORDER

	if ( (fabs(px) / mass) >= epsfp)
	{
		return px / mass;
	}
	if ( (fabs(px) / mass) <= epsrp)
	{
		return 0.0;
	}

	double b = epsfp;
	double a = epsrp;
	double x = px;
	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x3*x;
	double x5 = x4*x;

	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;

	double b2 = b*b;
	double b3 = b2*b;
	double b4 = b3*b;

	if (px < 0){



		return (a5 * x - 5 * a4 * b*x - 6 * a3 * b3 - 8 * a3 * b2 * x - 18 * a3 * b*x2 - 6 * a3 * x3 - 18 * a2 * b3 * x - 24 * a2 * b2 * x2 - 24 * a2 * b*x3 - 8 * a2 * x4 - 18 * a*b3 * x2 - 24 * a*b2 * x3 - 14 * a*b*x4 - 3 * a*x5 - 6 * b3 * x3 - 8 * b2 * x4 - 3 * b*x5) / ((a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4)*(a - b)*mass);
			// (a5 * x - 5 * a4 * b*x - 6 * a3 * b3 - 8 * a3 * b2 * x - 18 * a3 * b*x2 - 6 * a3 * x3 - 18 * a2 * b3 * x - 24 * a2 * b2 * x2 - 24 * a2 * b*x3 - 8 * a2 * x4 - 18 * a*b3 * x2 - 24 * a*b2 * x3 - 14 * a*b*x4 - 3 * a*x5 - 6 * b3 * x3 - 8 * b2 * x4 - 3 * b*x5) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
			//(10 * a5 * x - 50 * a4 * b*x - 60 * a3 * b3 - 80 * a3 * b2 * x - 180 * a3 * b*x2 - 60 * a3 * x3 - 180 * a2 * b3 * x - 240 * a2 * b2 * x2 - 240 * a2 * b*x3 - 80 * a2 * x4 - 180 * a*b3 * x2 - 240 * a*b2 * x3 - 140 * a*b*x4 - 30 * a*x5 - 60 * b3 * x3 - 80 * b2 * x4 - 30 * b*x5) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
		//(1 / 12)*(6 * a3 * x + 24 * a2 * b2 + 6 * a2 * b*x + 8 * a2 * x2 + 24 * a*b2 * x + 8 * a*b*x2 + 3 * a*x3 + 8 * b2 * x2 + 3 * b*x3) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(6 * a3 + 6 * a2 * b + 16 * a2 * x + 24 * a*b2 + 16 * a*b*x + 9 * a*x2 + 16 * b2 * x + 9 * b*x2) / (mass*(a - b)*(a2 - 2 * a*b + b2));
	}
	else{



		return (a5 * x - 5 * a4 * b*x + 6 * a3 * b3 - 8 * a3 * b2 * x + 18 * a3 * b*x2 - 6 * a3 * x3 - 18 * a2 * b3 * x + 24 * a2 * b2 * x2 - 24 * a2 * b*x3 + 8 * a2 * x4 + 18 * a*b3 * x2 - 24 * a*b2 * x3 + 14 * a*b*x4 - 3 * a*x5 - 6 * b3 * x3 + 8 * b2 * x4 - 3 * b*x5) / ((a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4)*(a - b)*mass);

			//(a5 * x - 5 * a4 * b*x + 6 * a3 * b3 - 8 * a3 * b2 * x + 18 * a3 * b*x2 - 6 * a3 * x3 - 18 * a2 * b3 * x + 24 * a2 * b2 * x2 - 24 * a2 * b*x3 + 8 * a2 * x4 + 18 * a*b3 * x2 - 24 * a*b2 * x3 + 14 * a*b*x4 - 3 * a*x5 - 6 * b3 * x3 + 8 * b2 * x4 - 3 * b*x5) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
			//(10 * a5 * x - 50 * a4 * b*x + 60 * a3 * b3 - 80 * a3 * b2 * x + 180 * a3 * b*x2 - 60 * a3 * x3 - 180 * a2 * b3 * x + 240 * a2 * b2 * x2 - 240 * a2 * b*x3 + 80 * a2 * x4 + 180 * a*b3 * x2 - 240 * a*b2 * x3 + 140 * a*b*x4 - 30 * a*x5 - 60 * b3 * x3 + 80 * b2 * x4 - 30 * b*x5) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));



	}
#endif

};


double Simulation::dHdp3(double px, double py, double pz, double mass, double epsrp, double epsfp){

#ifdef SPLINE_FIRST_ORDER
	//std::cout << "Error in dHdpNablaLaplacexx: linear interpolation does not have third derivative!!! " << std::endl;
	return 0.0;
#endif

#ifdef SPLINE_SECOND_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return 0.0;
	}
	if (fabs(px)/mass <= epsrp)
	{
		return 0.0;
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return (1 / 6)*(12 * a*mass - 12 * b*mass + 12 * b) / (mass*(a *a - b *b));
	}
	else{

		double b = epsfp;
		double a = epsrp;
		double x = px;

		return (1 / 6)*(-12 * a*mass + 12 * b*mass - 12 * b) / (mass*(a *a - b *b));


	}
#endif


#ifdef SPLINE_THIRD_ORDER

	if (fabs(px)/mass >= epsfp)
	{
		return 0.0;
	}
	if (fabs(px) /mass <= epsrp)
	{
		return 0.0;
	}

	if (px < 0){

		double b = epsfp;
		double a = epsrp;
		double x = px;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 4)*(16 * a2 + 16 * a*b + 18 * a*x + 16 * b2 + 18 * b*x) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(18 * a + 18 * b) / (mass*(a - b)*(a2 - 2 * a*b + b2));
		//(1 / 6)*(12 * a*mass - 12 * b*mass + 12 * b) / (mass*(a *a - b *b));
	}
	else{

		double b = epsfp;
		double a = epsrp;
		double x = px;

		double a2 = a*a;
		double a3 = a2*a;

		double b2 = b*b;
		double b3 = b2*b;

		return (1 / 4)*(-16 * a2 - 16 * a*b + 18 * a*x - 16 * b2 + 18 * b*x) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(18 * a + 18 * b) / (mass*(a - b)*(a2 - 2 * a*b + b2));
		//(1 / 6)*(-12 * a*mass + 12 * b*mass - 12 * b) / (mass*(a *a - b *b));


	}
#endif


#ifdef SPLINE_FOURTH_ORDER

	if (fabs(px) / mass >= epsfp)
	{
		return 0.0;
	}
	if (fabs(px) / mass <= epsrp)
	{
		return 0.0;
	}

	double b = epsfp;
	double a = epsrp;
	double x = px;
	double x2 = x*x;
	double x3 = x2*x;

	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;

	double b2 = b*b;
	double b3 = b2*b;
	double b4 = b3*b;

	if (px < 0){


		return -(12 * (3 * a3 * b + 3 * a3 * x + 4 * a2 * b2 + 12 * a2 * b*x + 8 * a2 * x2 + 3 * a*b3 + 12 * a*b2 * x + 14 * a*b*x2 + 5 * a*x3 + 3 * b3 * x + 8 * b2 * x2 + 5 * b*x3)) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

			//-(12 * (3 * a3 * b + 3 * a3 * x + 4 * a2 * b2 + 12 * a2 * b*x + 8 * a2 * x2 + 3 * a*b3 + 12 * a*b2 * x + 14 * a*b*x2 + 5 * a*x3 + 3 * b3 * x + 8 * b2 * x2 + 5 * b*x3)) / ((a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4)*(a - b)*mass);

			//-(12 * (3 * a3 * b + 3 * a3 * x + 4 * a2 * b2 + 12 * a2 * b*x + 8 * a2 * x2 + 3 * a*b3 + 12 * a*b2 * x + 14 * a*b*x2 + 5 * a*x3 + 3 * b3 * x + 8 * b2 * x2 + 5 * b*x3)) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
			//(1 / 10)*(-360 * a3 * b - 360 * a3 * x - 480 * a2 * b2 - 1440 * a2 * b*x - 960 * a2 * x2 - 360 * a*b3 - 1440 * a*b2 * x - 1680 * a*b*x2 - 600 * a*x3 - 360 * b3 * x - 960 * b2 * x2 - 600 * b*x3) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

	}
	else{



		return (12 * (3 * a3 * b - 3 * a3 * x + 4 * a2 * b2 - 12 * a2 * b*x + 8 * a2 * x2 + 3 * a*b3 - 12 * a*b2 * x + 14 * a*b*x2 - 5 * a*x3 - 3 * b3 * x + 8 * b2 * x2 - 5 * b*x3)) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));

			//(12 * (3 * a3 * b - 3 * a3 * x + 4 * a2 * b2 - 12 * a2 * b*x + 8 * a2 * x2 + 3 * a*b3 - 12 * a*b2 * x + 14 * a*b*x2 - 5 * a*x3 - 3 * b3 * x + 8 * b2 * x2 - 5 * b*x3)) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
			//(1 / 10)*(360 * a3 * b - 360 * a3 * x + 480 * a2 * b2 - 1440 * a2 * b*x + 960 * a2 * x2 + 360 * a*b3 - 1440 * a*b2 * x + 1680 * a*b*x2 - 600 * a*x3 - 360 * b3 * x + 960 * b2 * x2 - 600 * b*x3) / ((a - b)*mass*(a4 - 4 * a3 * b + 6 * a2 * b2 - 4 * a*b3 + b4));
		//			(1 / 4)*(-16 * a2 - 16 * a*b + 18 * a*x - 16 * b2 + 18 * b*x) / (mass*(a - b)*(a2 - 2 * a*b + b2)) + (1 / 12)*x*(18 * a + 18 * b) / (mass*(a - b)*(a2 - 2 * a*b + b2));
		//(1 / 6)*(-12 * a*mass + 12 * b*mass - 12 * b) / (mass*(a *a - b *b));


	}
#endif

};


double Simulation::dHdp2xx(double px, double py, double pz,double mass, double epsrp, double epsfp){


	double invMass = 1.0 / mass;

	return  invMass*(1 - Rho_p(px, py, pz, epsrp, epsfp, mass) - 2 * px*dRho_p2xx(px, py, pz, epsrp, epsfp, mass) - 0.5*(px*px + py*py + pz*pz)*dRho_p2xx(px, py, pz, epsrp, epsfp, mass));

		//px*(1.0 - Rho_p(px, py, epsrp, epsfp, mass)) / mass - px*(-dRho_p(px, py, epsrp, epsfp, mass)*px / mass) / mass - (px)*dRho_p(px, py, epsrp, epsfp, mass) / (mass)-(px*px)*dRho_p2(px, py, epsrp, epsfp, mass) / (2.0*mass)*px / mass;

};

double Simulation::dHdp2xy(double px, double py, double pz, double mass, double epsrp, double epsfp){

	double invMass = 1.0 / mass;

	return  invMass*(-py*dRho_p(px, py, pz, epsrp, epsfp, mass) - px*dRho_p(py, px, pz, epsrp, epsfp, mass) - 0.5*(px*px + py*py + pz*pz)*dRho_p2xy(px, py, pz, epsrp, epsfp, mass));

	//px*(1.0 - Rho_p(px, py, epsrp, epsfp, mass)) / mass - px*(-dRho_p(px, py, epsrp, epsfp, mass)*px / mass) / mass - (px)*dRho_p(px, py, epsrp, epsfp, mass) / (mass)-(px*px)*dRho_p2(px, py, epsrp, epsfp, mass) / (2.0*mass)*px / mass;

};

double	Simulation::dHdpNablaLaplacexx(double px, double py, double pz, double mass, double epsrp, double epsfp){

#ifdef SPLINE_FIRST_ORDER
	std::cout << "Error in dHdpNablaLaplacexx: linear interpolation does not have third derivative!!! " << std::endl;
#endif
	// du^3 U in direction px



	return dHdp3(px, py, pz, mass, epsrp, epsfp);
		//invMass*(-3*dRho_p(px,py,pz,epsrp,epsfp,mass)-3*px*dRho_p2xx(px,py,pz,epsrp,epsfp,mass)-0.5*(px*px+py*py+pz*pz)*dRho_p3xx(px,py,pz,epsrp,epsfp,mass));

}
double	Simulation::dHdpNablaLaplacexy(double px, double py, double pz, double mass, double epsrp, double epsfp){

	//dpy dp^2xx

	//double invMass = 1.0 / mass;

	return 0;//
	//invMass*(-dRho_p(py, px, pz, epsrp, epsfp, mass) - 2 * py * dRho_p2xx(px, py, pz, epsrp, epsfp, mass)-0.5*(px*px+py*py+pz*pz)*dRho_p3xy(px,py,pz,epsrp,epsfp,mass));


}





double Simulation::HkineticEnergyFunction(double px, double py, double pz, double mass, double epsrp, double epsfp){


#ifndef OTHER_KINETIC_ENERGY

	return  (px*px + py*py+pz*pz)*(1 - dRho_p(px, py, pz,epsrp, epsfp, mass)) / (2.0*mass);
#endif

#ifdef OTHER_KINETIC_ENERGY

	return   (px*px + py*py+pz*pz)* (px*px + py*py+pz*pz)/4.0;


#endif

};

void Simulation::SaveCurrentPositionsAsOld(){

	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->setPositionOld();
		//remember current positions

	}

};


void Simulation::SaveCurrentSimulation(){


	std::ofstream filePos("CURRENT_POSITIONS");
	std::ofstream fileMom("CURRENT_MOMENTA");

#ifdef MEASURE_TIME_PER_TIMESTEP
	std::ofstream fileTime("CURRENT_TIME_PER_TIMESTEP");
#endif
	//*** save positions

	for (int i = 0; i < NumberOfParticles; i++){

		filePos << pParticles[i]->getPositionX() << " ";

	}

	filePos << std::endl;

	for (int i = 0; i < NumberOfParticles; i++){

		filePos << pParticles[i]->getPositionY() << " ";

	}

	for (int i = 0; i < NumberOfParticles; i++){

		filePos << pParticles[i]->getPositionZ() << " ";

	}

	filePos << std::endl;

	//*** save momenta

	for (int i = 0; i < NumberOfParticles; i++){

		fileMom << pParticles[i]->getMomentumX() << " ";

	}

	fileMom << std::endl;

	for (int i = 0; i < NumberOfParticles; i++){

		fileMom << pParticles[i]->getMomentumY() << " ";

	}

	for (int i = 0; i < NumberOfParticles; i++){

		fileMom << pParticles[i]->getMomentumZ() << " ";

	}

	fileMom << std::endl;


#ifdef BLOCK_AVERAGING
	double* sigmaBAtmp = BlockAveraging();
	double sigmaBA1 = sigmaBAtmp[0]; // dimerDistance variance
	double sigmaBA2 = sigmaBAtmp[1]; // dimerPotential variance

	//	fileAver << sigmaBA1 << std::endl;
	//	fileAver << sigmaBA2 << std::endl;

#endif

#ifdef MEASURE_TIME_PER_TIMESTEP
	// save average time per time step

	for (int i = 0; i < 100; i++) fileTime<< restrainedPartPercentageTime[i].getAverage()<<"\t";
#endif



};

#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS

void Simulation::SaveCurrentAverages(){


	int aveTmp = HOW_MANY_PHYSICAL_TIMES;

	if (AverageCounter > aveTmp){

		std::cout << "Number of averages per subinterval no correct!" << std::endl;

	}


	std::ofstream fileAver("CURRENT_AVERAGES_EVOLUTION.m");


	AverageMatrix[0][AverageCounter] = (double)(current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN);
	AverageMatrix[1][AverageCounter] = Temperature.getAverage();
	AverageMatrix[2][AverageCounter] = Pressure.getAverage();
	AverageMatrix[3][AverageCounter] = Position.getAverage();
	AverageMatrix[4][AverageCounter] = Potential.getAverage();
	AverageMatrix[5][AverageCounter] = RestrainedParticles.getAverage();
	AverageMatrix[6][AverageCounter] = VarianceOne.getAverage();
	AverageMatrix[7][AverageCounter] = VarianceTwo.getAverage();
	AverageMatrix[8][AverageCounter] = VarianceThree.getAverage();
	AverageMatrix[9][AverageCounter] = ForceAppliedOnDimer.getAverage();
	AverageMatrix[10][AverageCounter] = NumberOfInteractions.getAverage();

	fileAver << "% current time step, temperature, pressure, position (dimer distance), potential V, percentage restr particles, observable 1 (solvent solvent pot), observable 2 (dimer potential), observable 3 (dimer solvent pot), force applied on dimer,number of interactions" << std::endl;
	fileAver << "A=[\n";

	for (int j = 0; j <= AverageCounter; j++){

		fileAver << AverageMatrix[0][j] << "\t";
		fileAver << AverageMatrix[1][j] << "\t";
		fileAver << AverageMatrix[2][j] << "\t";
		fileAver << AverageMatrix[3][j] << "\t";
		fileAver << AverageMatrix[4][j] << "\t";
		fileAver << AverageMatrix[5][j] << "\t";
		fileAver << AverageMatrix[6][j] << "\t";
		fileAver << AverageMatrix[7][j] << "\t";
		fileAver << AverageMatrix[8][j] << "\t";
		fileAver << AverageMatrix[9][j] << "\t";
		fileAver << AverageMatrix[10][j] << "\t";
		fileAver << "\n";
	}

	fileAver << "];\n";
	AverageCounter = AverageCounter + 1;
};



void Simulation::SaveCurrentAveragesOnlyInstantOnes(){

#ifdef RUN_REFERENCE_SIMULATION_WITH_STD
	std::ofstream fileAver("CURRENT_AVERAGES_STD.m");
#endif
#ifndef RUN_REFERENCE_SIMULATION_WITH_STD
	std::ofstream fileAver("CURRENT_AVERAGES.m");
#endif
	fileAver << "% current time step, temperature, pressure, position (dimer distance), potential V, percentage restr particles, observable 1 (solvent solvent pot), observable 2 (dimer potential), observable 3 (dimer solvent pot), config temperature,number of interactions" << std::endl;
	fileAver << "A=[\n";

	fileAver << (double)(current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN) << "\t";
	fileAver << Temperature.getAverage() << "\t";
	fileAver << Pressure.getAverage() << "\t";
	fileAver << Position.getAverage() << "\t";
	fileAver << Potential.getAverage() << "\t";
	fileAver << RestrainedParticles.getAverage() << "\t";
	fileAver << VarianceOne.getAverage() << "\t";
	fileAver << VarianceTwo.getAverage() << "\t";
	fileAver << VarianceThree.getAverage() << "\t";
//	fileAver << ConfigTemperature.getAverage() << "\t";
	fileAver << NumberOfInteractions.getAverage() << "\t";



#ifdef METROPOLIZATION
	fileAver << RejectionRate.getAverage() << "\t";
	fileAver << RejectionRateFD.getAverage() << "\t";

#endif

	fileAver << TransitionTime.getAverage() << "\t";

	fileAver << "\n";


	fileAver << "];\n";


#ifdef WRITE_AVER_TIME_PER_TIME_STEP


	std::ofstream fileAverTimeForce("TIME_PER_FORCEUPDATE_inst.m");
	fileAverTimeForce << "timePerForceUpdate=" << timePerForceUpdate.getAverage() << ";\n" << std::endl;

	fileAverTimeForce << "timePerTimeStep=[" << std::endl;

	for (int i = 0; i < 100; i++) fileAverTimeForce << restrainedPartPercentageTime[i].getAverage() << "\t";
	fileAverTimeForce << "];" << std::endl;

	fileAverTimeForce << "nrOfInteractions=[" << std::endl;
	for (int i = 0; i < 100; i++) fileAverTimeForce << numberOfInteractionsPerRestrPart[i].getAverage() << "\t";
	fileAverTimeForce << "];" << std::endl;

	//***********************************************

	std::ofstream fileAverTime("TIME_PER_TIMESTEP_inst.m");
	fileAverTime << "timePerTimeStep=" << timePerTimeStep.getAverage() << ";\n" << std::endl;

	fileAverTime << "timePerTimeStepInPercRestr=[" << std::endl;

	for (int i = 0; i < 100; i++) fileAverTime << restrainedPartPercentageTimePerTimeStep[i].getAverage() << "\t";
	fileAverTime << "];" << std::endl;


#endif
};
#endif

void Simulation::SaveCurrentMomentaAsOld(){

	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->setMomentumOld();
		//remember current momenta

	}

};


void Simulation::setCurrentTimeInSeconds(double seconds){

	currentTimeInSeconds = seconds;
	ConvertTimeUnits(currentTimeInSeconds);

};


void Simulation::setNumberStepsFORloop(double physTime, double dtMinimal, int NumberOfDifferentDt, double deltaDT){


	physicalTime = physTime;
	dtMinimalLoop = dtMinimal;
	nrOfDTloop = NumberOfDifferentDt;
	deltaTimeStepSize = deltaDT;


};

void Simulation::setSimulationParameters(double NrOfTimeSteps, double TimeStepSize, double er, double ef, int nrDimer, int NrRepl, double initialDensity, double nrEquil){

#ifdef GAMMA_ANALYSIS

	gamma = 1;// dt * 1000;

	std::cout << "\nGamma = " << gamma << std::endl;

#endif

	NumberOfTimeSteps = NrOfTimeSteps;
	NumberOfDimerParticles = nrDimer;
	NumberOfReplicas = NrRepl;
	physicalTime = NrOfTimeSteps*TimeStepSize;
	NumberOfEquilibrationSteps = nrEquil;
	dt = TimeStepSize;

	// BoxLength must be chosen as a multiple of cutoff because of grid cells !

	// choose boxlength according to density
	BoxLength = pow(((double)(NumberOfParticles - 2) / initialDensity),1.0/d);
	//adjust it according to cutoff -> better bigger to avoid explosion (...+1)
	BoxLength = cut_off*(round(BoxLength / cut_off) + 1);

	std::cout << "cut_off = "<< cut_off <<"\t BoxLength = "<< BoxLength << std::endl;

// #ifdef 	STANDARD_UPDATE_FORCES
// 	if (er > 0 || ef > 0){
// 		std::cout << "\n \n \nMessage from setSimulationParameters: \nEpsR and EpsF must be set zero when using STANDARD_UPDATE_FORCES!" << std::endl;
// 		while (1);
// 	}
// #endif

	// this is the old version where the boxlength was adjusted for given density:
	//if (NumberOfParticles > 2)	BoxLength = sqrt((double)(NumberOfParticles - 2) / initialDensity);//sqrt(NumberOfParticles*k*T / alpha);
	//else BoxLength = 2 * cut_off + 1;
	//if (BoxLength < 2 * cut_off) {
	//	std::cout << "In simulation parameters: boxlength must be bigger than 2*cut off!" << std::endl;
	//	while (1);
	//}





	//// BoxLength must be chosen as a multiple of cutoff because of grid cells !
	//// choose boxlength according to density
	//BoxLength = pow(((double)(NumberOfParticles - 2) / initialDensity), 1.0 / d);
	////adjust it according to cutoff -> better bigger to avoid explosion (...+1)
	//BoxLength = cut_off*(round(BoxLength / cut_off) + 1);


	std::cout << "cut_off = " << cut_off << "\t BoxLength = " << BoxLength << std::endl;

	std::cout << "Epsilon parameters chosen for dimer and solvent separately"<< std::endl;
	for (int i = 0; i < NumberOfParticles; i++){

		if (i < NumberOfDimerParticles){

			epsr[i] = 0.0;// ((double)i) / 10;//0;
			epsf[i] = er;// ((doublFe)i) * 3 / 10;// 0;

		}
		else{


			epsr[i] = 0.0;// ((double)i) / 10;//er;
			epsf[i] = ef;// ((double)i) * 3 / 10;// ef;
		}

	}


// uncomment for ARPS
	// for (int i = 0; i < NumberOfParticles; i++){
	//
	// 	if (i < NumberOfDimerParticles){
	//
	// 		epsr[i] = 0;// ((double)i) / 10;//0;
	// 		epsf[i] = 0;// ((double)i) * 3 / 10;// 0;
	//
	// 	}
	// 	else{
	//
	//
	// 		epsr[i] = er;// ((double)i) / 10;//er;
	// 		epsf[i] = ef;// ((double)i) * 3 / 10;// ef;
	// 	}
	//
	// }

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

	double physicalTimeACF1 = TIME_AUTOCORRELATION_1;
	double physicalTimeACF2 = TIME_AUTOCORRELATION_2;
	double physicalTimeACF3 = TIME_AUTOCORRELATION_3;

	NumberOfStepsInTimeInterval1 = (int)(physicalTimeACF1 / dt);
	NumberOfStepsInTimeInterval2 = (int)(physicalTimeACF2 / dt);
	NumberOfStepsInTimeInterval3 = (int)(physicalTimeACF3 / dt);

	std::cout << "Variance by Autocorrelation function " << std::endl;
	std::cout << "Observable 1: physical time = " << physicalTimeACF1 << std::endl;
	std::cout << "Observable 2: physical time = " << physicalTimeACF2 << std::endl;
	std::cout << "Observable 3: physical time = " << physicalTimeACF3 << std::endl;
	std::cout << "Observable 1: nr of steps = " << NumberOfStepsInTimeInterval1 << std::endl;
	std::cout << "Observable 2: nr of steps = " << NumberOfStepsInTimeInterval2 << std::endl;
	std::cout << "Observable 3: nr of steps = " << NumberOfStepsInTimeInterval3 << std::endl;
	std::cout << std::endl;


#endif


#ifdef VARIANCE_IR

	LengthOfPartialTrajectoryN = (int)(LENGTH_N_IN_REPLICAS / dt);
	std::cout << "Variance by Independent Replicas " << std::endl;
	std::cout << "Observable 2 only: physical time (from acf corr time) = " << LENGTH_N_IN_REPLICAS << std::endl;
	std::cout << "Number of steps N = " << LengthOfPartialTrajectoryN << std::endl;

#endif

#ifdef MORE_VARIANCE_IR

	LengthOfPartialTrajectoryN2 = (int)(LENGTH_N_IN_REPLICAS2 / dt);
	std::cout << "Variance by Independent Replicas " << std::endl;
	std::cout << "Observable 2 only: physical time (from acf corr time) = " << LENGTH_N_IN_REPLICAS << std::endl;
	std::cout << "Number of steps N = " << LengthOfPartialTrajectoryN << std::endl;

#endif


#ifdef NO_FORCES
	std::cout << "\n  ********* NO FORCES ************** " << std::endl;
#endif

	std::cout << "\n************** SIMULATION PARAMETERS ************** " << std::endl;
	std::cout << "dt=  " << dt << std::endl;
	std::cout << "Number of time steps= " << NumberOfTimeSteps << std::endl;
	std::cout << "Number of equilibration time steps= " << NumberOfEquilibrationSteps << std::endl;
	std::cout << "Number of replicas= " << NumberOfReplicas << std::endl;
	//std::cout << "Number of particles= " << NumberOfParticles << " from which number of dimer= " << NumberOfDimerParticles << std::endl;
	std::cout << "BoxLength= " << BoxLength << " in dimension  " << d << std::endl;
	std::cout << "Force parameters: \n cut_off=" << cut_off << "\n Sigma_LJ= " << sigma_LJ << "\t Epsilon_LJ=" << epsilon_LJ << std::endl;
	std::cout << "Dimer: hDimer=" << hDimer << "\t wDimer= " << wDimer << std::endl;
	//std::cout << "Density=" << initialDensity << std::endl;
	std::cout << "Density " << NumberOfParticles / pow(BoxLength, d) << std::endl;
	std::cout << "***************************************************" << std::endl;


	BA_count = 0;
	measureTime = 0.0;
	current_n = 0;

};

double** Simulation::readMatrixFromFile(int nrRows, int nrColumns){

	std::cout << "reading matrix from file.." << std::endl;

	std::ifstream file("MATRIX");

	double ** array_2d = new double*[nrRows];


#ifdef READ_INITIAL_CONDITION_FROM_FILE

	//  0  0.01   0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000		0.8000    0.9000    1.0000
	double* epsArray = new double[11];
	epsArray[0] = 0;
	for (int i = 1; i < 11; i++) {
		epsArray[i] = epsArray[i - 1] + 0.1;


	}


	if (((double)abs(epsr[NumberOfDimerParticles] - 0.01) < pow(10, -14))){
		std::ifstream file("MATRIX_0d01");
		if (!file) std::cout << "error file not found" << std::endl;
		std::cout << "reading matrix for eps= " << 0.01 << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[0]) < pow(10, -14))){
		std::ifstream file("MATRIX");
		if (!file) std::cout << "error file not found" << std::endl;
		std::cout << "reading matrix for eps= " << epsArray[0] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[1]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_2");
		std::cout << "reading matrix for eps= " << epsArray[1] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[2]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_3");
		std::cout << "reading matrix for eps= " << epsArray[2] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[3]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_4");
		std::cout << "reading matrix for eps= " << epsArray[3] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[4]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_5");
		std::cout << "reading matrix for eps= " << epsArray[4] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[5]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_6");
		std::cout << "reading matrix for eps= " << epsArray[5] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[6]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_7");
		std::cout << "reading matrix for eps= " << epsArray[6] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[7]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_8");
		std::cout << "reading matrix for eps= " << epsArray[7] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[8]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_9");
		std::cout << "reading matrix for eps= " << epsArray[8] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[9]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_10");
		std::cout << "reading matrix for eps= " << epsArray[9] << std::endl;
	}
	if (((double)abs(epsr[NumberOfDimerParticles] - epsArray[10]) < pow(10, -14))){
		std::ifstream file("MATRIX_ARPS_11");
		std::cout << "reading matrix for eps= " << epsArray[10] << std::endl;

	}
	/*else{
		std::ifstream file("MATRIX");
		std::cout << "reading matrix for eps= " << epsArray[0] << std::endl;
		}
		*/

#endif

	if (!file){
		std::cout << "Error: cannot find file MATRIX!" << std::endl;
		while (1);
	}

	for (unsigned int i = 0; i < nrRows; i++) {
		array_2d[i] = new double[nrColumns];

		for (unsigned int j = 0; j < nrColumns; j++) {
			file >> array_2d[i][j];
		}
	}





	/*for (unsigned int i = 0; i < nrRows; i++) {


		for (unsigned int j = 0; j < nrColumns; j++) {
		std::cout << array_2d[i][j];
		}
		std::cout << std::endl;
		}*/

	return array_2d;


};


int Simulation::roundMy(double d1)
{
	if (d1 >= 0) {
		if (d1 - ((int)d1) < 0.5)
			return (int)d1;
		return (int)d1 + 1;
	}
	else {
		if (d1 - ((int)d1) > -0.5)
			return (int)d1;
		return (int)d1 - 1;
	}
}


int Simulation::roundUp(double x)
{
	return (int)x + 1;
}

double Simulation::Rho_p(double px, double py, double pz, double epsrp, double epsfp, double mass){

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta3 = 0;
	double eta4 = 0;
	double eta5 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp)
	{


#ifdef KINETIC_ENERGY_DOMAIN_COUNTERS
		//counters for domain identification - in which part of the kinetic energy I am the most of the time
		FullDynamicsDomain.addSample(1.0);
		RestrainedDynamicsDomain.addSample(0.0);
		SplineDomain.addSample(0.0);
#endif


		return 0.0;
	}
	if (e_kin <= epsrp)
	{


#ifdef KINETIC_ENERGY_DOMAIN_COUNTERS
		//counters for domain identification - in which part of the kinetic energy I am the most of the time
		FullDynamicsDomain.addSample(0.0);
		RestrainedDynamicsDomain.addSample(1.0);
		SplineDomain.addSample(0.0);
#endif

		return 1.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta3 = eta*eta*eta;
	eta4 = eta3*eta;
	eta5 = eta4*eta;



#ifdef KINETIC_ENERGY_DOMAIN_COUNTERS
	//counters for domain identification - in which part of the kinetic energy I am the most of the time
	FullDynamicsDomain.addSample(0.0);
	RestrainedDynamicsDomain.addSample(0.0);
	SplineDomain.addSample(1.0);
#endif

	return (-6.0*eta5 + 15.0*eta4 - 10.0*eta3 + 1.0);

};



double Simulation::dRho_p(double px, double py, double pz, double epsrp, double epsfp, double mass){

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta2 = 0;
	double eta3 = 0;
	double eta4 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp || e_kin <= epsrp)
	{
		return 0.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta2 = eta*eta;
	eta3 = eta2*eta;
	eta4 = eta3*eta;

	return (-30.0*eta4 + 60.0*eta3 - 30.0*eta2)*(1 / delta)*(px / (mass));

};


double Simulation::dRho_p2xx(double px, double py, double pz, double epsrp, double epsfp, double mass){

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta2 = 0;
	double eta3 = 0;
	double eta4 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp || e_kin <= epsrp)
	{
		return 0.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta2 = eta*eta;
	eta3 = eta2*eta;
	eta4 = eta3*eta;

	double dK = ((1 / delta)*(px / (mass)));

	return (-30.0 * 4 * eta3 + 3 * 60.0*eta2 - 2 * 30.0*eta)*dK*dK + (-30.0*eta4 + 60.0*eta3 - 30.0*eta2)/mass;

};


double Simulation::dRho_p2xy(double px, double py, double pz, double epsrp, double epsfp, double mass){

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta2 = 0;
	double eta3 = 0;
	double eta4 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp || e_kin <= epsrp)
	{
		return 0.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta2 = eta*eta;
	eta3 = eta2*eta;
	eta4 = eta3*eta;

	double dK = px*py/(delta*delta*mass*mass);

	return (-30.0 * 4 * eta3 + 3 * 60.0*eta2 - 2 * 30.0*eta)*dK;

};


double Simulation::dRho_p3xx(double px, double py, double pz, double epsrp, double epsfp, double mass){

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta2 = 0;
	double eta3 = 0;
	double eta4 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp || e_kin <= epsrp)
	{
		return 0.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta2 = eta*eta;
	eta3 = eta2*eta;
	eta4 = eta3*eta;

	double dK = px/ (delta*mass);
	double dK2 = px/(delta*mass*delta*mass);
	double dK3 = dK*dK*dK;


	return (-3*30.0 * 4 * eta2 + 2*3 * 60.0*eta - 2 * 30.0)*dK3 +3*(-30.0 * 4 * eta3 + 3 * 60.0*eta2 - 2 * 30.0*eta)*dK2;

};

double Simulation::dRho_p3xy(double px, double py, double pz, double epsrp, double epsfp, double mass){

	//dpy d^2dpx

	double delta = epsfp - epsrp;
	double eta = 0;
	double eta2 = 0;
	double eta3 = 0;
	double eta4 = 0;

	double e_kin = (px*px + py*py + pz*pz) / (2.0*mass);

	if (e_kin >= epsfp || e_kin <= epsrp)
	{
		return 0.0;
	}
	eta = (e_kin - epsrp) / delta;

	eta2 = eta*eta;
	eta3 = eta2*eta;
	eta4 = eta3*eta;

	double dm = 1.0 / (delta*mass);
	double dm2 = dm*dm;
	double dm3 = dm2*dm;

	double dK3 = px*px*py *dm3;
	double dK2 = py *dm2;

	return (-3 * 30.0 * 4 * eta2 + 2 * 3 * 60.0*eta - 2 * 30.0)*dK3 +  (-30.0 * 4 * eta3 + 3 * 60.0*eta2 - 2 * 30.0*eta)*dK2;

};


double	Simulation::getCurrentTimeInSeconds() const{

	return currentTimeInSeconds;

};

double* Simulation::getArrayTimePerTimeStep() const{

	double* averTimePerTimeStepDependingOnPercRestr = new double[100];
	for (int i = 0; i < 100; i++) averTimePerTimeStepDependingOnPercRestr[i] = restrainedPartPercentageTime[i].getAverage();
	return averTimePerTimeStepDependingOnPercRestr;

};

double Simulation::getObservableForVariance(int whichVariance){

	if (whichVariance == 1)

		return VarianceOne.getAverage();
	if (whichVariance == 2)
		return VarianceTwo.getAverage();
	if (whichVariance == 3)
		return VarianceThree.getAverage();
}

double* Simulation::getSamples() {

	double *samples = new double[NUMBER_OF_AVERAGES];

	/*samples[0]=Temperature.getAverage();
	samples[1]=Pressure.getAverage();
	samples[2]=ForceAppliedOnDimer.getAverage();
	samples[3]=Potential.getAverage();
	samples[4]=standardDeviation;
	samples[5] = standardDeviation2;
	samples[6]=RestrainedParticles.getAverage();
	samples[7]=VarianceOne.getAverage();
	samples[8] = VarianceTwo.getAverage();
	*/
	samples[0] = (double)(current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN);
	samples[1] = Temperature.getAverage();
	samples[2] = Pressure.getAverage();
	samples[3] = Position.getAverage();
	samples[4] = Potential.getAverage();
	samples[5] = RestrainedParticles.getAverage();
	samples[6] = VarianceOne.getAverage();
	samples[7] = VarianceTwo.getAverage();
	samples[8] = VarianceThree.getAverage();
#ifdef COMPUTE_MOMENTA_SQUARE_AVERAGE
	samples[9] = MomentaSquare.getAverage();
#endif
#ifndef COMPUTE_MOMENTA_SQUARE_AVERAGE
	samples[9] = ForceAppliedOnDimer.getAverage();
#endif
	samples[10] = NumberOfInteractions.getAverage();


	std::cout << "\nNumber of interactions is " << NumberOfInteractions.getAverage() << std::endl;

	return samples;
};

double Simulation::getVar_mu() const{

	return Var_mu;
};

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION
double Simulation::getStandardDeviation() {

	//double AAverEnd2 = VarianceTwo.getAverage();

	//file2 << "A2=[";

	double A2 = 0;
	for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {

		double tmp = (AutocorrelationArrayReplica2[i].getAverage() - AveragesArrayReplica2[i].getAverage()*AveragesArrayReplica2[0].getAverage());;

		if (i == 0 || i == NumberOfStepsInTimeInterval2 - 1)
			A2 = A2 + 0.5*tmp;
		else
			A2 = A2 + tmp;

	}
	A2 = 2 * dt*A2;


	/*file2 << "AAver2=[";
	for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {

		file2 << AveragesArrayReplica2[i].getAverage() << " ";
	}
	file2 << "\n];" << std::endl;
*/

	/*Variance1B = 2 * dt*trapz(A1 - AAver1.*AAver1(1));
	Variance2B = 2 * dt*trapz(A2 - AAver2.*AAver2(1));
	Variance3B = 2 * dt*trapz(A3 - AAver3.*AAver3(1));
*/

	double standardDeviation = sqrt(A2);
	return standardDeviation;


};

#endif

double Simulation::getValues(int i) const {

	switch (i){

	case 1: return NumberOfParticles;

	case 2: return NumberOfTimeSteps;

	case 3: return dt;


	}


};

void Simulation::InitialCondition(){

	errorIndicator = false;

	Particle* ptr;
	pParticles = new Particle*[NumberOfParticles];



	//activeParticleArray = new Particle*[NumberOfParticles];
//	indexArrayActiveParticles = new int[NumberOfParticles];
	//indexArrayRestrainedParticles = new int[NumberOfParticles];



	double px, py,pz;
	int i = 0;
	double e_kin = 0.0;
	double p_totx = 0.0;
	double p_toty = 0.0;
	double p_totz = 0.0;

	//	for (int j=1; j<=NumberOfParticles; j++)
	//	{


	int NumberOfParticlesPerAxe = (int)pow(NumberOfParticles, 1.0 / (double)d);


	NumberOfParticles = pow(NumberOfParticlesPerAxe, d);
	std::cout << "Number of particles on lattice " << NumberOfParticles << std::endl;

	int count = 0;

	for (int i = 1; i <= NumberOfParticlesPerAxe; i++){
		for (int j = 1; j <= NumberOfParticlesPerAxe; j++){
			for (int k = 1; k <= NumberOfParticlesPerAxe; k++){

#ifdef ONE_PARTICLE_ONLY

				px = 0.0;
				py = 0.01;
				pz = 0.0;

#endif

#ifndef ONE_PARTICLE_ONLY


				px = ComputeBoxMuller(0, 1);
				py = ComputeBoxMuller(0, 1);
				pz = ComputeBoxMuller(0, 1);


#endif
				double xPosition = i*BoxLength / (NumberOfParticlesPerAxe + 1);
				double yPosition = j*BoxLength / (NumberOfParticlesPerAxe + 1);
				double zPosition = k*BoxLength / (NumberOfParticlesPerAxe + 1);

				//	std::cout<<xPosition<<"\t "<<yPosition<<std::endl;


				//if(count == 0)	ptr=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);
				//else ptr=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);

				//pParticles[count]=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);//ptr;

				pParticles[count] = new Particle();
				pParticles[count]->setMomentum(px, py, pz);
				pParticles[count]->setPosition(xPosition, yPosition, zPosition);

				if (count < NumberOfDimerParticles) pParticles[count]->setParticleType(1);
				else pParticles[count]->setParticleType(0);

				pParticles[count]->setEps(epsr[count], epsf[count]);

				pParticles[count]->setParticleIndex(count);



				//	if(count == 0)	pParticles[count]->setMomentum(px,py);//=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);//ptr;
				//	else pParticles[count]=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);

				//	pParticles[count]=new Particle(px,py,xPosition,yPosition,1,epsr[count],epsf[count]);//ptr;

				p_totx = p_totx + pParticles[count]->getMomentumX();
				p_toty = p_toty + pParticles[count]->getMomentumY();
				p_totz = p_totz + pParticles[count]->getMomentumZ();



				count++;

			}
		}

	}

	NumberOfParticles = count;
	std::cout << "Number of Particles after box creation is " << NumberOfParticles << std::endl;

	//if (NumberOfDimerParticles> 1)	pParticles[1]->setPosition(pParticles[1]->getPositionX() + cut_off, pParticles[0]->getPositionY(),  pParticles[1]->getPositionZ());


	//if(pParticles[0]->getEpsr() !=0 && pParticles[0]->getEpsf() !=0 ) std::cout<< "Main particle with non-zero epsr/epsf !!!!"<<std::endl;


	p_totx = p_totx / NumberOfParticles;
	p_toty = p_toty / NumberOfParticles;
	p_totz = p_totz / NumberOfParticles;

	if (NumberOfParticles != 1){
		for (int i = 0; i < NumberOfParticles; i++)
		{
			pParticles[i]->setMomentum(pParticles[i]->getMomentumX() - p_totx, pParticles[i]->getMomentumY() - p_toty, pParticles[i]->getMomentumZ() - p_totz);

		}
	}



	double particlesEnergy = 0;


	for (int i = 0; i < NumberOfParticles; i++) {

		double m = pParticles[i]->getMass();
		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();
		double epsfi = pParticles[i]->getEpsf();
		double epsri = pParticles[i]->getEpsr();

		particlesEnergy += (px*px + py*py +pz*pz);//*(-Rho_p(px,py,epsri,epsfi,m))/m- (px*px+py*py)*((dRho_p(px,py,epsri,epsfi,m)*px+(dRho_p(py,px,epsri,epsfi,m)*py)))/(2.0*m);

	}

	//particlesEnergy=2.0*ComputeKineticEnergy();

	double T_0 = particlesEnergy / (NumberOfParticles*d*k);

	//double T_0=ComputeTemperature();

	double scale = 1;
	if (fabs(T_0) > 0.001) scale = sqrt(T / T_0);//sqrt((T+particlesEnergy/(NumberOfParticles*d*k))/T_0);

	if (NumberOfParticles != 1){
		for (int i = 0; i < NumberOfParticles; i++)
		{
			pParticles[i]->setMomentum(pParticles[i]->getMomentumX()*scale, pParticles[i]->getMomentumY()*scale, pParticles[i]->getMomentumZ()*scale);
			#ifndef EXP_PERTUBATION
			std::cout << "Particle " << i << "\t Epsr= " << pParticles[i]->getEpsr() << "\t Epsf= " << pParticles[i]->getEpsf() << std::endl;
			#endif
			#ifdef EXP_PERTUBATION
			std::cout << "Particle " << i << std::endl;
			#endif
		}
	}


#ifdef READ_INITIAL_CONDITION_FROM_FILE

	std::ofstream file("MATRIX");
	for (int i = 0; i < NumberOfParticles; i++){

		file <<	pParticles[i]->getMomentumX()<<" ";

	}

	file << std::endl;

	for (int i = 0; i < NumberOfParticles; i++){

		file << pParticles[i]->getMomentumY() << " ";

	}

	for (int i = 0; i < NumberOfParticles; i++){

		file <<	pParticles[i]->getMomentumZ()<<" ";

	}
	//while (1);

	int matrixSizeComlumn = NumberOfParticles;
	int matrixSizeRow = d;

	double** newMatrix = new double*[matrixSizeRow];
	for (int i = 0; i < matrixSizeRow; i++){

		newMatrix[i] = new double[matrixSizeComlumn];
	}

	newMatrix = readMatrixFromFile(matrixSizeRow, matrixSizeComlumn);

	/*for (int i = 0; i < matrixSizeRow; i++){
		for (int j = 0; j < matrixSizeComlumn; j++){

		std::cout << newMatrix[i][j]<<" ";
		}
		std::cout<<std::endl;
		}

		std::cout << newMatrix[0][1]<<std::endl;
		*/
	//	while(1);
#endif

	std::cout << "Initial temperature is " << ComputeTemperature() << "\n" << std::endl;


	std::cout << "Updating forces in initial condition..." << std::endl;

#ifndef NO_FORCES
#ifndef ONE_PARTICLE
	//UpdateForcesParticles();

#endif
#endif
	ActiveStatus();

#ifdef NVE

	std::cout << "\nNVE SIMULATION! " << std::endl;

	intialEnergyNVE = ComputeEnergy();
	//std::cout << "Initial momenta in NVE: "<<std::endl;

	//write initial conditions
	/*for (int i = 0; i < NumberOfParticles; i++) {


		std::cout << "px(" << i << ")=" << pParticles[i]->getMomentumX() << "\t py(" << i << ") = " << pParticles[i]->getMomentumY() << "\t pz(" << i << ") = " << pParticles[i]->getMomentumZ() << std::endl;

	}

	std::cout << "\n";

	for (int i = 0; i < NumberOfParticles; i++) {


		std::cout << "qx(" << i << ")=" << pParticles[i]->getPositionX() << "\t qy(" << i << ") = " << pParticles[i]->getPositionY() << "\t qz(" << i << ") = " << pParticles[i]->getPositionZ() << std::endl;

	}*/

#endif



	ValueObservable_IC = VarianceOne.getAverage();
	if (NumberOfDimerParticles > 0) std::cout << "Initial value for Dimer Distance is " << DimerDistance << std::endl;

	std::vector<Particle*> particleVector;
	for (int i = 0; i < NumberOfParticles; i++) particleVector.push_back(pParticles[i]);
	IAVector PBC(0, BoxLength, 0, BoxLength, 0, BoxLength);

	grid = new NeighborSearchGridPBC(particleVector, PBC, cut_off/*+0.01*/);
	grid->initializeNeighborLists();



	//grid->print();

	//AddSamples();
	std::cout << "Update Forces in intial conditions \n" << std::endl;

	UpdateForcesParticlesAdaptivelyNewton();

	potentialEnergy = ComputePotentialEnergy();

	std::cout << "Number of particles= " << NumberOfParticles << " from which number of dimer= " << NumberOfDimerParticles << std::endl;
	std::cout << "Initial potential energy " << potentialEnergy << std::endl;

};



//***************** Metropolization


double Simulation::MetropolisHamiltonianARPS(double E_old){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "MetropolisHamiltonianARPS " << std::endl;
#endif

	int rejectionIndicator = 0;

#ifdef METROPOLIS_FDR_ADAPTIVE
	double E_new = potentialEnergy+ComputeKineticEnergy();
#endif
#ifndef METROPOLIS_FDR_ADAPTIVE
	double E_new =  ComputeEnergy();
#endif

	double alpha = rMetropolis->randDouble2();

	double alphaDt = beta*(E_new - E_old); // exp( -beta*(new_energy-X.old_energy(k)) );
	double acceptanceRate = exp(-alphaDt);


	if (acceptanceRate < alpha){

		rejectionIndicator = 1;


		for (int i = 0; i<NumberOfParticles; i++) {

			pParticles[i]->setMomentum(-(pParticles[i]->getMomentumXOld()), -(pParticles[i]->getMomentumYOld()), -(pParticles[i]->getMomentumZOld()));
			pParticles[i]->setPosition(pParticles[i]->getPositionXOld(), pParticles[i]->getPositionYOld(), (pParticles[i]->getPositionZOld()));
			pParticles[i]->setFOld();
			E_new = E_old;

#ifdef METROPOLIS_FDR_ADAPTIVE
			potentialEnergy = potentialEnergyOld;
#endif

			grid->onPositionChanged(pParticles[i]);
		}


	}
#ifdef METROPOLIS_FDR_ADAPTIVE
	else{

		grid->updateNeighborLists();

	}
#endif

	RejectionRate.addSample(rejectionIndicator);


#ifdef WRITE_INSTANT_REJECTION_RATE

	instantRejectionHam = rejectionIndicator;

#endif

	return E_new;
}





double Simulation::TproposalARPS(double pOldx, double pOldy, double pOldz, double pNewx, double pNewy, double pNewz, double m, double epsrn, double epsfn, double t){


	//double Tx = (pNewx - pOldx + gamma*t*dHdp(pOldx,pOldy, m, epsrn, epsfn))*(pNewx - pOldx + gamma*t*dHdp(pOldx,pOldy, m, epsrn, epsfn)) / (2 * gamma*t);//

	double h = sqrt(2 * gamma*t);

	double Tx = pNewx - pOldx + gamma*t*dHdp(pOldx, pOldy, pOldz, m, epsrn, epsfn);
	double Ty = pNewy - pOldy + gamma*t*dHdp(pOldy, pOldx, pOldz, m, epsrn, epsfn);
	double Tz = pNewz - pOldz + gamma*t*dHdp(pOldz, pOldy, pOldx, m, epsrn, epsfn);


	double Tabs = (Tx*Tx+Ty*Ty) / (2 * gamma*t);//

	//double T=((pNew - pOld+gamma*dt*pOld/m)*(pNew-pOld+gamma*dt*pOld/m))/(2*gamma*dt);

	//std::cout << Tboth << std::endl;

	return exp(-(1.0 / kT)*Tabs / 2.0);
}
//
//void Simulation::TproposalAlternative(double Told, double t){
//
//
//	double Tnew = ComputeKineticEnergy();
//
//	int rejectionIndicator = 0;
//
//	double Tprop_tmp1 = 0.0;
//	double Tprop_tmp2 = 0.0;
//
//	for (int i = 0; i < NumberOfParticles; i++){
//
//		//	numberOfMetropolisSteps++;
//
//		double px = pParticles[i]->getMomentumX();
//		double py = pParticles[i]->getMomentumY();
//
//		double pOldx = pParticles[i]->getMomentumXOld();
//		double pOldy = pParticles[i]->getMomentumYOld();
//
//		double m = pParticles[i]->getMass();
//		double epsri = pParticles[i]->getEpsr();
//		double epsfi = pParticles[i]->getEpsf();
//
//		//
//
//
//		double Tx = px - pOldx + gamma*t*dHdp(pOldx, pOldy, m, epsri, epsfi);
//		double Ty = py - pOldy + gamma*t*dHdp(pOldy, pOldx, m, epsri, epsfi);
//
//
//		Tprop_tmp1 = Tprop_tmp1 + (Tx*Tx + Ty*Ty) / (2 * gamma*t);//
//
//
//		double Tx2 = pOldx - px + gamma*t*dHdp(px, py, m, epsri, epsfi);
//		double Ty2 = pOldy - py + gamma*t*dHdp(py, px, m, epsri, epsfi);
//
//
//		Tprop_tmp2 = Tprop_tmp2 + (Tx2*Tx2 + Ty2*Ty2) / (2 * gamma*t);//
//
//		//double T=((pNew - pOld+gamma*dt*pOld/m)*(pNew-pOld+gamma*dt*pOld/m))/(2*gamma*dt);
//
//		//std::cout << Tboth << std::endl;
//
//
//
//
//	}
//
//	double Tproposal1 = exp(-(1.0 / kT)*Tprop_tmp1 / 2.0);
//	double Tproposal2 = exp(-(1.0 / kT)*Tprop_tmp2 / 2.0);
//
//	//SBRandom r((unsigned long)SBCTimeClock());
//	//double alpha = r.randDouble2(); // random number in [0,1)
//	double alpha = rMetropolis->randDouble2(); // random number in [0,1)
//
//	//		double T_pToPold = TproposalARPS(px, py, pOldx, pOldy, m, epsri, epsfi, t);
//	//		double T_pOld_To_p = TproposalARPS(pOldx, pOldy, px, py, m, epsri, epsfi, t);
//
//	//	std::cout << Tnew - Told << std::endl;
//
//	double acceptanceProba = 0;
//	//	if (T_pOld_To_p > 0)
//	acceptanceProba = exp(-(1.0 / kT)*(Tnew - Told))*Tproposal2 / Tproposal1;
//
//	//	std::cout<< Tproposal(pOld,p,m,epsri,epsfi)/Tproposal(p,pOld,m,epsri,epsfi)<<std::endl;
//
//
//
//	if (alpha > acceptanceProba){
//
//		//rejection_rate++;
//		rejectionIndicator = 1;
//
//		for (int i = 0; i < NumberOfParticles; i++)
//			pParticles[i]->setMomentum(pParticles[i]->getMomentumXOld(), pParticles[i]->getMomentumYOld());
//
//
//	}
//
//
//
//
//	RejectionRate.addSample(rejectionIndicator);
//
//
//};

void Simulation::MetropolisStochARPS(double Told, double t){

	double Tnew = ComputeKineticEnergy();

	int rejectionIndicator = 0;

	for (int i = 0; i<NumberOfParticles; i++){

		//	numberOfMetropolisSteps++;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double pOldx = pParticles[i]->getMomentumXOld();
		double pOldy = pParticles[i]->getMomentumYOld();
		double pOldz = pParticles[i]->getMomentumZOld();

		double m = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		//

		//SBRandom r((unsigned long)SBCTimeClock());
		//double alpha = r.randDouble2(); // random number in [0,1)
		double alpha = rMetropolis->randDouble2(); // random number in [0,1)

		double T_pToPold=TproposalARPS(px, py, pz,pOldx, pOldy,pOldz, m, epsri, epsfi, t);
		double T_pOld_To_p = TproposalARPS(pOldx, pOldy, pOldz,px, py, pz, m, epsri, epsfi, t);

		std::cout << Tnew - Told << std::endl;

		double acceptanceProba = 0;
		if (T_pOld_To_p > 0)
			double acceptanceProba = exp(-(1.0 / kT)*(Tnew - Told))*T_pToPold / T_pOld_To_p;

		//	std::cout<< Tproposal(pOld,p,m,epsri,epsfi)/Tproposal(p,pOld,m,epsri,epsfi)<<std::endl;



		if (alpha > acceptanceProba){

			//rejection_rate++;

			//for(int i=0; i<NUMBER_OF_PARTICLES; i++)
			pParticles[i]->setMomentum(pParticles[i]->getMomentumXOld(), pParticles[i]->getMomentumYOld(), pParticles[i]->getMomentumZOld());

			rejectionIndicator = 1;
		}

	}

	RejectionRate.addSample(rejectionIndicator);

};


void Simulation::MetropolisEulerMaruyama(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "MetropolisEulerMaruyama "<<t << std::endl;
#endif



	//metropolize particle par particle

	double rejectionIndicator = 0;
	SBRandom r((unsigned long)SBCTimeClock());

	for (int i = 0; i<NumberOfParticles; i++){

		double T_new = pParticles[i]->getEkin();
		double T_old = pParticles[i]->getEkinOld();

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double pxOld = pParticles[i]->getMomentumXOld();
		double pyOld = pParticles[i]->getMomentumYOld();
		double pzOld = pParticles[i]->getMomentumZOld();

		double mass = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();


		double alpha = r.randDouble2(); // random number in [0,1)

				double alphaDtX = (T_old - T_new) / kT + ( TproposalEulerMaruyama(px, py, pz, pxOld, pyOld, pzOld, mass, epsri, epsfi, t) - TproposalEulerMaruyama(pxOld, pyOld, pzOld, px, py, pz, mass, epsri, epsfi, t)) / (2 * kT);
				double acceptanceProbaX = exp(-alphaDtX);

				if (alpha > acceptanceProbaX){

					rejectionIndicator = 1;

					px = pxOld;

					//for(int i=0; i<NUMBER_OF_PARTICLES; i++)

				}

				double alphaDtY = (T_old - T_new) / kT + (TproposalEulerMaruyama(py, px, pz, pyOld, pxOld, pzOld, mass, epsri, epsfi, t) - TproposalEulerMaruyama(pyOld, pxOld, pzOld, py, px, pz, mass, epsri, epsfi, t)) / (2 * kT);
				double acceptanceProbaY = exp(-alphaDtY);

				alpha = r.randDouble2();

				if (alpha > acceptanceProbaY){

					rejectionIndicator ++;

					py = pyOld;



				}


				double alphaDtZ = (T_old - T_new) / kT + (TproposalEulerMaruyama(pz, px, py, pyOld, pzOld, pyOld, mass, epsri, epsfi, t) - TproposalEulerMaruyama(pzOld, pxOld, pyOld, pz, px, py, mass, epsri, epsfi, t)) / (2 * kT);
				double acceptanceProbaZ = exp(-alphaDtZ);

				alpha = r.randDouble2();

				if (alpha > acceptanceProbaZ){

					rejectionIndicator++;

					pz = pzOld;
				}

				pParticles[i]->setMomentum(px, py, pz);
			}

	RejectionRate.addSample(rejectionIndicator );
};








// ***************** end metropolization

double Simulation::TproposalEulerMaruyama(double pxNew, double pyNew, double pzNew, double pxOld, double pyOld, double pzOld, double m, double epsrn, double epsfn, double t){

	double dHpOld = dHdp(pxOld, pyOld, pzOld,  m, epsrn, epsfn);
	double Gtmp = (pxNew - pxOld + gamma*t*dHpOld);

	return Gtmp*Gtmp / (2 * gamma*t);//


};


double Simulation::MidPoint(double x, double y){

	return (x + y) / 2.0;

}


void Simulation::PerformTimeStep(){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "****start time step****" << std::endl;
#endif

#ifdef FD_ONLY

	//UpdateMomentaStochOrder2_CN(dt);
	//UpdateMomentaStochModifiedEqCN(dt);

	UpdateMomentaStochSecondOrder(dt);

	//UpdateMomentaStochLangevin(dt);

	ActiveStatus();

	//UpdateMomentaStochExact(0.5*dt);
	//UpdateMomentaARPSperturbation(dt);
	//UpdateMomentaStochExact(0.5*dt);

	//UpdateMomentaStochOrder1(dt);

#endif

#ifdef METROPOLIZATION


#ifdef METROPOLIS_FDR

	SaveCurrentMomentaAsOld();
	SaveCurrentPositionsAsOld();

	double oldKinetEnergy = ComputeKineticEnergy();
	double oldPotentialEnergy = potentialEnergy;

	UpdateMomentaParticles(0.5*dt);
	ActiveStatus();

	//potentialEnergy = ComputePotentialEnergyAdaptivelySubtract(potentialEnergy);
	//UpdateForcesParticlesAdaptivelyARPS_SUBTRACT();

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	UpdateForcesParticlesNewton();

	//grid->updateNeighborLists();
	//UpdateForcesParticlesAdaptivelyNewton();

	UpdateMomentaParticles(0.5*dt);

	//potentialEnergy = ComputePotentialEnergy();

	double energy = MetropolisHamiltonianARPS(oldKinetEnergy + oldPotentialEnergy);

	PeriodicBoundaryCondition();
	grid->updateNeighborLists();

	//ActiveStatus();

	SampleG();
	SaveCurrentMomentaAsOld();


	//double* oldKinEn = ComputeEnergyFDR();

	//UpdateMomentaVerlet_FDR(h);

	//UpdateMomentaStochOrder2_CN(dt);
	UpdateMomentaVerlet_FDR_MetropolizePerParticle(dt);


#endif



#ifdef METROPOLIS_FDR_ADAPTIVE

	SaveCurrentMomentaAsOld();
	SaveCurrentPositionsAsOld();

	double oldKinetEnergy = ComputeKineticEnergy(); //linear compl.
	double oldPotentialEnergy = potentialEnergy;

	UpdateMomentaParticles(0.5*dt);
	ActiveStatus();

	//potentialEnergy = ComputePotentialEnergyAdaptivelySubtract(potentialEnergy);
	//UpdateForcesParticlesAdaptivelyARPS_SUBTRACT();

	doAddSubst = 0;
	double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;

	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();


	if (doAddSubst == 1){


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();
		potentialEnergy = ComputePotentialEnergy();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}



	//grid->updateNeighborLists();
	//UpdateForcesParticlesAdaptivelyNewton();

	UpdateMomentaParticles(0.5*dt);

	//potentialEnergy = ComputePotentialEnergy();

	double energy = MetropolisHamiltonianARPS(oldKinetEnergy + oldPotentialEnergy);

	//PeriodicBoundaryCondition();
	//grid->updateNeighborLists();

	//ActiveStatus();

	SampleG();
	SaveCurrentMomentaAsOld();


	//double* oldKinEn = ComputeEnergyFDR();

	//UpdateMomentaVerlet_FDR(h);

	//UpdateMomentaStochOrder2_CN(dt);
	UpdateMomentaVerlet_FDR_MetropolizePerParticle(dt);

//	if ((current_n != 0 && (current_n % 10000 == 0)) )
//	std::cout << "Check of computation of potential energy: " << potentialEnergy - ComputePotentialEnergy() << std::endl;

#endif



#ifdef METROPOLIS_LANGEVIN


	SaveCurrentMomentaAsOld();
	SaveCurrentPositionsAsOld();

	double oldKinetEnergy = ComputeKineticEnergy();
	double oldPotentialEnergy = ComputePotentialEnergy();//potentialEnergy;

	UpdateMomentaParticles(0.5*dt);
	ActiveStatus();

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	UpdateForcesParticlesNewton();

	//grid->updateNeighborLists();

	//UpdateForcesParticlesAdaptivelyNewton();
	UpdateMomentaParticles(0.5*dt);

	double energy = MetropolisHamiltonianARPS(oldKinetEnergy + oldPotentialEnergy);

	//PeriodicBoundaryCondition();
	//grid->updateNeighborLists();


	//UpdateMomentaStochOrder2_CN(dt);

	MetropolisStochLangevin(dt);
	//MetropolisStochLangevinSeparableSpace(dt);
	//UpdateMomentaVerlet_FDR_MetropolizePerParticle(dt);

#endif

#endif

#ifdef ADAPTIV_UPDATE_FORCES

#ifdef Z_BAOAB

	UpdateMomentaParticles(0.5*dt);
	UpdateMomentaARPSperturbation(0.5* dt);

//	UpdatePositionsParticles(0.5*dt);
//	PeriodicBoundaryCondition();



	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;

	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}


	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}


	UpdateMomentaARPSperturbation(0.5* dt);
	UpdateMomentaParticles(0.5*dt);


	UpdateMomentaStochExact(dt);

#endif

#ifdef BAOAB

	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	UpdatePositionsParticles(0.5*dt);
	PeriodicBoundaryCondition();

	UpdateMomentaStochOrder2_CN(dt);

	ActiveStatus();


	doAddSubst = 0;
	double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;

	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}


	UpdatePositionsParticles(0.5*dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}

	UpdateMomentaParticles(0.5*dt);


#endif


#ifdef OBABO

#ifndef NVE
//	UpdateMomentaStochOrder2_CN(0.5*dt);
//	UpdateMomentaStochModifiedEqCN(0.5*dt);

	UpdateMomentaStochSecondOrder(0.5*dt);

#endif
	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;

	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();


	if (doAddSubst == 1){


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}


	UpdateMomentaParticles(0.5*dt);

#ifndef NVE

//	UpdateMomentaStochOrder2_CN(0.5*dt);
//	UpdateMomentaStochModifiedEqCN(0.5*dt);

	//UpdateMomentaStochOrder2_CN(dt);
	//UpdateMomentaStochModifiedEqCN(dt);

	UpdateMomentaStochSecondOrder(0.5*dt);

#endif

#endif


#ifdef OBABO_noNBL

	//UpdateMomentaStochOrder2_CN(0.5*dt);
	//UpdateMomentaStochModifiedEqCN(0.5*dt);



	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif
		UpdateForcesParticlesNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif

		UpdateMomentaParticles(0.5*dt);

	//UpdateMomentaStochOrder2_CN(0.5*dt);
	//UpdateMomentaStochModifiedEqCN(0.5*dt);

	UpdateMomentaStochOrder2_CN(dt);
	UpdateMomentaStochModifiedEqCN(dt);

#endif

#ifdef ABOBA



		//UpdatePositionsParticles(0.5*dt);
		//PeriodicBoundaryCondition();

	UpdatePositionsParticles(0.5*dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();


	if (doAddSubst == 1){


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}


		UpdateMomentaParticles(0.5*dt);

		UpdateMomentaStochOrder2_CN(dt);

		UpdateMomentaParticles(0.5*dt);

		ActiveStatus();

		doAddSubst = 0;
		double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;

		if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

			doAddSubst = 1;

			//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
			UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
		}



		UpdatePositionsParticles(0.5*dt);
		PeriodicBoundaryCondition();

		grid->updateNeighborLists();






#endif

#ifdef BABM


		UpdateMomentaStochOrder1(0.5*dt);
		UpdateMomentaStochModifiedEq(0.5*dt);


	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate= CONDITION_ADAPTIVE_FORCES_UPDATE;


	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){

		//UpdateForcesParticlesAdaptivelyARPS_ADD();


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}

	UpdateMomentaParticles(0.5*dt);

	UpdateMomentaStochOrder1(0.5*dt);
	UpdateMomentaStochModifiedEq(0.5*dt);


	/*UpdateMomentaARPSperturbation(0.5*dt);
	UpdateMomentaStochExact(dt);
	UpdateMomentaARPSperturbation(0.5*dt);
*/
#endif



#ifdef MBABM


	//UpdateMomentaFD_BDV(0.5*dt);

#ifndef NVE
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate= CONDITION_ADAPTIVE_FORCES_UPDATE;


	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){

		//UpdateForcesParticlesAdaptivelyARPS_ADD();


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}

	UpdateMomentaParticles(0.5*dt);

	//UpdateMomentaFD_BDV(0.5*dt);

#ifndef NVE
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

	//SaveCurrentMomentaAsOld();
	//	UpdateMomentaVerlet_FDR_MetropolizePerParticle(dt);

	/*UpdateMomentaARPSperturbation(0.5*dt);
	UpdateMomentaStochExact(dt);
	UpdateMomentaARPSperturbation(0.5*dt);
	*/
#endif

#ifdef BABE

#ifndef NVE
	UpdateMomentaStochOrder2_CN(0.5*dt);
	UpdateMomentaStochModifiedEqCN(0.5* dt);
#endif

#ifndef ONLY_FLUCTUATION_DISSIPATION_PART
	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate= CONDITION_ADAPTIVE_FORCES_UPDATE;


	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){

		//UpdateForcesParticlesAdaptivelyARPS_ADD();


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}

	UpdateMomentaParticles(0.5*dt);

#endif

#ifndef NVE
	//UpdateMomentaStochExactARPS(dt);
	UpdateMomentaStochModifiedEqCN(0.5* dt);
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

#endif

	//--------------------------------------------------


#ifdef BAB_MOMENTA_LANGEVIN


#ifndef NVE
	UpdateMomentaStochLangevin(0.5*dt);
	//UpdateMomentaStochExact(dt);
#endif

#ifndef ONLY_FLUCTUATION_DISSIPATION_PART
	UpdateMomentaParticles(0.5*dt);

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate = CONDITION_ADAPTIVE_FORCES_UPDATE;


	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();

	if (doAddSubst == 1){

		//UpdateForcesParticlesAdaptivelyARPS_ADD();


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}

	UpdateMomentaParticles(0.5*dt);

#endif

#ifndef NVE
	UpdateMomentaStochLangevin(0.5*dt);
#endif

#endif

	//-------------------------------------------------

#ifdef FD_EXACT

	UpdateMomentaStochExactARPS(dt);

	//UpdateMomentaStochOrder2_CN(dt);

#endif


#ifndef BAOAB_SCHEMES

#ifdef FIRST_ORDER_SPLITTING

	UpdateMomentaParticles(dt);

#endif

#ifdef SECOND_ORDER_SPLITTING

#ifdef SECOND_ORDER_STRANG

#ifndef NVE
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

#endif

	UpdateMomentaParticles(0.5*dt);

#endif

	ActiveStatus();

	doAddSubst = 0;
	double ConditionUpdate= CONDITION_ADAPTIVE_FORCES_UPDATE;


	if ((double)Particle::howManyRestrained / (double)NumberOfParticles > ConditionUpdate)	{// 0.75)	{

		doAddSubst = 1;

		//UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
		UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT();
	}

	UpdatePositionsParticles(dt);
	PeriodicBoundaryCondition();

	grid->updateNeighborLists();



	if (doAddSubst == 1){

		//UpdateForcesParticlesAdaptivelyARPS_ADD();


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewtonARPS_ADD();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}
	else{


#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep0 = SBCTimeClock();
#endif

		UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
		tPerTimeStep1 = SBCTimeClock();
#endif
	}


#ifdef SECOND_ORDER_SPLITTING

	UpdateMomentaParticles(0.5*dt);

#endif

#ifndef NVE

#ifdef FIRST_ORDER_SPLITTING

	UpdateMomentaStochOrder1(dt);

#endif

#ifdef SECOND_ORDER_SPLITTING

#ifdef SECOND_ORDER_STRANG
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

#ifdef SECOND_ORDER_FD_MIDPOINT
	UpdateMomentaStochOrder2(dt);
#endif

#ifdef SECOND_ORDER_FD_CN
	UpdateMomentaStochOrder2_CN(dt);
#endif

#endif


#endif

#endif


#endif


	//*************************************************************



#ifdef STANDARD_UPDATE_FORCES

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "****start time step****" << std::endl;
#endif

#ifdef FIRST_ORDER_SPLITTING

	UpdateMomentaParticles(dt);

#endif

#ifdef SECOND_ORDER_SPLITTING

#ifdef SECOND_ORDER_STRANG

#ifndef NVE
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif
#endif


	UpdateMomentaParticles(0.5*dt);

#endif

	//ActiveStatus();

	UpdatePositionsParticlesStandard(dt);
	PeriodicBoundaryCondition();

	//tPerTimeStep0 = SBCTimeClock();

	grid->updateNeighborLists();

	//tPerTimeStep1 = SBCTimeClock();

	// sveta
	//UpdateForcesParticlesAdaptively();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
	tPerTimeStep0 = SBCTimeClock();
#endif
	UpdateForcesParticlesAdaptivelyNewton();

#ifdef MEASURE_TIME_PER_FORCE_UPDATE
	tPerTimeStep1 = SBCTimeClock();
#endif
	//UpdateForcesParticles();


#ifdef SECOND_ORDER_SPLITTING

	UpdateMomentaParticles(0.5*dt);

#endif


#ifndef NVE

#ifdef FIRST_ORDER_SPLITTING

	UpdateMomentaStochOrder1(dt);

#endif

#ifdef SECOND_ORDER_SPLITTING

#ifdef SECOND_ORDER_STRANG
	UpdateMomentaStochOrder2_CN(0.5*dt);
#endif

#ifdef SECOND_ORDER_FD_MIDPOINT
	UpdateMomentaStochOrder2(dt);
#endif

#ifdef SECOND_ORDER_FD_CN
	UpdateMomentaStochOrder2_CN(dt);
#endif

#endif

#endif
#endif


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "****end time step****" << std::endl;
#endif

};


//
//void Simulation::PerformTimeStep(){
//
//
//#ifdef CHECK_SPLITTING_SCHEME
//	if (current_n == 0) std::cout << "****start time step****" << std::endl;
//#endif
//
//#ifdef ADAPTIV_UPDATE_FORCES
//
//#ifdef FIRST_ORDER_SPLITTING
//
//	UpdateMomentaParticles(dt);
//
//#endif
//
//#ifdef SECOND_ORDER_SPLITTING
//
//#ifdef SECOND_ORDER_STRANG
//
//#ifndef NVE
//	UpdateMomentaStochOrder2(0.5*dt);
//#endif
//
//#endif
//
//	UpdateMomentaParticles(0.5*dt);
//
//#endif
//
//	ActiveStatus();
//
//	//tPerTimeStep0 = SBCTimeClock();
//
//	UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT();
//
//	//tPerTimeStep1 = SBCTimeClock();
//
//	UpdatePositionsParticles(dt);
//	PeriodicBoundaryCondition();
//
//	//tPerTimeStep0 = SBCTimeClock();
//
//	grid->updateNeighborLists();
//
//	//tPerTimeStep1 = SBCTimeClock();
//	/*
//	if (current_n == 1000){
//		grid->print();
//
//	}
//*/
//	//no condition like * needed with neighbor lists!!!!
//	//*if ((double)Particle::howManyRestrained / (double)NumberOfParticles >= 0)	{
//	//*********************************
//
//	//UpdateForcesParticlesAdaptively();
//
//
//	//UpdateForcesParticlesAdaptivelyNewton();
//
//
//	//UpdateForcesParticles();
//
//	//tPerTimeStep0 = SBCTimeClock();
//
//	UpdateForcesParticlesAdaptivelyARPS_ADD();
//
//	//tPerTimeStep1 = SBCTimeClock();
//
//#ifdef SECOND_ORDER_SPLITTING
//
//	UpdateMomentaParticles(0.5*dt);
//
//#endif
//
//#ifndef NVE
//
//#ifdef FIRST_ORDER_SPLITTING
//
//	UpdateMomentaStochOrder1(dt);
//
//#endif
//
//#ifdef SECOND_ORDER_SPLITTING
//
//#ifdef SECOND_ORDER_STRANG
//	UpdateMomentaStochOrder2(0.5*dt);
//#endif
//
//#ifdef SECOND_ORDER_FD_MIDPOINT
//	UpdateMomentaStochOrder2(dt);
//#endif
//
//#ifdef SECOND_ORDER_FD_CN
//	UpdateMomentaStochOrder2_CN(dt);
//#endif
//
//#endif
//
//
//#endif
//
//#endif
//
//
////*************************************************************
//
//
//
//#ifdef STANDARD_UPDATE_FORCES
//
//#ifdef CHECK_SPLITTING_SCHEME
//	if (current_n == 0) std::cout << "****start time step****" << std::endl;
//#endif
//
//#ifdef FIRST_ORDER_SPLITTING
//
//	UpdateMomentaParticles(dt);
//
//#endif
//
//#ifdef SECOND_ORDER_SPLITTING
//
//#ifdef SECOND_ORDER_STRANG
//
//#ifndef NVE
//	UpdateMomentaStochOrder2(0.5*dt);
//#endif
//#endif
//
//
//	UpdateMomentaParticles(0.5*dt);
//
//#endif
//
//	//ActiveStatus();
//
//	UpdatePositionsParticlesStandard(dt);
//	PeriodicBoundaryCondition();
//
//	//tPerTimeStep0 = SBCTimeClock();
//
//	grid->updateNeighborLists();
//
//	//tPerTimeStep1 = SBCTimeClock();
//
//	// sveta
//	//UpdateForcesParticlesAdaptively();
//
//	//tPerTimeStep0 = SBCTimeClock();
//
//	UpdateForcesParticlesAdaptivelyNewton();
//
//	//tPerTimeStep1 = SBCTimeClock();
//
//	//UpdateForcesParticles();
//
//
//#ifdef SECOND_ORDER_SPLITTING
//
//	UpdateMomentaParticles(0.5*dt);
//
//#endif
//
//
//#ifndef NVE
//
//#ifdef FIRST_ORDER_SPLITTING
//
//	UpdateMomentaStochOrder1(dt);
//
//#endif
//
//#ifdef SECOND_ORDER_SPLITTING
//
//#ifdef SECOND_ORDER_STRANG
//	UpdateMomentaStochOrder2(0.5*dt);
//#endif
//
//#ifdef SECOND_ORDER_FD_MIDPOINT
//	UpdateMomentaStochOrder2(dt);
//#endif
//
//#ifdef SECOND_ORDER_FD_CN
//	UpdateMomentaStochOrder2_CN(dt);
//#endif
//
//#endif
//
//#endif
//#endif
//
//
//#ifdef CHECK_SPLITTING_SCHEME
//	if (current_n == 0) std::cout << "****end time step****" << std::endl;
//#endif
//
//};

void Simulation::PerformTrajectory(){




#ifdef WRITE_DIMER_DISTANCE_MATLAB


	std::ofstream FileDimerDistance("dimer_distance.m");
	FileDimerDistance << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n";
	FileDimerDistance << "sigma_LJ=" << sigma_LJ << ";\n epsilon_LJ=" << epsilon_LJ << ";\n hDimer=" << hDimer << ";\n wDimer=" << wDimer << ";\n cut_off=" << cut_off << ";\n";

	FileDimerDistance << "d=[";
#endif

#ifdef WRITE_POSITIONS_MATLAB

	std::ofstream FilePositions("positions_particles.m");


#ifdef ONE_PARTICLE_ONLY

	FilePositions << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "\n q=[" << std::endl;
#endif

#endif

#ifdef WRITE_INST_ENERGY

	std::ofstream FileEnergy("energy.m");

	FileEnergy << "E=[" << std::endl;

#endif

#ifdef WRITE_MOMENTA_MATLAB

	std::ofstream FileMomenta("momenta_particles.m");



#ifdef ONE_PARTICLE_ONLY

	FileMomenta << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "\n p=[" << std::endl;
#endif
#endif

#ifdef WRITE_MOMENTA_MATLAB_HISTOGRAM
	std::ofstream FileMomentaHist("momenta_particles_histogram.m");
	FileMomentaHist << "p=[" << std::endl;
#endif


#ifdef WRITE_P_SQUARE_2D

	std::ofstream FilePSquare("pSquare.m");
	FilePSquare <<  "P2=[\n";
#endif


#ifdef WRITE_ERROR_LANGEVIN

	std::ofstream FileError("errorAverage.m");
	std::ofstream FileErrorInst("errorInst.m");
	FileError << "Eaver=[\n";
	FileErrorInst << "Einst=[\n";

#endif




#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS

	double howManyTimesAverages = HOW_MANY_PHYSICAL_TIMES;
	double writingIntervalAverages = NumberOfTimeSteps / howManyTimesAverages;
	std::cout << "\n Writing interval for averages =" << writingIntervalAverages << std::endl;

#endif


#ifdef NVE

	PerformTimeStep();
	intialEnergyNVE = ComputeEnergy();// ComputeKineticEnergy();
	std::cout << "Initial energy in NVE is " << intialEnergyNVE << std::endl;


#endif

	current_n = 0;
	measureTime = 0.0;

	while (current_n < NumberOfEquilibrationSteps){

		if (errorIndicator == true) current_n = NumberOfEquilibrationSteps;


		//tPerTimeStep0 = SBCTimeClock();

		PerformTimeStep();

		//tPerTimeStep1 = SBCTimeClock();

		current_n++;

	}

#ifdef NVE

	std::cout << "computing the initial energy after equilibration\n" << std::endl;

	intialEnergyNVE = ComputeEnergy();

#endif

	current_n = 0;// 3 * MaxLengthOfCurrentN;// 2.1*pow(10.0, 9);


	BA_count = 0;

	while ((((double)current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN)) < NumberOfTimeSteps){


		if (errorIndicator == true)  NumberOfTimeSteps = 0.0;

		Tick();

#ifdef MEASURE_TIME_PER_TIMESTEP
		tPerTimeStep3 = SBCTimeClock();
#endif

		PerformTimeStep();

#ifdef MEASURE_TIME_PER_TIMESTEP
		tPerTimeStep4 = SBCTimeClock();
#endif

		DisplayTimeLeft();



#ifdef WRITE_DIMER_DISTANCE_MATLAB

		//if ((int)current_n % 10 == 0)
		FileDimerDistance << AutocorrelationAverage.getAverage() << "\n";//FileDimerDistance  << DimerDistance  << "\n";//

#endif

#ifdef WRITE_POSITIONS_MATLAB

#ifndef ONE_PARTICLE_ONLY
		int writeEverryN=1;

		if ((current_n % writeEverryN == 0) && current_n < MaxLengthOfCurrentN){

			//	FilePositions << "px(" << current_n + 1 << ")=" << pParticles[0]->getMomentumX() << ";" << std::endl;

			//if (NumberOfParticles == 1)			FilePositions << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "q(" << (current_n + 1)/writeEverryN << ",:)=[\n";
			//else
			FilePositions << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "q(:,:,:," << (current_n + 1) / writeEverryN << ")=[\n";


			for (int i = 0; i < NumberOfParticles; i++)
			{

				double qx = pParticles[i]->getPositionX();
				double qy = pParticles[i]->getPositionY();
				double qz = pParticles[i]->getPositionZ();
				unsigned int isRestrained = pParticles[i]->getActiveStatus();

				FilePositions << qx << "\t" << qy <<"\t" << qz <<  "\t"<<isRestrained<< "\n";//"\n";//

			}
			FilePositions << "];\n" << "\n";

		}

#endif

#ifdef ONE_PARTICLE_ONLY



		if (current_n < MaxLengthOfCurrentN){

			for (int i = 0; i < NumberOfParticles; i++)
			{

				double qx = pParticles[i]->getPositionX();
				double qy = pParticles[i]->getPositionY();
				double qz = pParticles[i]->getPositionY();


				FilePositions << qx << "\t" << qy "\t" << qz  << "\n";//"\n";//

			}


		}

#endif

#endif


#ifdef WRITE_MOMENTA_MATLAB

#ifndef WRITE_POSITIONS_MATLAB
		int writeEverryN=1;

#endif
#ifndef ONE_PARTICLE_ONLY
		if ((current_n % writeEverryN == 0)&& current_n < MaxLengthOfCurrentN){

			FileMomenta << "epsR=" << pParticles[NumberOfParticles - 1]->getEpsr() << "; epsF=" << pParticles[NumberOfParticles - 1]->getEpsf() << ";" << std::endl;

			//	FilePositions << "px(" << current_n + 1 << ")=" << pParticles[0]->getMomentumX() << ";" << std::endl;

			FileMomenta << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "p(:,:," << (current_n + 1) / writeEverryN << ")=[\n";


			for (int i = 0; i < NumberOfParticles; i++)
			{

				double px = pParticles[i]->getMomentumX();
				double py = pParticles[i]->getMomentumY();
				double pz = pParticles[i]->getMomentumZ();

				FileMomenta << px << "\t" << py << "\t" << pz <<  "\n";//"\n";//

			}
			FileMomenta << "];\n" << "\n";

		}

#endif

#ifdef ONE_PARTICLE_ONLY

		for (int i = 0; i < NumberOfParticles; i++)
		{

			double px = pParticles[i]->getMomentumX();
			double py = pParticles[i]->getPositionY();
			double pz = pParticles[i]->getPositionZ();

#ifdef WRITE_P_SQUARE_2D

			double p2 = px*px + py * py+ pz * pz;
			FilePSquare << p2 << "\n";
#endif

			FileMomenta << px << "\t" << py"\t" << pz << "\n";//"\n";//

		}

#endif

#endif

#ifdef WRITE_MOMENTA_MATLAB_HISTOGRAM

		int writeEverryNhist = 1;

		if ((current_n % writeEverryNhist == 0) && current_n < MaxLengthOfCurrentN){


			for (int i = 0; i < NumberOfParticles; i++)
			{

				double px = pParticles[i]->getMomentumX();
				double py = pParticles[i]->getMomentumY();

				FileMomentaHist << px << "\t" << py"\t" << pz << "\t" ;//"\n";//

			}


		}





#endif

#ifdef WRITE_INST_ENERGY


		FileEnergy << ComputeEnergy() << "\n";

#endif


#ifdef WRITE_ERROR_LANGEVIN

#ifndef NVE
		std::cout << "NVE turned off! Write_error_langevin working only with NVE\n"<<std::endl;
#endif
		//double ekin = ComputeEnergy();// ComputeKineticEnergy();

		//if ((int)current_n % 10 == 0)
		//FileError << 100 * abs((ekin - oldEnergy) / ekin) << "\n";


		double eCurrent = ComputeEnergy();// ComputeKineticEnergy();
		double errorEnergyPercentage = 100 * abs((eCurrent - intialEnergyNVE)) / abs(intialEnergyNVE);


		if (current_n != 0 && (current_n % 1000 == 0)) FileError << EnergyErrorNVE.getAverage() << "\n";
		if (current_n != 0 && (current_n % 1000 == 0)) FileErrorInst << errorEnergyPercentage << "\n";

#endif

#ifdef HIST_MAX_NABLA_U



		double pdelta = pSTEP;

		if ((current_n % WritingPeriodFile == 0) && current_n < MaxLengthOfCurrentN){

			std::ofstream FileHistMaxNablaU("histMaxNablaU.m");

			FileHistMaxNablaU << "pMax=" << pMAX << ";" << std::endl;

			FileHistMaxNablaU << "H=[";

			for (int i = 0; i < numberOfBinsHistpMax; i++){

				FileHistMaxNablaU << HistNablaU[i].getAverage() << "\t";
			}
			FileHistMaxNablaU << "];" << std::endl;

			FileHistMaxNablaU << "pHist=[";
			for (int i = 0; i < numberOfBinsHistpMax; i++){

				FileHistMaxNablaU << (double)i*pdelta << "\t";
			}
			FileHistMaxNablaU << "];" << std::endl;


			FileHistMaxNablaU << "NablUInstant=[";
			for (int i = 0; i < numberOfBinsHistpMax; i++){

				FileHistMaxNablaU << NablaUArray[i] << "\t";
			}
			FileHistMaxNablaU << "];" << std::endl;


			FileHistMaxNablaU << "NablUCounter=[";
			for (int i = 0; i < numberOfBinsHistpMax; i++){

				FileHistMaxNablaU << NablaUCounterArray[i].getAverage() << "\t";
			}
			FileHistMaxNablaU << "];" << std::endl;



			std::ofstream FileHistProbaMatrixNablaU("probaMatrixNablaU.m");

			FileHistProbaMatrixNablaU << "delta=" << pSTEP << ";\n";
			FileHistProbaMatrixNablaU << "pMax=" << pMAX << ";\n";
			FileHistProbaMatrixNablaU << "nablaUMax=" << NABLA_U_MAX << ";\n";

			FileHistProbaMatrixNablaU << "P=[\n";

			for (int i = 0; i < numberOfBinsHistNablaU; i++){

				for (int j = 0; j < numberOfBinsHistpMax; j++){


					FileHistProbaMatrixNablaU << ProbaMatrixNablaU[i][j] / probaMatrixNablaUCounterTotal << "\t";


				}
				FileHistProbaMatrixNablaU << "\n";
			}

			FileHistProbaMatrixNablaU << "];" << std::endl;

			//std::cout << "here" << std::endl;
			//while (1);

		}

#endif

#ifdef ADD_SAMPLES



		AddSamples();


#ifdef GAMMA_ANALYSIS

		if (current_n > 100000){



			unsigned int hamIndic = 0;
			unsigned int fdIndic = 0;

//			fixedRejectionRateHam = 0.1;

			double diffRejRate = ((RejectionRate.getAverage() - fixedRejectionRateHam) / fixedRejectionRateHam);
			double diffRejRateFD = (RejectionRateFD.getAverage() - fixedRejectionRateFD) / fixedRejectionRateFD;

			//std::cout << "\n diff rej rate = " << diffRejRate << std::endl;


			if (fabs(diffRejRate) < toleranceRejectionRate && fabs(diffRejRateFD)<toleranceRejectionRate){



				NumberOfTimeSteps = 0.0;

				std::ofstream MaxDt("maxDt.m");

				MaxDt << "EpsF = " << pParticles[0]->getEpsf() << ";\n" << std::endl;
				MaxDt << "EpsR = " << pParticles[0]->getEpsr() << ";\n" << std::endl;
				MaxDt << "RejRate = " << RejectionRate.getAverage() << ";\n" << std::endl;
				MaxDt << "RejRateFD = " << RejectionRateFD.getAverage() << ";\n" << std::endl;
				MaxDt << "mxDt= " << dt << ";\n" << std::endl;
				MaxDt << "gamma= " << gamma << ";\n" << std::endl;


#ifdef WRITE_INSTANT_REJECTION_RATE


					CleanSamples();

					std::ofstream MaxDtInst("rejRateInst.m");
					std::ofstream MaxDtInstFD("rejRateInstFD.m");

					MaxDtInst << "R =[ " << std::endl;
					MaxDtInstFD << "RFD =[ " << std::endl;

					for (int n = 0; n < 10000; n++){

						PerformTimeStep();

						MaxDtInst << instantRejectionHam << std::endl;
						MaxDtInstFD << instantRejectionFD << std::endl;

					}
					MaxDtInst <<  "];\n" << std::endl;
					MaxDtInstFD << "];\n" << std::endl;



#endif

			}
			else{



				if (diffRejRate > toleranceRejectionRate){

					dt = 0.9*dt;

					hamIndic = 1;

					std::cout << "\n\n\n" << (RejectionRate.getAverage() - fixedRejectionRateHam) / fixedRejectionRateHam << ">" << toleranceRejectionRate << std::endl;
					std::cout << "\nRejection Rate Ham " << RejectionRate.getAverage() << std::endl;
					std::cout << "\n decreasing time step size \ndt =" << dt << std::endl;

					CleanSamples();
					current_n = 0;
				}
				if (hamIndic == 0 && diffRejRate < -toleranceRejectionRate){

					dt = 1.1*dt;

					fdIndic = 1;


					std::cout << "\nRejection Rate Ham " << RejectionRate.getAverage() << std::endl;
					std::cout << "\n increasing time step size \ndt =" << dt << std::endl;

					CleanSamples();
					current_n = 0;
				}

				if (fabs(diffRejRate)<toleranceRejectionRate){



					if (diffRejRateFD > toleranceRejectionRate){

						gamma = 0.9*gamma;

						std::cout << "\nRejection Rate FD " << RejectionRateFD.getAverage() << std::endl;
						std::cout << "\n decreasing gamma \ngamma =" << gamma << std::endl;

						CleanSamples();
						current_n = 0;

					}

					if (fdIndic == 0 && diffRejRateFD < -toleranceRejectionRate){

						gamma = 1.1*gamma;

						std::cout << "\nRejection Rate FD " << RejectionRateFD.getAverage() << std::endl;
						std::cout << "\n increasing gamma \ngamma =" << gamma << std::endl;

						CleanSamples();
						current_n = 0;

					}
				}
			}


		}
#endif




#endif

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

		WriteAutocorrelationAveragedInReplicas();

#endif

#ifdef VARIANCE_IR
		SampleVarianceIR();
#endif

#ifdef MORE_VARIANCE_IR
		SampleVarianceIR2();
#endif



		//if (current_n != 0 && ((int)current_n % 1000000) == 0)		SaveCurrentSimulation();
		//if (current_n != 0 && (current_n % WritingPeriodFile == 0))	{ //1000000

		//	SaveCurrentSimulation();
		//}
#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS
		if (current_n != 0 && (current_n % WritingPeriodFile == 0))	{
			SaveCurrentAveragesOnlyInstantOnes();
			//SaveCurrentAverages();
#endif
		}



		WriteOnScreen();

		Tock();

		DisplayTimeLeft();

		if (current_n < MaxLengthOfCurrentN)		current_n = current_n + 1.0;
		else {
			current_n = 0.0;
			numberOfCurrent_n = numberOfCurrent_n + 1;
		}

		//}

	}

	SaveCurrentSimulation();
	//

#ifdef WRITE_DIMER_DISTANCE_MATLAB

	FileDimerDistance << "];\n";

#endif

#ifdef WRITE_ERROR_LANGEVIN


	FileError << "];\n";
	FileErrorInst << "];\n";

#endif

#ifdef BLOCK_AVERAGING

	double* sigmaBA=BlockAveraging();
	double sigmaBA1 = sigmaBA[0]; // dimerDistance variance
	double sigmaBA2 = sigmaBA[1]; // dimerPotential variance

	std::cout << "Discr Standard deviation of DimerDistance from BA is " << sigmaBA1 << std::endl;
	std::cout << "Cont Variance of DimerDistance from BA is " << sigmaBA1*sigmaBA1*dt << std::endl;

	std::cout << "Discr Standard deviation of DimerPotential from BA is " << sigmaBA2 << std::endl;
	std::cout << "Cont Variance of DimerPotential from BA is " << sigmaBA2*sigmaBA2*dt << std::endl;

	if (Var_mu != 0)	{
		std::cout << "Only for DimerDistance:" << std::endl;
		std::cout << "Correlation length N_corr is " << sigmaBA1*sigmaBA1 / (Var_mu*Var_mu) << std::endl;
		std::cout << "Correlation time T_corr is " << sigmaBA1*sigmaBA1*dt / (Var_mu*Var_mu) << std::endl;
	}
	else std::cout << "Var_mu is zero.. (DimerDistance)" << std::endl;
	//std::cout<< "Temperature and pressure written into file!"<<std::endl;

	standardDeviation=sigmaBA1;
	standardDeviation2 = sigmaBA2;

#endif

	if (errorIndicator == true) {
		std::cout << "ERROR IN PBC - TRAJECTORY NOT COMPUTED! CHECK INITIAL CONDITION IF USING LENNARD JONES..\n" << std::endl;
		while (1);
	}



#ifdef RDF

	std::ofstream filerdf("data/RDF");


	for (int i = 0; i<numberOfBins; i++)
	{
		filerdf << i*drRDF<<"\t";
		//fprintf(FileRadialDistribution, "%lf,\t", i*dr);//(i+0.5)*dr);
	}

	filerdf << "\n";
	//fprintf(FileRadialDistribution, "];\ng=[");
	double tmp_g = 0;

	int NrCalc = NumberOfTimeSteps;

	for (int i = 0; i<numberOfBins; i++)
	{

		tmp_g = gRDF[i];
		gRDF[i] = tmp_g / (2.0*3.14159265359*drRDF*NrCalc);//cut_off_rdf*NrCalc*NrRadiuses);//tmp_g/ ((4.0*3.14159265359/3.0)*(CUB(i+1)-CUB(i))*CUB(dr)*rho);//3.14159265359

		filerdf << gRDF[i] << "\t";
		//fprintf(FileRadialDistribution, "%lf,\t", gRDF[i]);

	}

	filerdf << "\n";

#endif

#ifdef DISTRIBUTION_RESTR_PART

	std::ofstream filedrp("data/DistrRestrPartDimerCenter");


	for (int i = 0; i<numberOfBins_DRP; i++)
	{
		filedrp << i*drDRP << "\t";
		//fprintf(FileRadialDistribution, "%lf,\t", i*dr);//(i+0.5)*dr);
	}

	filedrp << "\n";


	for (int i = 0; i<numberOfBins_DRP; i++)		filedrp << gDRP[i].getAverage() << "\t";

	filedrp << "\n";

	for (int i = 0; i<numberOfBins_DRP; i++)		filedrp << gDRPRestrained[i].getAverage() << "\t";

	filedrp << "\n";

	//fprintf(FileRadialDistribution, "];\ng=[");
	//double tmp_gDRP = 0;

	//int NrCalcDRP = NumberOfTimeSteps;

	//for (int i = 0; i<numberOfBins_DRP; i++)
	//{

	//	tmp_gDRP = gDRP[i];
	//	gDRP[i] = tmp_gDRP / (2.0*3.14159265359*drDRP*NrCalcDRP);//cut_off_drp*NrCalcDRP*NrRadiuses);//tmp_g/ ((4.0*3.14159265359/3.0)*(CUB(i+1)-CUB(i))*CUB(dr)*rho);//3.14159265359

	//	filedrp << gDRP[i] << "\t";
	//	//fprintf(FileRadialDistribution, "%lf,\t", gRDF[i]);

	//}



#endif

#ifdef WRITE_POSITIONS_MATLAB
#ifdef ONE_PARTICLE_ONLY

	FilePositions << "];\n" << std::endl;
#endif
#endif

#ifdef WRITE_MOMENTA_MATLAB
#ifdef ONE_PARTICLE_ONLY

	FileMomenta << "];\n" << std::endl;
#endif

#endif


#ifdef WRITE_MOMENTA_MATLAB_HISTOGRAM


	FileMomentaHist << "];\n" << std::endl;


#endif

#ifdef WRITE_INST_ENERGY


	FileEnergy << "];" << std::endl;

#endif


#ifdef WRITE_P_SQUARE_2D

	FilePSquare << "];\n";

#endif


#ifdef HISTOGRAM
	std::ofstream HistrogramP2("histogram_p2.m");

	HistrogramP2 << "H=[\n";



	for (int i = 0; i < numerOfBinsInHistogram; i++){


		HistrogramP2<< ((double)Histogram[i] / (double)histogramCounter)  <<"\n";


	}

	HistrogramP2 << "];";



#endif

#ifdef COMPUTE_MOMENTA_SQUARE_AVERAGE

	std::ofstream FileMomSquare("momentaSquareMean.m");

	FileMomSquare << "pSquare= " << MomentaSquare.getAverage() << ";\n" << std::endl;

#endif

};



void Simulation::PerformTrajectoryChangingParameters(){



#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS

	double howManyTimesAverages = HOW_MANY_PHYSICAL_TIMES;
	double writingIntervalAverages = NumberOfTimeSteps / howManyTimesAverages;
	std::cout << "\n Writing interval for averages =" << writingIntervalAverages << std::endl;

#endif


#ifdef NVE

	PerformTimeStep();
	intialEnergyNVE = ComputeEnergy();// ComputeKineticEnergy();
	std::cout << "Initial energy in NVE is " << intialEnergyNVE << std::endl;


#endif

	current_n = 0;
	measureTime = 0.0;

	while (current_n < NumberOfEquilibrationSteps){

		if (errorIndicator == true) current_n = NumberOfEquilibrationSteps;

		//tPerTimeStep0 = SBCTimeClock();

		PerformTimeStep();

		//tPerTimeStep1 = SBCTimeClock();

		current_n++;

	}

#ifdef NVE

	std::cout << "computing the initial energy after equilibration\n" << std::endl;

	intialEnergyNVE = ComputeEnergy();

#endif

	current_n = 0;// 3 * MaxLengthOfCurrentN;// 2.1*pow(10.0, 9);
	double epsrTMP = pParticles[0]->getEpsr();
	double epsfFixed = pParticles[0]->getEpsf();

	while (epsrTMP < 0.8*epsfFixed){

		for (int i = 0; i < NumberOfParticles; i++){

			pParticles[i]->setEps(epsrTMP, epsfFixed);

		}

		std::cout << "\nCHANGING EPSR \nEpsR = " << epsrTMP << std::endl;

		current_n = 0;

		while ((((double)current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN)) < NumberOfTimeSteps){


			if (errorIndicator == true)  NumberOfTimeSteps = 0.0;

			Tick();

#ifdef MEASURE_TIME_PER_TIMESTEP
			tPerTimeStep3 = SBCTimeClock();
#endif
			PerformTimeStep();

#ifdef MEASURE_TIME_PER_TIMESTEP
			tPerTimeStep4 = SBCTimeClock();
#endif
			DisplayTimeLeft();
			AddSamplesChangingParameters();

			if (current_n != 0 && (current_n % 10000 == 0))	{ //1000000

				SaveCurrentSimulation();
#ifdef SAVE_CURRENT_AVERAGES_FOR_TIME_STEP_ANALYSIS
				SaveCurrentAveragesOnlyInstantOnes();
#endif
			}

			WriteOnScreen();

			Tock();

			DisplayTimeLeft();

			if (current_n < MaxLengthOfCurrentN)		current_n = current_n + 1.0;
			else {
				current_n = 0.0;
				numberOfCurrent_n = numberOfCurrent_n + 1;
			}


		}

		epsrTMP = epsrTMP + 0.1*pParticles[0]->getEpsf();
	}



	if (errorIndicator == true) {
		std::cout << "ERROR IN PBC - TRAJECTORY NOT COMPUTED! CHECK INITIAL CONDITION IF USING LENNARD JONES..\n" << std::endl;
		while (1);
	}




};



void Simulation::PeriodicBoundaryCondition() {

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n==0) 	std::cout << "PeriodicBoundaryCondition" << std::endl;

#endif

	for (int i = 0; i < NumberOfParticles; i++) {

		if (pParticles[i]->getActiveStatus() == 0) continue;

		double qx = pParticles[i]->getPositionX();
		double qy = pParticles[i]->getPositionY();
		double qz = pParticles[i]->getPositionZ();

		if (qx <= 0.0)
		{
			qx = qx + BoxLength;

		}

		if (qx > BoxLength)
		{
			qx = qx - BoxLength;

		}

#ifdef ERROR_ALERT_PBC
		if (qx <= 0.0)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;

		}

		if (qx > BoxLength)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;

		}
#endif
		if (qy <= 0.0)
		{
			qy = qy + BoxLength;

		}

		if (qy > BoxLength)
		{
			qy = qy - BoxLength;

		}

#ifdef ERROR_ALERT_PBC

		if (qy <= 0.0)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;
		}

		if (qy > BoxLength)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;
		}
#endif
		if (qz <= 0.0)
		{
			qz = qz + BoxLength;

		}

		if (qz > BoxLength)
		{
			qz = qz - BoxLength;

		}

#ifdef ERROR_ALERT_PBC
		if (qz <= 0.0)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;
		}

		if (qz > BoxLength)
		{
			std::cout << "**************************************************************ERROR PBC" << std::endl;
			errorIndicator = true;

		}
#endif

		pParticles[i]->setPosition(qx, qy, qz);
		grid->onPositionChanged(pParticles[i]);


	}

	if (errorIndicator == true) {

		std::cout << "Error in PBC\n" << std::endl;
		WriteOnScreen();
		//	while(1);
	}

	if (errorIndicator == true){


		grid->print();

		std::cout << "In PeriodicBoundaryConditions(), after explosion: how many active= " << Particle::howManyActive << std::endl;
		while (1);

	}

}


void Simulation::Tick(){

	t0 = SBCTimeClock() / C_PROC;

};

void Simulation::Tock(){

	t1 = SBCTimeClock() / C_PROC;
	//measureTime=measureTime + (t1-t0);
	MeasureTime.addSample(t1 - t0);
};

void Simulation::UpdateMomentaParticles(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaParticles "<<t << std::endl;
#endif

	for (int i = 0; i < NumberOfParticles; i++){

		pParticles[i]->setMomentum(pParticles[i]->getMomentumX() - pParticles[i]->getFX() * t, pParticles[i]->getMomentumY() - pParticles[i]->getFY() * t, pParticles[i]->getMomentumZ() - pParticles[i]->getFZ() * t);


	}



};

void	Simulation::UpdatePositionsParticles(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdatePositionsParticles "<<t << std::endl;
#endif


	for (int i = 0; i < NumberOfParticles; i++){

		if (pParticles[i]->getActiveStatus() == 0) continue;


			double epsri = pParticles[i]->getEpsr();
			double epsfi = pParticles[i]->getEpsf();

			double qx = pParticles[i]->getPositionX();
			double qy = pParticles[i]->getPositionY();
			double qz = pParticles[i]->getPositionZ();

			double px = pParticles[i]->getMomentumX();
			double py = pParticles[i]->getMomentumY();
			double pz = pParticles[i]->getMomentumZ();




			double mass = pParticles[i]->getMass();

			double dHdpx = dHdp(px, py, pz, mass, epsri, epsfi);
			double dHdpy = dHdp(py, px, pz, mass, epsri, epsfi);
			double dHdpz = dHdp(pz, py, px, mass, epsri, epsfi);


			pParticles[i]->setPosition(qx + dHdpx * t, qy + dHdpy * t, qz + dHdpz * t);
			grid->onPositionChanged(pParticles[i]);

			//		temperatureMomenta = +px*dHdpx + py*dHdpy;


	}

};


void Simulation::UpdatePositionsParticlesStandard(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdatePositionsParticlesStandard "<<t << std::endl;
#endif


	for (int i = 0; i < NumberOfParticles; i++){

		double qx = pParticles[i]->getPositionX();
		double qy = pParticles[i]->getPositionY();
		double qz = pParticles[i]->getPositionZ();

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double mass = pParticles[i]->getMass();

		double dHdpx = px / mass;
		double dHdpy = py / mass;
		double dHdpz = pz / mass;

		//	double px=pParticles[i]->getPositionX() + dHdp(pParticles[i]->getMomentumX(),pParticles[i]->getMomentumY(),pParticles[i]->getMass(), pParticles[i]->getEpsr(),pParticles[i]->getEpsf())* t;
		//	double py=pParticles[i]->getPositionY() + dHdp(pParticles[i]->getMomentumY(),pParticles[i]->getMomentumX(),pParticles[i]->getMass(), pParticles[i]->getEpsr(),pParticles[i]->getEpsf())* t;

		pParticles[i]->setPosition(qx + dHdpx * t, qy + dHdpy * t, qz + dHdpz * t);
		grid->onPositionChanged(pParticles[i]);

		//		temperatureMomenta = +px*dHdpx + py*dHdpy;
	}



};




void Simulation::UpdateMomentaStochOrder1(double deltat){



#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochOrder1 "<<deltat << std::endl;
#endif

	temperatureMomenta = 0.0;

	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1);
		float Gy = ComputeBoxMuller(0, 1);
		float Gz = ComputeBoxMuller(0, 1);

		double c2x = sigma*Gx*sqrt(deltat);
		double c2y = sigma*Gy*sqrt(deltat);
		double c2z = sigma*Gz*sqrt(deltat);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();


		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();



		double dHdpx = dHdp(px, py, pz, mass, epsri, epsfi);
		double dHdpy = dHdp(py, px, pz, mass, epsri, epsfi);
		double dHdpz = dHdp(pz, px, py, mass, epsri, epsfi);


		px = px - gamma*dHdpx*deltat + c2x;
		py = py - gamma*dHdpy*deltat + c2y;
		pz = pz - gamma*dHdpz*deltat + c2z;

		temperatureMomenta = +px*dHdpx + py*dHdpy+ pz*dHdpz;



		pParticles[i]->setMomentum(px, py,pz);


	}


}


void Simulation::UpdateMomentaStochOrder2_CN(double deltat){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochOrder2_CN "<<deltat << std::endl;
#endif

	temperatureMomenta = 0.0;

	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1);
		float Gy = ComputeBoxMuller(0, 1);
		float Gz = ComputeBoxMuller(0, 1);

		double c2x = sigma*Gx*sqrt(deltat);
		double c2y = sigma*Gy*sqrt(deltat);
		double c2z = sigma*Gz*sqrt(deltat);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();


		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();

		double px_mid = px - gamma*0.5*deltat*dHdp(px, py, pz, mass, epsri, epsfi) + 0.5*c2x;
		double py_mid = py - gamma*0.5*deltat*dHdp(py, px, pz, mass, epsri, epsfi) + 0.5*c2y;
		double pz_mid = pz - gamma*0.5*deltat*dHdp(pz, px, pz, mass, epsri, epsfi) + 0.5*c2z;

		double dHdpx = dHdp(px_mid, py_mid, pz_mid, mass, epsri, epsfi);
		double dHdpy = dHdp(py_mid, px_mid, pz_mid, mass, epsri, epsfi);
		double dHdpz = dHdp(pz_mid, px_mid, py_mid, mass, epsri, epsfi);

		////********************************************
		//double dHdpx = px*(1.0 - Rho_p(px, py, epsri, epsfi, mass)) / mass - (px*px + py*py)*dRho_p(px, py, epsri, epsfi, mass) / (2.0*mass);
		//double dHdpy = py*(1.0 - Rho_p(px, py, epsri, epsfi, mass)) / mass - (px*px + py*py)*dRho_p(py, px, epsri, epsfi, mass) / (2.0*mass);
		//
		//double Simulation::dHdp(double px, double py, double mass, double epsrp, double epsfp){
		//
		//	//if(epsfp==0 && epsrp == 0 && Rho_p(px,epsrp,epsfp,mass) !=0) std::cout<<"Rho !=0\n";
		//
		//	return  px*(1.0 - Rho_p(px, py, epsrp, epsfp, mass)) / mass - (px*px + py*py)*dRho_p(px, py, epsrp, epsfp, mass) / (2.0*mass);// (px/mass)*(1.0-Rho_p(px,py,epsrp,epsfp,mass)) -((px*px+py*py)/(2.0*mass))*dRho_p(px,py,epsrp,epsfp,mass);
		//
		//};
		////********************************************

		px = px - gamma*dHdpx*deltat + c2x;
		py = py - gamma*dHdpy*deltat + c2y;
		pz = pz - gamma*dHdpz*deltat + c2z;

		//temperatureMomenta = +px*dHdpx + py*dHdpy+ pz*dHdpz;

		pParticles[i]->setMomentum(px, py, pz);


	}


}


void	Simulation::UpdateMomentaStochExact(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochExact "<< t << std::endl;
#endif



	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
		float Gy = ComputeBoxMuller(0, 1);
		float Gz = ComputeBoxMuller(0, 1);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();


		double mass = pParticles[i]->getMass();

		double alpha = exp(-gamma*t / mass);

		px = alpha*px + sqrt((1 - alpha*alpha)*kT)*Gx;
		py = alpha*py + sqrt((1 - alpha*alpha)*kT)*Gy;
		pz = alpha*pz + sqrt((1 - alpha*alpha)*kT)*Gz;


		pParticles[i]->setMomentum(px, py, pz);


	}


};

double Simulation::AnalyticalSolutionFD(double t, double px, double mass, double epsrp, double epsfp){


	float Gx = ComputeBoxMuller(0, 1);

	if (fabs(px) >= epsfp)
	{

		double alpha = exp(-gamma*t / mass);
		return alpha*px + sqrt((1 - alpha*alpha)*kT)*Gx;

	}
	if (fabs(px) <= epsrp)
	{

		return px + sqrt(2*gamma*kT*t)*Gx;
	}

	if (px < 0){


		double b = epsfp;
		double a = epsrp;

		double alpha = exp(-gamma*b*t / (mass*(b - a)));


		//return alpha*(px + a) + sqrt(mass*(b - a)*(1-alpha*alpha ) / (b*beta))*Gx ;
		return alpha*(px)+a*(alpha - 1) + sqrt(mass*(b - a)*(1 - alpha*alpha) / (b*beta))*Gx;
	}
	else
	{
		double b = epsfp;
		double a = epsrp;

		double alpha = exp(-gamma*b*t / (mass*(b - a)));

	//	return alpha*(px - a) + sqrt(mass*(b - a)*(1-alpha*alpha) / (b*beta))*Gx;
		return alpha*(px)-a*(alpha-1) + sqrt(mass*(b - a)*(1 - alpha*alpha) / (b*beta))*Gx;

	}



}

void	Simulation::UpdateMomentaStochExactARPS(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochExactARPS " << t << std::endl;
#endif



	for (int i = 0; i < NumberOfParticles; i++){


		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double mass = pParticles[i]->getMass();


		px = AnalyticalSolutionFD(t, px, mass, epsri, epsfi);
		py = AnalyticalSolutionFD(t, py, mass, epsri, epsfi);
		pz = AnalyticalSolutionFD(t, pz, mass, epsri, epsfi);

		pParticles[i]->setMomentum(px, py, pz);


	}


};



void Simulation::UpdateMomentaStochModifiedEqCN(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochModifiedEqCN" << std::endl;
#endif

	int *activeStatus = new int[3];
	double *mom = new double[3];

	/*-(1/4)*gamma^2*((D^2)(dU))(p)*dt^2/beta*/

	for (int i = 0; i < NumberOfParticles; i++){

		//pParticles[i]->updateActiveStatus();

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double mass  = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();


		activeStatus[0] = pParticles[i]->decideActiveNonActiveOneCoordinate(px);
		activeStatus[1] = pParticles[i]->decideActiveNonActiveOneCoordinate(py);
		activeStatus[2] = pParticles[i]->decideActiveNonActiveOneCoordinate(pz);


		mom[0] = px;
		mom[1] = py;
		mom[2] = pz;

		double additionalForcing = 0.0;

		for (int dim = 0; dim < d; dim++){

			if (activeStatus[dim] == 2){

				double F1x  = dHdp3(mom[dim], py, pz, mass, epsri, epsfi);

				mom[dim] = mom[dim]  -0.25*gamma *gamma  * F1x*t*t*kT;
						}
		}

		pParticles[i]->setMomentum(mom[0], mom[1], mom[2]);

	}

	delete mom;
	delete activeStatus;

}


void Simulation::UpdateMomentaStochSecondOrder(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochSecondOrder" << std::endl;
#endif


	int *activeStatus = new int[3];
	double *mom = new double[3];
	double *correction = new double[3];
	correction[0] = 0;
	correction[1] = 0;
	correction[2] = 0;

	//---------- EXPLICIT MIDPOINT


	temperatureMomenta = 0.0;

	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1);
		float Gy = ComputeBoxMuller(0, 1);
		float Gz = ComputeBoxMuller(0, 1);

		double c2x = sigma*Gx*sqrt(t);
		double c2y = sigma*Gy*sqrt(t);
		double c2z = sigma*Gz*sqrt(t);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();

		double px_mid = px - gamma*0.5*t*dHdp(px, py, pz, mass, epsri, epsfi) + 0.5*c2x;
		double py_mid = py - gamma*0.5*t*dHdp(py, px, pz, mass, epsri, epsfi) + 0.5*c2y;
		double pz_mid = pz - gamma*0.5*t*dHdp(pz, px, pz, mass, epsri, epsfi) + 0.5*c2z;

		double dHdpx = dHdp(px_mid, py_mid, pz_mid, mass, epsri, epsfi);
		double dHdpy = dHdp(py_mid, px_mid, pz_mid, mass, epsri, epsfi);
		double dHdpz = dHdp(pz_mid, px_mid, py_mid, mass, epsri, epsfi);


	//------SECOND ORDER CORRECTION TERM------------------------------------------------------------

	/*-(1/4)*gamma^2*((D^2)(dU))(p)*dt^2/beta*/


		activeStatus[0] = pParticles[i]->decideActiveNonActiveOneCoordinate(px);
		activeStatus[1] = pParticles[i]->decideActiveNonActiveOneCoordinate(py);
		activeStatus[2] = pParticles[i]->decideActiveNonActiveOneCoordinate(pz);


		mom[0] = px;
		mom[1] = py;
		mom[2] = pz;

		double additionalForcing = 0.0;

		for (int dim = 0; dim < d; dim++){

			if (activeStatus[dim] == 2){

				double F1x = dHdp3(mom[dim], py, pz, mass, epsri, epsfi);

				correction[dim] = - 0.25*gamma *gamma  * F1x*t*t*kT;

			}
		}

		px = px - gamma*dHdpx*t + c2x + correction[0];
		py = py - gamma*dHdpy*t + c2y + correction[1];
		pz = pz - gamma*dHdpz*t + c2z + correction[2];

		//temperatureMomenta = +px*dHdpx + py*dHdpy+ pz*dHdpz;

		pParticles[i]->setMomentum(px, py, pz);

	//	pParticles[i]->setMomentum(mom[0], mom[1], mom[2]);

	}

	delete mom;
	delete activeStatus;
	delete correction;

}



void	Simulation::UpdateMomentaARPSperturbation(double t){




#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaARPSperturbation "<<t << std::endl;
#endif


	for (int i = 0; i < NumberOfParticles; i++){

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();

		double dHdpx = dHdp(px, py, pz, mass, epsri, epsfi); //px*(1.0 - Rho_p(px, py,pz, epsri, epsfi, mass)) / mass - (px*px + py*py +pz*pz)*dRho_p(px, py, epsri, epsfi, mass) / (2.0*mass);
		double dHdpy = dHdp(py, px, pz, mass, epsri, epsfi);
		double dHdpz = dHdp(pz, px, py, mass, epsri, epsfi);

		px = px - gamma*(dHdpx - px)*t;
		py = py - gamma*(dHdpy - py)*t;
		pz = pz - gamma*(dHdpz - pz)*t;

		pParticles[i]->setMomentum(px, py, pz);


	}



};


void	Simulation::UpdateMomentaVerlet_FDR(double h){

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaVerlet_FDR "<< h << std::endl;
#endif


#ifndef METROPOLIZATION
	SampleG();
#endif

	for (int i = 0; i < NumberOfParticles; i++){

		double Gx = pParticles[i]->GetGx();//ComputeBoxMuller(0, 1);
		double Gy = pParticles[i]->GetGy();//ComputeBoxMuller(0, 1);
		double Gz = pParticles[i]->GetGz();//ComputeBoxMuller(0, 1);

		double Rx = sqrt(kT)*Gx;
		double Ry = sqrt(kT)*Gy;
		double Rz = sqrt(kT)*Gz;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double mass = pParticles[i]->getMass();
		double epsr = pParticles[i]->getEpsr();
		double epsf = pParticles[i]->getEpsf();

		px = px + 0.5*h*Rx;
		Rx = Rx - h*dHdp(px, py, pz, mass, epsr, epsf);
		px = px + 0.5*h*Rx;

		py = py + 0.5*h*Ry;
		Ry = Ry - h*dHdp(py, px, pz, mass, epsr, epsf);
		py = py + 0.5*h*Ry;

		pz = pz + 0.5*h*Rz;
		Rz = Rz - h*dHdp(pz, py, px, mass, epsr, epsf);
		pz = pz + 0.5*h*Rz;

		pParticles[i]->setMomentum(px,py,pz);
		pParticles[i]->setG(sqrt(beta)*Rx, sqrt(beta)*Ry, sqrt(beta)*Rz);
	}




};




void Simulation::MetropolisFDR(double* T_old, double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "MetropolisFDR "<< t << std::endl;
#endif


	double* T_new = ComputeEnergyFDR();


	for (int i = 0; i<NumberOfParticles; i++) {

	double rejectionIndicator = 0;

	double alpha = rMetropolis->randDouble2(); // random number in [0,1)

	double alphaDt = beta*(T_new[i] - T_old[i]); // exp( -beta*(new_energy-X.old_energy(k)) );
	double acceptanceRate = exp(-alphaDt);


	//absValueAlphaAverageMetropolis.addSample(abs(alphaDt));
	//alphaAverageMetropolis.addSample(alphaDt);

	//if (alpha > alphaDt/(1+alphaDt)){
	if (acceptanceRate < alpha){

		rejectionIndicator = 1;

		pParticles[i]->setMomentum(-pParticles[i]->getMomentumXOld(), -pParticles[i]->getMomentumYOld(), -pParticles[i]->getMomentumZOld());

		}

	RejectionRate.addSample(rejectionIndicator);

	}

};




void	Simulation::UpdateMomentaStochLangevin(double t){

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochLangevin" << std::endl;
#endif

	double h = t;// sqrt(2 * gamma*t);


	double rejectionIndicator = 0;

	for (int i = 0; i < NumberOfParticles; i++){

		//	for (int i = 0; i < 1; i++){

		rejectionIndicator = 0;


		double Gx = ComputeBoxMuller(0, 1);
		double Gy = ComputeBoxMuller(0, 1);
		double Gz = ComputeBoxMuller(0, 1);

		double Rx = pParticles[i]->getRX();//sqrt(kT)*Gx;
		double Ry = pParticles[i]->getRY();//sqrt(kT)*Gy;
		double Rz = pParticles[i]->getRZ();// sqrt(kT)*Gz;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double mass = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double alpha = exp(-gamma*t / mass);

		Rx = Rx - 0.5*h*dHdp(px,py,pz,mass,epsri,epsfi);
		px = px + h*Rx;
		Rx = Rx - 0.5*h*dHdp(px, py, pz, mass, epsri, epsfi);
		Rx = alpha*Rx + sqrt((1 - alpha*alpha)*kT)*Gx;

		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);
		py = py + h*Ry;
		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);
		Ry = alpha*Ry + sqrt((1 - alpha*alpha)*kT)*Gy;

		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);
		pz = pz + h*Rz;
		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);
		Rz = alpha*Rz + sqrt((1 - alpha*alpha)*kT)*Gz;

		pParticles[i]->setMomentum(px, py, pz);
		pParticles[i]->setR(Rx, Ry, Rz);
	}



}



void	Simulation::MetropolisStochLangevin(double t){

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "MetropolisStochLangevin "<< t << std::endl;
#endif

	double h = t;// sqrt(2 * gamma*t);




	for (int i = 0; i < NumberOfParticles; i++){


		int rejectionIndicator = 0;

		double Gx = ComputeBoxMuller(0, 1);
		double Gy = ComputeBoxMuller(0, 1);
		double Gz = ComputeBoxMuller(0, 1);

		double Rx = pParticles[i]->getRX();//sqrt(kT)*Gx;
		double Ry = pParticles[i]->getRY();//sqrt(kT)*Gy;
		double Rz = pParticles[i]->getRZ();// sqrt(kT)*Gz;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double RxOld = Rx;
		double RyOld = Ry;
		double RzOld = Rz;

		double pxOld = px;
		double pyOld = py;
		double pzOld = pz;

		double mass = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double EOld = Hp(px, py, pz, mass, epsri, epsfi) +0.5*(Rx*Rx+Ry*Ry+Rz*Rz);

		//Verlet

		Rx = Rx - 0.5*h*dHdp(px, py, pz, mass, epsri, epsfi);
		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);
		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);

		px = px + h*Rx;
		py = py + h*Ry;
		pz = pz + h*Rz;

		Rx = Rx - 0.5*h*dHdp(px, py, pz, mass, epsri, epsfi);
		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);
		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);

		pParticles[i]->setMomentum(px, py, pz);
		pParticles[i]->setR(Rx, Ry, Rz);

		double ENew = Hp(px, py, pz, mass, epsri, epsfi) + 0.5*(Rx*Rx + Ry*Ry + Rz*Rz);

		double alphaMetropolis = rMetropolis->randDouble2();

		double alphaDt = beta*(ENew - EOld); // exp( -beta*(new_energy-X.old_energy(k)) );
		double acceptanceRate = exp(-alphaDt);

		if (acceptanceRate < alphaMetropolis){

			rejectionIndicator = 1;

			pParticles[i]->setMomentum(pxOld, pyOld, pzOld);
			pParticles[i]->setR(-RxOld, -RyOld, -RzOld);

		}

		RejectionRateFD.addSample(rejectionIndicator);

		// exact integration

		double alpha = exp(-gamma*t / mass);

		Rx = pParticles[i]->getRX();
		Ry = pParticles[i]->getRY();
		Rz = pParticles[i]->getRZ();

		Rx = alpha*Rx + sqrt((1 - alpha*alpha)*kT)*Gx;
		Ry = alpha*Ry + sqrt((1 - alpha*alpha)*kT)*Gy;
		Rz = alpha*Rz + sqrt((1 - alpha*alpha)*kT)*Gz;

		pParticles[i]->setR(Rx, Ry, Rz);

	}



}


void	Simulation::MetropolisStochLangevinSeparableSpace(double t){

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "MetropolisStochLangevinSeparableSpace" << std::endl;
#endif

	double h = t;

	for (int i = 0; i < NumberOfParticles; i++){


		int rejectionIndicatorX = 0;
		int rejectionIndicatorY = 0;
		int rejectionIndicatorZ = 0;

		double Gx = ComputeBoxMuller(0, 1);
		double Gy = ComputeBoxMuller(0, 1);
		double Gz = ComputeBoxMuller(0, 1);

		double Rx = pParticles[i]->getRX();//sqrt(kT)*Gx;
		double Ry = pParticles[i]->getRY();//sqrt(kT)*Gy;
		double Rz = pParticles[i]->getRZ();// sqrt(kT)*Gz;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double RxOld = Rx;
		double RyOld = Ry;
		double RzOld = Rz;

		double pxOld = px;
		double pyOld = py;
		double pzOld = pz;

		double mass = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double EOldX = Hpx(px, mass, epsri, epsfi) + 0.5*(Rx*Rx);
		double EOldY = Hpx(py, mass, epsri, epsfi) + 0.5*(Ry*Ry);
		double EOldZ = Hpx(pz, mass, epsri, epsfi) + 0.5*(Rz*Rz);

		//Verlet

		Rx = Rx - 0.5*h*dHdp(px, py, pz, mass, epsri, epsfi);
		px = px + h*Rx;
		Rx = Rx - 0.5*h*dHdp(px, py, pz, mass, epsri, epsfi);


		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);
		py = py + h*Ry;
		Ry = Ry - 0.5*h*dHdp(py, px, pz, mass, epsri, epsfi);

		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);
		pz = pz + h*Rz;
		Rz = Rz - 0.5*h*dHdp(pz, py, px, mass, epsri, epsfi);

	//	pParticles[i]->setMomentum(px, py, pz);
	//	pParticles[i]->setR(Rx, Ry, Rz);

		double ENewX = Hpx(px,  mass, epsri, epsfi) + 0.5*(Rx*Rx);
		double ENewY = Hpx(py, mass, epsri, epsfi) + 0.5*(Ry*Ry);
		double ENewZ = Hpx(pz, mass, epsri, epsfi) + 0.5*(Rz*Rz);

		double alphaMetropolisX = rMetropolis->randDouble2();

		double alphaDtX = beta*(ENewX - EOldX); // exp( -beta*(new_energy-X.old_energy(k)) );
		double acceptanceRateX = exp(-alphaDtX);

		if (acceptanceRateX < alphaMetropolisX){

			rejectionIndicatorX = 1;

			px = pxOld;
			Rx = -RxOld;

		}

		double alphaMetropolisY = rMetropolis->randDouble2();

		double alphaDtY = beta*(ENewY - EOldY); // eYp( -beta*(new_energy-Y.old_energy(k)) );
		double acceptanceRateY = exp(-alphaDtY);

		if (acceptanceRateY < alphaMetropolisY){

			rejectionIndicatorY = 1;

			py = pyOld;
			Ry = -RyOld;

		}

		double alphaMetropolisZ = rMetropolis->randDouble2();

		double alphaDtZ = beta*(ENewZ - EOldZ); // eZp( -beta*(new_energy-Z.old_energy(k)) );
		double acceptanceRateZ = exp(-alphaDtZ);

		if (acceptanceRateZ < alphaMetropolisZ){

			rejectionIndicatorZ = 1;

			pz = pzOld;
			Rz = -RzOld;

		}

		pParticles[i]->setMomentum(px, py, pz);
		pParticles[i]->setR(Rx, Ry, Rz);


		RejectionRateFD.addSample(rejectionIndicatorX);
		RejectionRateFD.addSample(rejectionIndicatorY);
		RejectionRateFD.addSample(rejectionIndicatorZ);

		// exact integration

		double alpha = exp(-gamma*t / mass);

		Rx = pParticles[i]->getRX();
		Ry = pParticles[i]->getRY();
		Rz = pParticles[i]->getRZ();

		Rx = alpha*Rx + sqrt((1 - alpha*alpha)*kT)*Gx;
		Ry = alpha*Ry + sqrt((1 - alpha*alpha)*kT)*Gy;
		Rz = alpha*Rz + sqrt((1 - alpha*alpha)*kT)*Gz;

		pParticles[i]->setR(Rx, Ry, Rz);

	}





}



void	Simulation::UpdateMomentaVerlet_FDR_MetropolizePerParticle(double t){

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaVerlet_FDR_MetropolizePerParticle" << std::endl;
#endif

	double h = sqrt(2 * gamma*t);

#ifndef METROPOLIZATION
	SampleG();
#endif

	double rejectionIndicator = 0;

	for (int i = 0; i < NumberOfParticles; i++){

		//	for (int i = 0; i < 1; i++){

		rejectionIndicator = 0;


		double Gx = pParticles[i]->GetGx();//ComputeBoxMuller(0, 1);
		double Gy = pParticles[i]->GetGy();//ComputeBoxMuller(0, 1);
		double Gz = pParticles[i]->GetGz();//ComputeBoxMuller(0, 1);

		double Rx = sqrt(kT)*Gx;
		double Ry = sqrt(kT)*Gy;
		double Rz = sqrt(kT)*Gz;

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double mass = pParticles[i]->getMass();
		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();

		double pxOld = px;
		double pyOld = py;
		double pzOld = pz;

		double OldKinEn = Hp(px,py,pz,mass,epsri,epsfi) + 0.5*(Gx*Gx + Gy*Gy + Gz*Gz);

		px = px + 0.5*h*Rx;
		Rx = Rx - h*dHdp(px, py, pz, mass, epsri, epsfi);
		px = px + 0.5*h*Rx;

		py = py + 0.5*h*Ry;
		Ry = Ry - h*dHdp(py, px, pz, mass, epsri, epsfi);
		py = py + 0.5*h*Ry;

		pz = pz + 0.5*h*Rz;
		Rz = Rz - h*dHdp(pz, py, px, mass, epsri, epsfi);
		pz = pz + 0.5*h*Rz;

		pParticles[i]->setMomentum(px, py, pz);
		pParticles[i]->setG(sqrt(beta)*Rx, sqrt(beta)*Ry, sqrt(beta)*Ry);

		Gx = sqrt(beta)*Rx;
		Gy = sqrt(beta)*Ry;
		Gz = sqrt(beta)*Rz;


		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		double NewKinEn = Hp(px, py, pz, mass, epsri, epsfi) + 0.5*(Gx*Gx + Gy*Gy + Gz*Gz);

		double alpha = rMetropolis->randDouble2(); // random number in [0,1)

		double alphaDt = beta*(NewKinEn - OldKinEn); // exp( -beta*(new_energy-X.old_energy(k)) );
		double acceptanceRate = exp(-alphaDt);


		//absValueAlphaAverageMetropolis.addSample(abs(alphaDt));
		//alphaAverageMetropolis.addSample(alphaDt);

		//if (alpha > alphaDt/(1+alphaDt)){


		// why didnt I decouple the space directions? -> because of the original spline... in the new spline (gradient interpolation)
		// one could accept/reject each direction separately !!!!

		if (acceptanceRate < alpha){

			rejectionIndicator = 1;

			pParticles[i]->setMomentum(pxOld, pyOld, pzOld);
			//pParticles[i]->setMomentum(-pxOld, -pyOld, -pzOld);

		}

		RejectionRateFD.addSample(rejectionIndicator);


#ifdef WRITE_INSTANT_REJECTION_RATE

		instantRejectionFD = rejectionIndicator;

#endif

	}



}



void Simulation::SampleG(){

	for (int i = 0; i < NumberOfParticles; i++){

		pParticles[i]->setG(ComputeBoxMuller(0, 1), ComputeBoxMuller(0, 1), ComputeBoxMuller(0, 1));
		pParticles[i]->setGOld();
	}


};



void Simulation::UpdateMomentaFD_BDV(double t){



	// Schema for fluctuation-dissipation part from Fathi/Stoltz



#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaFD_BDV "<<t << std::endl;
#endif


	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1) *sqrt(kT); // G^n/sqrt(beta)
		float Gy = ComputeBoxMuller(0, 1) *sqrt(kT); // G^n/sqrt(beta)
		float Gz = ComputeBoxMuller(0, 1) *sqrt(kT); // G^n/sqrt(beta)

		double h = sqrt(2.0 * gamma*t);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();

		double pnewTmp_x = px + h*Gx / 2.0;
		double pnewTmp_y = py + h*Gy / 2.0;
		double pnewTmp_z = pz + h*Gy / 2.0;

		double dHdpx = dUdpx(pnewTmp_x, pnewTmp_y, pnewTmp_z, epsri, epsfi, mass);
		double dHdpy = dUdpy(pnewTmp_x, pnewTmp_y, pnewTmp_z, epsri, epsfi, mass);
		double dHdpz = dUdpz(pnewTmp_x, pnewTmp_y, pnewTmp_z, epsri, epsfi, mass);


		px = px + h*Gx - h*h*dHdpx / 2.0;
		py = py + h*Gy - h*h*dHdpy / 2.0;
		pz = pz + h*Gz - h*h*dHdpz / 2.0;


		pParticles[i]->setMomentum(px, py, pz);


	}

};

double Simulation::dUdpx(double px, double py, double pz, double epsri, double epsfi, double mass){

	return px*(1.0 - Rho_p(px, py, pz,epsri, epsfi, mass)) / mass - (px*px + py*py + pz*pz)*dRho_p(px, py,pz, epsri, epsfi, mass) / (2.0*mass);

};

double Simulation::dUdpy(double px, double py, double pz, double epsri, double epsfi, double mass){

	return py*(1.0 - Rho_p(px, py, pz, epsri, epsfi, mass)) / mass - (px*px + py*py +pz*pz)*dRho_p(py, px, pz,epsri, epsfi, mass) / (2.0*mass);

};


double Simulation::dUdpz(double px, double py, double pz, double epsri, double epsfi, double mass){

	return pz*(1.0 - Rho_p(px, py, pz, epsri, epsfi, mass)) / mass - (px*px + py*py + pz*pz)*dRho_p(pz, px, py, epsri, epsfi, mass) / (2.0*mass);

};


void Simulation::UpdateMomentaStochModifiedEq(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochModifiedEq "<<t << std::endl;
#endif

	std::cout << "UpdateMomentaStochModifiedEq has to rewritten for new arps definitionn\n " << t << std::endl;


	for (int i = 0; i < NumberOfParticles; i++){

		pParticles[i]->updateActiveStatus();
		unsigned int activeStatus = pParticles[i]->getActiveStatus();

		switch (activeStatus){

		case 0: break;

		case 1:{

				   float Gx = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
				   float Gy = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
				   float Gz = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)

				   double px = pParticles[i]->getMomentumX();
				   double py = pParticles[i]->getMomentumX();
				   double pz = pParticles[i]->getMomentumX();

				   double mass = pParticles[i]->getMass();
				   //double epsri = pParticles[i]->getEpsr();
				   //double epsfi = pParticles[i]->getEpsf();

				   double Ux = px / mass;// dHdp(px, py, pz, mass, epsri, epsfi);
				   double Uy = py / mass;// dHdp(py, px, pz, mass, epsri, epsfi);
				   double Uz = pz / mass;// dHdp(pz, py, px, mass, epsri, epsfi);


				   double Uxx = 1 / mass;// dHdp2xx(px, py, pz, mass, epsri, epsfi);
				   double Uxy = 0.0;// dHdp2xy(px, py, pz, mass, epsri, epsfi);
				   double Uxz = 0.0;// dHdp2xy(px, pz, py, mass, epsri, epsfi);
				   double Uyy = 1 / mass;// dHdp2xx(py, px, pz, mass, epsri, epsfi);
				   double Uzz = 1.0 / mass;// dHdp2xx(pz, py, px, mass, epsri, epsfi);
				   double Uyz = 0.0;// dHdp2xy(py, pz, px, mass, epsri, epsfi);


				   double F2x = Uxx*Ux;
				   double F2y = Uyy*Uy;
				   double F2z = Uzz*Uz;

				   //dpy dp^2xx in std case is 0



				   double F1x = 0.0;
				   double F1y = 0.0;
				   double F1z = 0.0;

				   double Fx = gamma*(F1x - gamma*F2x) / 6;
				   double Fy = gamma*(F1y - gamma*F2y) / 6;
				   double Fz = gamma*(F1z - gamma*F2z) / 6;


				   //c compute determinant of Id+dt*gamma*nabla2U/3

				   double t2 = t*t;

				   double invB = 1 / (1 + t*gamma / (mass * 3));

				   double sigmaMatrix11 = sqrt(invB) - 1;
				   double sigmaMatrix22 = sqrt(invB) - 1;
				   double sigmaMatrix33 = sqrt(invB) - 1;

				   double Sx = sigmaMatrix11*Gx;
				   double Sy = sigmaMatrix22*Gy;
				   double Sz = sigmaMatrix33*Gz;

				   double sq = sqrt(2 * gamma*t);


				   px = px + t2*Fx + sq*Sx;
				   py = py + t2*Fy + sq*Sy;
				   px = px + t2*Fz + sq*Sz;
		}
			break;

		case 2:{

				   float Gx = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
				   float Gy = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
				   float Gz = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)

				   double px = pParticles[i]->getMomentumX();
				   double py = pParticles[i]->getMomentumX();
				   double pz = pParticles[i]->getMomentumX();

				   double mass = pParticles[i]->getMass();
				   double epsri = pParticles[i]->getEpsr();
				   double epsfi = pParticles[i]->getEpsf();

				   double Ux = dHdp(px, py, pz, mass, epsri, epsfi);
				   double Uy = dHdp(py, px, pz, mass, epsri, epsfi);
				   double Uz = dHdp(pz, py, px, mass, epsri, epsfi);


				   double Uxx = dHdp2xx(px, py, pz, mass, epsri, epsfi);
				   double Uxy = dHdp2xy(px, py, pz, mass, epsri, epsfi);
				   double Uxz = dHdp2xy(px, pz, py, mass, epsri, epsfi);
				   double Uyy = dHdp2xx(py, px, pz, mass, epsri, epsfi);
				   double Uzz = dHdp2xx(pz, py, px, mass, epsri, epsfi);
				   double Uyz = dHdp2xy(py, pz, px, mass, epsri, epsfi);


				   double F2x = Uxx*Ux + Uxy*Uy + Uxz*Uz;
				   double F2y = Uxy*Ux + Uyy*Uy + Uyz*Uy;
				   double F2z = Uxz*Ux + Uyz*Uy + Uzz*Uz;

				   //dpy dp^2xx

				   double dxdxxU = dHdpNablaLaplacexx(px, py, pz, mass, epsri, epsfi);
				   double dxdyyU = dHdpNablaLaplacexy(py, px, pz, mass, epsri, epsfi);
				   double dxdzzU = dHdpNablaLaplacexy(pz, px, py, mass, epsri, epsfi);

				   double dydxxU = dHdpNablaLaplacexy(px, py, pz, mass, epsri, epsfi);
				   double dydyyU = dHdpNablaLaplacexx(py, px, pz, mass, epsri, epsfi);
				   double dydzzU = dHdpNablaLaplacexy(pz, py, px, mass, epsri, epsfi);

				   double dzdxxU = dHdpNablaLaplacexy(px, pz, py, mass, epsri, epsfi);
				   double dzdyyU = dHdpNablaLaplacexy(py, pz, px, mass, epsri, epsfi);
				   double dzdzzU = dHdpNablaLaplacexx(pz, px, py, mass, epsri, epsfi);


				   double F1x = dxdxxU + dxdyyU + dxdzzU;
				   double F1y = dydxxU + dydyyU + dydzzU;
				   double F1z = dzdxxU + dzdyyU + dzdzzU;

				   double Fx = gamma*(F1x - gamma*F2x) / 6;
				   double Fy = gamma*(F1y - gamma*F2y) / 6;
				   double Fz = gamma*(F1z - gamma*F2z) / 6;


				   //c compute determinant of Id+dt*gamma*nabla2U/3
				   // m11 m12 m13	Uxx Uxy Uxz
				   // m21 m22 m23	Uxy Uyy Uyz
				   // m31 m32 m33	Uxz Uyz Uzz

				   double m11 = Uxx;
				   double m12 = Uxy;
				   double m13 = Uxz;

				   double m21 = Uxy;
				   double m22 = Uyy;
				   double m23 = Uyz;

				   double m31 = Uxz;
				   double m32 = Uyz;
				   double m33 = Uzz;

				   //c compute determinant of B:=Id+dt*gamma*nabla2U/3
				   double t2 = t*t;
				   double gamma2 = gamma*gamma;

				   double t3 = t2*t;
				   double gamma3 = gamma2*gamma;

				   double determinant = 1.0 + (m11 * m22 * m33 - m11 * m23 * m32 - m21 * m12 * m33 + m21 * m13 * m32 + m31 * m12 * m23 - m31 * m13 * m22)*t3 * gamma3 / 27 + (m11 * m22 + m11 * m33 + m22 * m33 - m23 * m32 - m21 * m12 - m31 * m13)*t2 * gamma2 / 9 + (m11 + m22 + m33)*t*gamma / 3;

				   double invB11 = 1 + (m22 * m33 - m23 * m32)*gamma2 * t2 / 9 + (m22 + m33)*gamma*t / 3;
				   invB11 = invB11 / determinant;

				   double invB12 = -((m12 * m33 - m13 * m32))*gamma2 * t2 / 9 - m12 * t*gamma / 3;
				   invB12 = invB12 / determinant;

				   double invB13 = ((m12 * m23 - m13 * m22))*gamma2 * t2 / 9 - m13 * t*gamma / 3;
				   invB13 = invB13 / determinant;

				   double invB21 = -((m21 * m33 - m23 * m31))*gamma2 * t2 / 9 - m21 * t*gamma / 3;
				   invB21 = invB21 / determinant;

				   double invB22 = 1 + (m11 * m33 - m31 * m13)*gamma2 * t2 / 9 + (m11 + m33)*gamma*t / 3;
				   invB22 = invB22 / determinant;

				   double invB23 = -((m11 * m23 - m13 * m21))*gamma2 * t2 / 9 - m23 * t*gamma / 3;
				   invB23 = invB23 / determinant;

				   double invB31 = ((m21 * m32 - m22 * m31))*gamma2 * t2 / 9 - m31 * t*gamma / 3;
				   invB31 = invB31 / determinant;

				   double invB32 = -((m11 * m32 - m12 * m31))*gamma2 * t2 / 9 - m32 * t*gamma / 3;
				   invB32 = invB32 / determinant;

				   double invB33 = 1 + (m11 * m22 - m12 * m21)*gamma2 * t2 / 9 + (m11 + m22)*gamma*t / 3;
				   invB33 = invB33 / determinant;

				   double sigmaMatrix11 = sqrt(invB11) - 1;
				   double sigmaMatrix12 = sqrt(invB12);
				   double sigmaMatrix13 = sqrt(invB13);

				   double sigmaMatrix21 = sqrt(invB21);
				   double sigmaMatrix22 = sqrt(invB22) - 1;
				   double sigmaMatrix23 = sqrt(invB23);

				   double sigmaMatrix31 = sqrt(invB11);
				   double sigmaMatrix32 = sqrt(invB12);
				   double sigmaMatrix33 = sqrt(invB13) - 1;

				   double Sx = sigmaMatrix11*Gx + sigmaMatrix12*Gy + sigmaMatrix13*Gz;
				   double Sy = sigmaMatrix21*Gx + sigmaMatrix22*Gy + sigmaMatrix23*Gz;
				   double Sz = sigmaMatrix31*Gx + sigmaMatrix32*Gy + sigmaMatrix33*Gz;

				   double sq = sqrt(2 * gamma*t);


				   px = px + t2*Fx + sq*Sx;
				   py = py + t2*Fy + sq*Sy;
				   px = px + t2*Fz + sq*Sz;
		}

		break;

		}

	}

}


void	UpdateMomentaStochMidpointExplicit(double t);

double Simulation::computeDeterminant2(double** m){


	return m[1][ 1] * m[2][ 2] - m[1][ 2] * m[2][ 1];


}

double Simulation::computeDeterminant3(double** m){

	return m[1][1] * m[2][2] * m[3][3] - m[1][1] * m[2][3] * m[3][2] - m[1][2] * m[2][1] * m[3][3] + m[1][2] * m[2][3] * m[3][1] + m[1][3] * m[2][1] * m[3][2] - m[1][3] * m[2][2] * m[3][1];


}


void Simulation::ComputeAdditionalTermU3(double t){


		for (int i = 0; i < NumberOfParticles; i++){

			pParticles[i]->updateActiveStatus();
			unsigned int activeStatus = pParticles[i]->getActiveStatus();

			switch (activeStatus){

			case 0: break;

			case 1:{

					   float Gx = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
					   float Gy = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
					   float Gz = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)

					   double px = pParticles[i]->getMomentumX();
					   double py = pParticles[i]->getMomentumX();
					   double pz = pParticles[i]->getMomentumX();

					   double mass = pParticles[i]->getMass();
					   //double epsri = pParticles[i]->getEpsr();
					   //double epsfi = pParticles[i]->getEpsf();

					   double Ux = px / mass;// dHdp(px, py, pz, mass, epsri, epsfi);
					   double Uy = py / mass;// dHdp(py, px, pz, mass, epsri, epsfi);
					   double Uz = pz / mass;// dHdp(pz, py, px, mass, epsri, epsfi);


					   double Uxx = 1 / mass;// dHdp2xx(px, py, pz, mass, epsri, epsfi);
					   double Uxy = 0.0;// dHdp2xy(px, py, pz, mass, epsri, epsfi);
					   double Uxz = 0.0;// dHdp2xy(px, pz, py, mass, epsri, epsfi);
					   double Uyy = 1 / mass;// dHdp2xx(py, px, pz, mass, epsri, epsfi);
					   double Uzz = 1.0 / mass;// dHdp2xx(pz, py, px, mass, epsri, epsfi);
					   double Uyz = 0.0;// dHdp2xy(py, pz, px, mass, epsri, epsfi);


					   double F2x = Uxx*Ux;
					   double F2y = Uyy*Uy;
					   double F2z = Uzz*Uz;

					   //dpy dp^2xx in std case is 0



					   double F1x = 0.0;
					   double F1y = 0.0;
					   double F1z = 0.0;

					   double Fx = gamma*(F1x - gamma*F2x) / 6;
					   double Fy = gamma*(F1y - gamma*F2y) / 6;
					   double Fz = gamma*(F1z - gamma*F2z) / 6;


					   //c compute determinant of Id+dt*gamma*nabla2U/3

					   double t2 = t*t;

					   double invB = 1 / (1 + t*gamma / (mass * 3));

					   double sigmaMatrix11 = sqrt(invB) - 1;
					   double sigmaMatrix22 = sqrt(invB) - 1;
					   double sigmaMatrix33 = sqrt(invB) - 1;

					   double Sx = sigmaMatrix11*Gx;
					   double Sy = sigmaMatrix22*Gy;
					   double Sz = sigmaMatrix33*Gz;

					   double sq = sqrt(2 * gamma*t);


					   px = px + t2*Fx + sq*Sx;
					   py = py + t2*Fy + sq*Sy;
					   px = px + t2*Fz + sq*Sz;
			}
				break;

			case 2:{

					   float Gx = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
					   float Gy = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)
					   float Gz = ComputeBoxMuller(0, 1); // G^n/sqrt(beta)

					   double px = pParticles[i]->getMomentumX();
					   double py = pParticles[i]->getMomentumX();
					   double pz = pParticles[i]->getMomentumX();

					   double mass = pParticles[i]->getMass();
					   double epsri = pParticles[i]->getEpsr();
					   double epsfi = pParticles[i]->getEpsf();

					   double Ux = dHdp(px, py, pz, mass, epsri, epsfi);
					   double Uy = dHdp(py, px, pz, mass, epsri, epsfi);
					   double Uz = dHdp(pz, py, px, mass, epsri, epsfi);


					   double Uxx = dHdp2xx(px, py, pz, mass, epsri, epsfi);
					   double Uxy = dHdp2xy(px, py, pz, mass, epsri, epsfi);
					   double Uxz = dHdp2xy(px, pz, py, mass, epsri, epsfi);
					   double Uyy = dHdp2xx(py, px, pz, mass, epsri, epsfi);
					   double Uzz = dHdp2xx(pz, py, px, mass, epsri, epsfi);
					   double Uyz = dHdp2xy(py, pz, px, mass, epsri, epsfi);


					   double F2x = Uxx*Ux + Uxy*Uy + Uxz*Uz;
					   double F2y = Uxy*Ux + Uyy*Uy + Uyz*Uy;
					   double F2z = Uxz*Ux + Uyz*Uy + Uzz*Uz;

					   //dpy dp^2xx

					   double dxdxxU = dHdpNablaLaplacexx(px, py, pz, mass, epsri, epsfi);
					   double dxdyyU = dHdpNablaLaplacexy(py, px, pz, mass, epsri, epsfi);
					   double dxdzzU = dHdpNablaLaplacexy(pz, px, py, mass, epsri, epsfi);

					   double dydxxU = dHdpNablaLaplacexy(px, py, pz, mass, epsri, epsfi);
					   double dydyyU = dHdpNablaLaplacexx(py, px, pz, mass, epsri, epsfi);
					   double dydzzU = dHdpNablaLaplacexy(pz, py, px, mass, epsri, epsfi);

					   double dzdxxU = dHdpNablaLaplacexy(px, pz, py, mass, epsri, epsfi);
					   double dzdyyU = dHdpNablaLaplacexy(py, pz, px, mass, epsri, epsfi);
					   double dzdzzU = dHdpNablaLaplacexx(pz, px, py, mass, epsri, epsfi);


					   double F1x = dxdxxU + dxdyyU + dxdzzU;
					   double F1y = dydxxU + dydyyU + dydzzU;
					   double F1z = dzdxxU + dzdyyU + dzdzzU;

					   double Fx = gamma*(F1x - gamma*F2x) / 6;
					   double Fy = gamma*(F1y - gamma*F2y) / 6;
					   double Fz = gamma*(F1z - gamma*F2z) / 6;


					   //c compute determinant of Id+dt*gamma*nabla2U/3
					   // m11 m12 m13	Uxx Uxy Uxz
					   // m21 m22 m23	Uxy Uyy Uyz
					   // m31 m32 m33	Uxz Uyz Uzz

					   double m11 = Uxx;
					   double m12 = Uxy;
					   double m13 = Uxz;

					   double m21 = Uxy;
					   double m22 = Uyy;
					   double m23 = Uyz;

					   double m31 = Uxz;
					   double m32 = Uyz;
					   double m33 = Uzz;

					   //c compute determinant of B:=Id+dt*gamma*nabla2U/3
					   double t2 = t*t;
					   double gamma2 = gamma*gamma;

					   double t3 = t2*t;
					   double gamma3 = gamma2*gamma;

					   double determinant = 1.0 + (m11 * m22 * m33 - m11 * m23 * m32 - m21 * m12 * m33 + m21 * m13 * m32 + m31 * m12 * m23 - m31 * m13 * m22)*t3 * gamma3 / 27 + (m11 * m22 + m11 * m33 + m22 * m33 - m23 * m32 - m21 * m12 - m31 * m13)*t2 * gamma2 / 9 + (m11 + m22 + m33)*t*gamma / 3;

					   double invB11 = 1 + (m22 * m33 - m23 * m32)*gamma2 * t2 / 9 + (m22 + m33)*gamma*t / 3;
					   invB11 = invB11 / determinant;

					   double invB12 = -((m12 * m33 - m13 * m32))*gamma2 * t2 / 9 - m12 * t*gamma / 3;
					   invB12 = invB12 / determinant;

					   double invB13 = ((m12 * m23 - m13 * m22))*gamma2 * t2 / 9 - m13 * t*gamma / 3;
					   invB13 = invB13 / determinant;

					   double invB21 = -((m21 * m33 - m23 * m31))*gamma2 * t2 / 9 - m21 * t*gamma / 3;
					   invB21 = invB21 / determinant;

					   double invB22 = 1 + (m11 * m33 - m31 * m13)*gamma2 * t2 / 9 + (m11 + m33)*gamma*t / 3;
					   invB22 = invB22 / determinant;

					   double invB23 = -((m11 * m23 - m13 * m21))*gamma2 * t2 / 9 - m23 * t*gamma / 3;
					   invB23 = invB23 / determinant;

					   double invB31 = ((m21 * m32 - m22 * m31))*gamma2 * t2 / 9 - m31 * t*gamma / 3;
					   invB31 = invB31 / determinant;

					   double invB32 = -((m11 * m32 - m12 * m31))*gamma2 * t2 / 9 - m32 * t*gamma / 3;
					   invB32 = invB32 / determinant;

					   double invB33 = 1 + (m11 * m22 - m12 * m21)*gamma2 * t2 / 9 + (m11 + m22)*gamma*t / 3;
					   invB33 = invB33 / determinant;

					   double sigmaMatrix11 = sqrt(invB11) - 1;
					   double sigmaMatrix12 = sqrt(invB12);
					   double sigmaMatrix13 = sqrt(invB13);

					   double sigmaMatrix21 = sqrt(invB21);
					   double sigmaMatrix22 = sqrt(invB22) - 1;
					   double sigmaMatrix23 = sqrt(invB23);

					   double sigmaMatrix31 = sqrt(invB11);
					   double sigmaMatrix32 = sqrt(invB12);
					   double sigmaMatrix33 = sqrt(invB13) - 1;

					   double Sx = sigmaMatrix11*Gx + sigmaMatrix12*Gy + sigmaMatrix13*Gz;
					   double Sy = sigmaMatrix21*Gx + sigmaMatrix22*Gy + sigmaMatrix23*Gz;
					   double Sz = sigmaMatrix31*Gx + sigmaMatrix32*Gy + sigmaMatrix33*Gz;

					   double sq = sqrt(2 * gamma*t);


					   px = px + t2*Fx + sq*Sx;
					   py = py + t2*Fy + sq*Sy;
					   px = px + t2*Fz + sq*Sz;
			}

				break;

			}

		}

	}


	void Simulation::UpdateMomentaStochExtendedMidpoint(double t){


#ifdef CHECK_SPLITTING_SCHEME
		if (current_n == 0) std::cout << "UpdateMomentaStochExtendedMidpoint " << t << std::endl;
#endif

		for (int i = 0; i < NumberOfParticles; i++){

			float Gx = ComputeBoxMuller(0, 1);
			float Gy = ComputeBoxMuller(0, 1);
			float Gz = ComputeBoxMuller(0, 1);

			double c2x = sigma*Gx*sqrt(t);
			double c2y = sigma*Gy*sqrt(t);
			double c2z = sigma*Gz*sqrt(t);

			double px = pParticles[i]->getMomentumX();
			double py = pParticles[i]->getMomentumY();
			double pz = pParticles[i]->getMomentumZ();

			double epsri = pParticles[i]->getEpsr();
			double epsfi = pParticles[i]->getEpsf();
			double mass = pParticles[i]->getMass();

			double px_tmp = dHdp(px, py, pz, mass, epsri, epsfi);
			double py_tmp = dHdp(py, px, pz, mass, epsri, epsfi);
			double pz_tmp = dHdp(pz, px, py, mass, epsri, epsfi);

			double px_tmp_old = px_tmp;
			double py_tmp_old = py_tmp;
			double pz_tmp_old = pz_tmp;

			unsigned int iteration_counter = 0;
			do{

				px_tmp_old = px_tmp;
				py_tmp_old = py_tmp;
				pz_tmp_old = pz_tmp;

				px_tmp = px - gamma*dHdp(MidPoint(px_tmp, px), py, pz, mass, epsri, epsfi)*t + c2x;
				py_tmp = py - gamma*dHdp(MidPoint(py_tmp, py), px, pz, mass, epsri, epsfi)*t + c2y;
				pz_tmp = pz - gamma*dHdp(MidPoint(pz_tmp, pz), px, py, mass, epsri, epsfi)*t + c2z;

				iteration_counter++;

			} while (fabs(px_tmp - px_tmp_old) > fix_point_tolerance || fabs(py_tmp - py_tmp_old) > fix_point_tolerance || fabs(pz_tmp - pz_tmp_old) > fix_point_tolerance);

			//std::cout<<"Fixpoint iteration:"<< iteration_counter<<std::endl;
			pParticles[i]->setMomentum(px_tmp, py_tmp, pz_tmp);
		}


	};

	void Simulation::UpdateMomentaStochOrder2(double t){


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateMomentaStochOrder2 "<<t << std::endl;
#endif

	for (int i = 0; i < NumberOfParticles; i++){

		float Gx = ComputeBoxMuller(0, 1);
		float Gy = ComputeBoxMuller(0, 1);
		float Gz = ComputeBoxMuller(0, 1);

		double c2x = sigma*Gx*sqrt(t);
		double c2y = sigma*Gy*sqrt(t);
		double c2z = sigma*Gz*sqrt(t);

		double px = pParticles[i]->getMomentumX();
		double py = pParticles[i]->getMomentumY();
		double pz = pParticles[i]->getMomentumZ();

		double epsri = pParticles[i]->getEpsr();
		double epsfi = pParticles[i]->getEpsf();
		double mass = pParticles[i]->getMass();

		double px_tmp = dHdp(px, py, pz, mass, epsri, epsfi);
		double py_tmp = dHdp(py, px, pz, mass, epsri, epsfi);
		double pz_tmp = dHdp(pz, px, py, mass, epsri, epsfi);

		double px_tmp_old = px_tmp;
		double py_tmp_old = py_tmp;
		double pz_tmp_old = pz_tmp;

		unsigned int iteration_counter = 0;
		do{

			px_tmp_old = px_tmp;
			py_tmp_old = py_tmp;
			pz_tmp_old = pz_tmp;

			px_tmp = px - gamma*dHdp(MidPoint(px_tmp, px), py, pz, mass, epsri, epsfi)*t + c2x;
			py_tmp = py - gamma*dHdp(MidPoint(py_tmp, py), px, pz, mass, epsri, epsfi)*t + c2y;
			pz_tmp = pz - gamma*dHdp(MidPoint(pz_tmp, pz), px, py, mass, epsri, epsfi)*t + c2z;

			iteration_counter++;

		} while (fabs(px_tmp - px_tmp_old) > fix_point_tolerance || fabs(py_tmp - py_tmp_old) > fix_point_tolerance || fabs(pz_tmp - pz_tmp_old) > fix_point_tolerance);

		//std::cout<<"Fixpoint iteration:"<< iteration_counter<<std::endl;
		pParticles[i]->setMomentum(px_tmp, py_tmp, pz_tmp);
	}


};



void Simulation::UpdateForcesParticles() {

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticles" << std::endl;
#endif


	//potentialEnergy = 0.0;
	DimerDistance = 0.0;
	DimerPotential = 0.0;

	interactionsCounter = 0.0;

	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->setF(0.0, 0.0 ,0.0);

	}

	//std::cout<<NumberOfParticles<<std::endl;
	//int nPairs = 0;

	/*if (current_n == 1000)
		int tmp = 123;*/

	for (int i = 0; i < NumberOfParticles; i++){

		// solvent type 0, dimer type 1
		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (int j = 0; j < NumberOfParticles; j++){

			if (i == j) continue;

			interactionsCounter++;

			unsigned int typeParticle_j = pParticles[j]->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			if (typeParticle_i + typeParticle_j == 2) {

				DimerDistance = DimerDistance + r;
				DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + 0.5*InteractionPotential(r, typeParticle_i, typeParticle_j);



			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);

			//if (act != 0) {
			//	double tmp = InteractionPotential(r, typeParticle_i, typeParticle_j);
			//	nPairs++;
			//}

		}
	}

	/*potentialEnergy /= 2;
	if (nPairs >= 1)
		int tmp = 123;*/

	//std::cout << current_n << "\tnPairs\t" << nPairs << std::endl;

	NumberOfInteractions.addSample(interactionsCounter);

};


void Simulation::UpdateForcesParticlesNewton() {


#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesNewton" << std::endl;
#endif

	potentialEnergy = 0.0;
	DimerDistance = 0.0;
	DimerPotential = 0.0;
	DimerSolventPotential = 0.0;
	SolventSolventPotential = 0.0;

	interactionsCounter = 0.0;


	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->setF(0.0, 0.0, 0.0);

	}


	for (int i = 0; i < NumberOfParticles; i++){

#ifdef ONE_PARTICLE_ONLY


		double qijx = pParticles[i]->getPositionX();
		double qijy = pParticles[i]->getPositionY();
		double qijz = pParticles[i]->getPositionZ();

		if (PBC == 1){

			qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
			qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
			qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
		}

		double r =  sqrt(qijx*qijx + qijy*qijy+ qijz*qijz);

		double dx = qijx;// / r;
		double dy = qijy;// / r;
		double dz = qijz;// / r;

		//double act = r;
		// harmonic oscillator V=0.5kr^2 -> d/dx V =kx

		potentialEnergy = potentialEnergy + 0.5*r*r;
		DimerPotential = DimerPotential + 0.5*r*r;

		pParticles[i]->setF(dx, dy, dz);


#endif
#ifndef ONE_PARTICLE_ONLY

		// solvent type 0, dimer type 1
		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (int j = i + 1; j < NumberOfParticles; j++){

			unsigned int typeParticle_j = pParticles[j]->getParticleType();


			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			if (typeParticle_i + typeParticle_j == 2) 			 {

				DimerDistance = DimerDistance + r;
				DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}
			if (typeParticle_i + typeParticle_j == 1) 			 {

				DimerSolventPotential = DimerSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

				SolventSolventPotential = SolventSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}


			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);

			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
			pParticles[j]->setF(pParticles[j]->getFX() - dx*act, pParticles[j]->getFY() - dy*act, pParticles[j]->getFZ() - dz*act);

		}

#endif
	}

	//std::cout << abs(potentialEnergy - DimerPotential - DimerSolventPotential - SolventSolventPotential) << std::endl;

	NumberOfInteractions.addSample(interactionsCounter);

};


void Simulation::UpdateForcesParticlesAdaptivelyNewton() {

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptively" << std::endl;
#endif

	//potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!
	//DimerDistance = 0.0;
	//DimerPotential = 0.0;

	//int numberOfActiveParticles = Particle::howManyActive;
	//int numberOfRestrainedParticles = Particle::howManyRestrained;

	interactionsCounter = 0.0;

	for (int i = 0; i < NumberOfParticles; i++) {

		pParticles[i]->setF(0.0, 0.0, 0.0);

	}

	//----------- dimer

	for (unsigned int i = 0; i<NumberOfDimerParticles; i++) {


		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = i + 1; j<NumberOfDimerParticles; j++) {


			interactionsCounter++;


			unsigned int typeParticle_j = pParticles[j]->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			/*	if (typeParticle_i + typeParticle_j == 2) {

			DimerDistance = DimerDistance + r;
			DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

			SolventSolventPotential = SolventSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
			*/
			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
			pParticles[j]->setF(pParticles[j]->getFX() - dx*act, pParticles[j]->getFY() - dy*act, pParticles[j]->getFZ() - dz*act);

		}

	}



	//-----------------------

	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);

		//if (pParticles[i]->getActiveStatus() == 0) continue;

		//std::cout << neighborIndexI->size() << std::endl;
		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			//int indexOfParticleJ = neighborIndexI->getIndex(neighborIndexI->getObject(j));


			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));

			if (particleJ->getParticleIndex() > i) continue;

			//std::cout << indexOfParticleJ << " \t" << particleJ->getParticleIndex()<<std::endl;

			interactionsCounter++;

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			if (typeParticle_j == 1 && typeParticle_i == 1) continue;

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

		/*	if (typeParticle_i + typeParticle_j == 2) {

				DimerDistance = DimerDistance + r;
				DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

				SolventSolventPotential = SolventSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
*/
			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
			particleJ->setF(particleJ->getFX() - dx*act, particleJ->getFY() - dy*act, particleJ->getFZ() - dz*act);

		}

	}

	NumberOfInteractions.addSample(interactionsCounter);


};




void Simulation::UpdateForcesParticlesAdaptively() {

	//unsigned int nPairs=0;

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptively" << std::endl;
#endif

	potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!
	DimerDistance = 0.0;
	DimerPotential = 0.0;

	//int numberOfActiveParticles = Particle::howManyActive;
	//int numberOfRestrainedParticles = Particle::howManyRestrained;

	interactionsCounter = 0.0;

	for (int i = 0; i < NumberOfParticles ; i++) {

		pParticles[i]->setF(0.0, 0.0, 0.0);

	}

	/*if (current_n == 8)
		int ttt = 123;*/

	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);

		//if (pParticles[i]->getActiveStatus() == 0) continue;

		//std::cout << neighborIndexI->size() << std::endl;
		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			interactionsCounter++;

			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			if (typeParticle_i + typeParticle_j == 2) {

				DimerDistance = DimerDistance + r;
				DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);

			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
			//
			/*if (act != 0) {
				double tmp = InteractionPotential(r, typeParticle_i, typeParticle_j);
				nPairs++;
			}*/

		}

	}

	// we take everything twice
	//potentialEnergy /= 2;

	/*if (nPairs >= 5)
		int tmp = 123;*/
	/*if (current_n == 8)
		grid->print();*/

	//std::cout << current_n << "\tnPairs\t" << nPairs << std::endl;

	NumberOfInteractions.addSample(interactionsCounter);

};

void Simulation::UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT() {

	//unsigned int nPairs=0;

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptivelyARPS_SUBSTRACT" << std::endl;
#endif

	//std::cout << "potetial energy " << potentialEnergy << std::endl;
	//while (1);

	//potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!

	//no need to compute these two in substract since dimer is always active
	//DimerDistance = 0.0;
	//DimerPotential = 0.0;

	int numberOfActiveParticles = Particle::howManyActive;
	int numberOfRestrainedParticles = Particle::howManyRestrained;

	interactionsCounter = 0.0;

	/*if (current_n == 4858)
	int ttt = 123;*/

	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		//if particle is restrained -> continue
		if (pParticles[i]->getActiveStatus() == 0)
			continue;

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);

		//std::cout << neighborIndexI->size() << std::endl;

		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));


			interactionsCounter++;

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			/*double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
			double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
			double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());*/

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);

				//qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
				//qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
				//qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);

			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
			//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			//double dxOld = qijxOld / rOld;
			//double dyOld = qijyOld / rOld;
			//double dzOld = qijzOld / rOld;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);
			//double actOld =  InteractionForce(rOld, typeParticle_i, typeParticle_j);

			/*double interactionPot = InteractionPotential(r, typeParticle_i, typeParticle_j);

			if (particleJ->getActiveStatus() > 0){
				potentialEnergy = potentialEnergy - 0.5*interactionPot;
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

				SolventSolventPotential = SolventSolventPotential - interactionPot;
			}*/

						// substract force for the active particle
			pParticles[i]->setF(pParticles[i]->getFX() - dx*act, pParticles[i]->getFY() - dy*act, pParticles[i]->getFZ() - dz*act);
			if (particleJ->getActiveStatus() == 0){
				// if the other one is restrained
				particleJ->setF(particleJ->getFX() + dx*act, particleJ->getFY() + dy*act, particleJ->getFZ() + dz*act);

				//potentialEnergy = potentialEnergy - interactionPot;
			}

		}

	}

	// we take everything twice
	//potentialEnergy /= 2;


	//NumberOfInteractions.addSample(interactionsCounter);

};


void Simulation::UpdateForcesParticlesAdaptivelyARPS_ADD() {

	//unsigned int nPairs=0;

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptivelyARPS_ADD" << std::endl;
#endif

	//potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!
	//DimerDistance = 0.0;
	//DimerPotential = 0.0;

	int numberOfActiveParticles = Particle::howManyActive;
	int numberOfRestrainedParticles = Particle::howManyRestrained;

	//interactionsCounter = 0.0;

	/*if (current_n == 4858)
	int ttt = 123;*/

	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		if (pParticles[i]->getActiveStatus() == 0)
			continue;

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);

		//std::cout << neighborIndexI->size() << std::endl;

		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));



			interactionsCounter++;

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			//double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
			//double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
			//double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);

				//qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
				//qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
				//qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);

			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
			//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			/*double dxOld = qijxOld / rOld;
			double dyOld = qijyOld / rOld;
			double dzOld = qijzOld / rOld;*/

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			double interactionPot = InteractionPotential(r, typeParticle_i, typeParticle_j);

			/*if (particleJ->getActiveStatus() > 0){
				potentialEnergy = potentialEnergy + 0.5*interactionPot;
			}


			if (typeParticle_i + typeParticle_j == 0) 			 {

				SolventSolventPotential = SolventSolventPotential + interactionPot;
			}*/

			//add force for the active particle
			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);

			if (particleJ->getActiveStatus() == 0)
				// if the other one is restrained
				particleJ->setF(particleJ->getFX() - dx*act, particleJ->getFY() - dy*act, particleJ->getFZ() - dz*act);

		//	potentialEnergy = potentialEnergy + interactionPot;


		}

	}

	// we take everything twice
	//potentialEnergy /= 2;


	NumberOfInteractions.addSample(interactionsCounter);

};



void Simulation::UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT() {

	//unsigned int nPairs=0;

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptivelyNewtonARPS_SUBSTRACT" << std::endl;
#endif

	//std::cout << "potetial energy " << potentialEnergy << std::endl;
	//while (1);

	//potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!

	//no need to compute these two in substract since dimer is always active
	//DimerDistance = 0.0;
	//DimerPotential = 0.0;

	int numberOfActiveParticles = Particle::howManyActive;
	int numberOfRestrainedParticles = Particle::howManyRestrained;

	interactionsCounter = 0.0;

	/*if (current_n == 4858)
	int ttt = 123;*/

	//----------- dimer

	for (unsigned int i = 0; i<NumberOfDimerParticles; i++) {


		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = i + 1; j<NumberOfDimerParticles; j++) {


			interactionsCounter++;


			unsigned int typeParticle_j = pParticles[j]->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);
			potentialEnergy = potentialEnergy - InteractionPotential(r, typeParticle_i, typeParticle_j);
			/*	if (typeParticle_i + typeParticle_j == 2) {

			DimerDistance = DimerDistance + r;
			DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

			SolventSolventPotential = SolventSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
			*/
			pParticles[i]->setF(pParticles[i]->getFX() - dx*act, pParticles[i]->getFY() - dy*act, pParticles[i]->getFZ() - dz*act);
			pParticles[j]->setF(pParticles[j]->getFX() + dx*act, pParticles[j]->getFY() + dy*act, pParticles[j]->getFZ() + dz*act);

		}

	}


	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		//if particle is restrained -> continue
		if (pParticles[i]->getActiveStatus() == 0)
			continue;



		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);

		//std::cout << neighborIndexI->size() << std::endl;

		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));

			int realIndexJ = particleJ->getParticleIndex();

			if ((particleJ->getActiveStatus() > 0) && realIndexJ >= i)
				continue;

			interactionsCounter++;

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			//skip dimer dimer interaction
			if (typeParticle_i == 1 && typeParticle_j == 1) continue;

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			/*double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
			double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
			double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());*/

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);

				//qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
				//qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
				//qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);

			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
			//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			//double dxOld = qijxOld / rOld;
			//double dyOld = qijyOld / rOld;
			//double dzOld = qijzOld / rOld;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);
			potentialEnergy = potentialEnergy - InteractionPotential(r, typeParticle_i, typeParticle_j);
			//double actOld =  InteractionForce(rOld, typeParticle_i, typeParticle_j);

			/*double interactionPot = InteractionPotential(r, typeParticle_i, typeParticle_j);

			if (particleJ->getActiveStatus() > 0){
			potentialEnergy = potentialEnergy - 0.5*interactionPot;
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

			SolventSolventPotential = SolventSolventPotential - interactionPot;
			}*/

			// subtract force for the active particle
			pParticles[i]->setF(pParticles[i]->getFX() - dx*act, pParticles[i]->getFY() - dy*act, pParticles[i]->getFZ() - dz*act);
			//if (particleJ->getActiveStatus() == 0){
				// if the other one is restrained
				particleJ->setF(particleJ->getFX() + dx*act, particleJ->getFY() + dy*act, particleJ->getFZ() + dz*act);


			//}

		}

	}

	// we take everything twice
	//potentialEnergy /= 2;


	//NumberOfInteractions.addSample(interactionsCounter);

};


void Simulation::UpdateForcesParticlesAdaptivelyNewtonARPS_ADD() {

	//unsigned int nPairs=0;

#ifdef CHECK_SPLITTING_SCHEME
	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptivelyNewtonARPS_ADD" << std::endl;
#endif

	//potentialEnergy = 0.0;
	// set dimer functions zero because dimer is always active!!
	//DimerDistance = 0.0;
	//DimerPotential = 0.0;

	int numberOfActiveParticles = Particle::howManyActive;
	int numberOfRestrainedParticles = Particle::howManyRestrained;

	//interactionsCounter = 0.0;

	/*if (current_n == 4858)
	int ttt = 123;*/

	//----------- dimer

	for (unsigned int i = 0; i<NumberOfDimerParticles; i++) {


		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = i + 1; j<NumberOfDimerParticles; j++) {


			interactionsCounter++;


			unsigned int typeParticle_j = pParticles[j]->getParticleType();

			double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);
			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
			/*	if (typeParticle_i + typeParticle_j == 2) {

			DimerDistance = DimerDistance + r;
			DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			if (typeParticle_i + typeParticle_j == 0) 			 {

			SolventSolventPotential = SolventSolventPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
			}

			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
			*/
			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
			pParticles[j]->setF(pParticles[j]->getFX() - dx*act, pParticles[j]->getFY() - dy*act, pParticles[j]->getFZ() - dz*act);

		}

	}


	for (unsigned int i = 0; i<NumberOfParticles; i++) {

		if (pParticles[i]->getActiveStatus() == 0)
			continue;

		SBCContainerIndex<Particle*> const* neighborIndexI = grid->getNeighborIndex(i);



		//std::cout << neighborIndexI->size() << std::endl;

		unsigned int typeParticle_i = pParticles[i]->getParticleType();

		for (unsigned int j = 0; j<neighborIndexI->size(); j++) {

			Particle* particleJ = static_cast<Particle*>(neighborIndexI->getObject(j));

			int realIndexJ = particleJ->getParticleIndex();

			if ((particleJ->getActiveStatus() >0) && realIndexJ >= i)
				continue;


			interactionsCounter++;

			//	std::cout << "Neighbor pair: " << i << " " << static_cast<Particle*>(neighborIndexI->getObject(j)) << std::endl;
			//nPairs++;
			//ForcesCalculation(pParticles[i], particleJ);

			unsigned int typeParticle_j = particleJ->getParticleType();

			if (typeParticle_i == 1 && typeParticle_j == 1) continue;

			double qijx = (pParticles[i]->getPositionX() - particleJ->getPositionX());
			double qijy = (pParticles[i]->getPositionY() - particleJ->getPositionY());
			double qijz = (pParticles[i]->getPositionZ() - particleJ->getPositionZ());

			//double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
			//double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
			//double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());

			if (PBC == 1){
				qijx = qijx - BoxLength*roundMy(qijx / BoxLength);//   minus(qij,multiplyByScalar( roundMy( divideByScalar(qij,BoxLength)) , BoxLength));
				qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
				qijz = qijz - BoxLength*roundMy(qijz / BoxLength);

				//qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
				//qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
				//qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);

			}

			double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
			//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);

			double dx = qijx / r;
			double dy = qijy / r;
			double dz = qijz / r;

			/*double dxOld = qijxOld / rOld;
			double dyOld = qijyOld / rOld;
			double dzOld = qijzOld / rOld;*/

			double act = InteractionForce(r, typeParticle_i, typeParticle_j);

			double interactionPot = InteractionPotential(r, typeParticle_i, typeParticle_j);
			potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);
			/*if (particleJ->getActiveStatus() > 0){
			potentialEnergy = potentialEnergy + 0.5*interactionPot;
			}


			if (typeParticle_i + typeParticle_j == 0) 			 {

			SolventSolventPotential = SolventSolventPotential + interactionPot;
			}*/

			//add force for the active particle
			pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);

		//	if (particleJ->getActiveStatus() == 0)
				// if the other one is restrained
				particleJ->setF(particleJ->getFX() - dx*act, particleJ->getFY() - dy*act, particleJ->getFZ() - dz*act);

			//	potentialEnergy = potentialEnergy + interactionPot;


		}

	}

	// we take everything twice
	//potentialEnergy /= 2;


	NumberOfInteractions.addSample(interactionsCounter);

};






//
//
//void Simulation::UpdateForcesParticlesAdaptively() {
//
//
//
//#ifdef CHECK_SPLITTING_SCHEME
//	if (current_n == 0) std::cout << "UpdateForcesParticlesAdaptively" << std::endl;
//#endif
//
//	// set dimer functions zero because dimer is always active!!
//	DimerDistance = 0.0;
//	DimerPotential = 0.0;
//
//
//	int numberOfActiveParticles = Particle::howManyActive;
//	int numberOfRestrainedParticles = Particle::howManyRestrained;
//
//
//	interactionsCounter = 0.0;
//
//	//for (int i = 0; i < numberOfRestrainedParticles; i++) {
//
//	//	//std::cout << Particle::howManyRestrained << std::endl;
//	//	std::cout<< indexArrayRestrainedParticles[i]<<std::endl;
//
//	//
//	//}
//
//	//while (1);
//
//	for (int i = 0; i < numberOfActiveParticles; i++) {
//
//		int indexI = indexArrayActiveParticles[i];
//
//		pParticles[indexI]->setF(0.0, 0.0, 0.0);
//
//	}
//
//
//	//double fActiveX = 0;
//	//double fActiveY = 0;
//
//	// for i in active particles- i=1:K
//	for (int i = 0; i < numberOfActiveParticles; i++) {
//
//		int indexI = indexArrayActiveParticles[i];
//
//
//
//		//for (int j = i + 1; j<numberOfActiveParticles; j++){
//		for (int indexJ = 0; indexJ < NumberOfParticles; indexJ++){
//
//
//			if (indexI == indexJ) continue;
//
//
//			//	int indexJ = indexArrayActiveParticles[j];
//
//			//ForcesCalculation(indexI, indexJ);
//			ForcesCalculationActiveParticles(indexI, indexJ);
//
//			//	std::cout << "I="<<indexI<<"J="<<indexJ << std::endl;
//		}
//
//	}
//
//	//---------- restrained particles loop -------------------------------
//	for (int i = 0; i < numberOfRestrainedParticles; i++) {
//
//		int indexI = indexArrayRestrainedParticles[i];
//
//		for (int j = 0; j < numberOfActiveParticles; j++){
//
//			int indexJ = indexArrayActiveParticles[j];
//
//			if (indexI == indexJ)continue;
//
//			//	ForcesCalculation(indexI, indexJ);
//			ForcesCalculationRestrainedParticles(indexI, indexJ);
//		}
//
//	}
//
//	//for (int i = 0; i < numberOfRestrainedParticles; i++) {
//
//	//	int indexI = indexArrayRestrainedParticles[i];
//
//	//	//-------- active particles loop
//
//	//	for (int j = i + 1; j < numberOfRestrainedParticles; j++){
//
//	//		int indexJ = indexArrayRestrainedParticles[j];
//
//	//		ForcesCalculation(indexI, indexJ);
//
//	//	}
//	//}
//
//
//	NumberOfInteractions.addSample(interactionsCounter);
//	//std::cout << "\n " << interactionsCounter << std::endl;
//
//	/*	for (int i=0; i<NumberOfParticles; i++)
//	{
//	std::cout<< "F(x,y,z)="<<pParticles[i]->getF().getX()<< "\t"<<pParticles[i]->getF().getY()<<"\t"<<pParticles[i]->getF().getZ()<<std::endl;
//
//	}*/
//
//	//while(1);
//
//};

//
//void Simulation::ForcesCalculationActiveParticles(int i, int j){
//
//
//	unsigned int typeParticle_i = pParticles[i]->getParticleType();
//	unsigned int typeParticle_j = pParticles[j]->getParticleType();
//
//	double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
//	double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
//	double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());
//
//	//double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
//	//double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	//double dxOld = qijxOld / rOld;
//	//double dyOld = qijyOld / rOld;
//
//	double act = 0.0;
//	//double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	//potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);// -InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
//
//
//
//	//return 1.0;
//};
//
//void Simulation::ForcesCalculationActiveParticles(Particle* particleI, Particle* particleJ){
//
//
//	unsigned int typeParticle_i = particleI->getParticleType();
//	unsigned int typeParticle_j = particleJ->getParticleType();
//
//	double qijx = (particleI->getPositionX() - particleJ->getPositionX());
//	double qijy = (particleI->getPositionY() - particleJ->getPositionY());
//	double qijz = (particleI->getPositionZ() - particleJ->getPositionZ());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	//double dxOld = qijxOld / rOld;
//	//double dyOld = qijyOld / rOld;
//
//	double act = 0.0;
//	//double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);// -InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	particleI->setF(particleI->getFX() + dx*act, particleI->getFY() + dy*act, particleI->getFZ() + dz*act);
//
//
//};
//
//
//void Simulation::ForcesCalculationActiveParticlesNewton(int i, int j){
//
//
//	unsigned int typeParticle_i = pParticles[i]->getParticleType();
//	unsigned int typeParticle_j = pParticles[j]->getParticleType();
//
//	double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
//	double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
//	double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());
//
//	//double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
//	//double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//		//qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
//		//qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	//double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	//double dxOld = qijxOld / rOld;
//	//double dyOld = qijyOld / rOld;
//
//	double act = 0.0;
//	//double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j);// -InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	//pParticles[i]->setF(pParticles[i]->getFX() - dxOld*actOld + dx*act, pParticles[i]->getFY() - dyOld*actOld + dy*act);
//	//pParticles[j]->setF(pParticles[j]->getFX() - (-dxOld*actOld + dx*act), pParticles[j]->getFY() - (-dyOld*actOld + dy*act));
//
//	pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act, pParticles[i]->getFZ() + dz*act);
//	pParticles[j]->setF(pParticles[j]->getFX() - dx*act, pParticles[j]->getFY() - dy*act, pParticles[j]->getFZ() - dz*act);
//
//
//	//return 1.0;
//};
//
//
//
//
//void Simulation::ForcesCalculation(int i, int j){
//
//
//	unsigned int typeParticle_i = pParticles[i]->getParticleType();
//	unsigned int typeParticle_j = pParticles[j]->getParticleType();
//
//	double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
//	double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
//	double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());
//
//	double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
//	double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
//	double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//		qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
//		qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
//		qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	double dxOld = qijxOld / rOld;
//	double dyOld = qijyOld / rOld;
//	double dzOld = qijzOld / rOld;
//
//	double act = 0.0;
//	double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j) - InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	pParticles[i]->setF(pParticles[i]->getFX() - dxOld*actOld + dx*act, pParticles[i]->getFY() - dyOld*actOld + dy*act, pParticles[i]->getFZ() - dzOld*actOld + dz*act);
//	pParticles[j]->setF(pParticles[j]->getFX() - (-dxOld*actOld + dx*act), pParticles[j]->getFY() - (-dyOld*actOld + dy*act), pParticles[j]->getFZ() - (-dzOld*actOld + dz*act));
//
//	//pParticles[i]->setF(pParticles[i]->getFX() + dx*act, pParticles[i]->getFY() + dy*act);
//	//pParticles[j]->setF(pParticles[j]->getFX() - dx*act, pParticles[j]->getFY() - dy*act);
//
//
//	//return 1.0;
//};
//
//
//void Simulation::ForcesCalculation(Particle* particleI, Particle* particleJ){
//
//	unsigned int typeParticle_i = particleI->getParticleType();
//	unsigned int typeParticle_j = particleJ->getParticleType();
//
//	double qijx = (particleI->getPositionX() - particleJ->getPositionX());
//	double qijy = (particleI->getPositionY() - particleJ->getPositionY());
//	double qijz = (particleI->getPositionZ() - particleJ->getPositionZ());
//
//	double qijxOld = (particleI->getPositionXOld() - particleJ->getPositionXOld());
//	double qijyOld = (particleI->getPositionYOld() - particleJ->getPositionYOld());
//	double qijzOld = (particleI->getPositionZOld() - particleJ->getPositionZOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//		qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
//		qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
//		qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	double dxOld = qijxOld / rOld;
//	double dyOld = qijyOld / rOld;
//	double dzOld = qijzOld / rOld;
//
//	double act = 0.0;
//	double actOld = 0.0;
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j) /*- InteractionPotential(rOld, typeParticle_i, typeParticle_j)*/;
//
//	//particleI->setF(particleI->getFX() - dxOld*actOld + dx*act, particleI->getFY() - dyOld*actOld + dy*act, particleI->getFZ() - dzOld*actOld + dz*act);
//	//particleJ->setF(particleJ->getFX() - (-dxOld*actOld + dx*act), particleJ->getFY() - (-dyOld*actOld + dy*act), particleJ->getFZ() - (-dzOld*actOld + dz*act));
//
//	particleI->setF(particleI->getFX() + dx*act, particleI->getFY() + dy*act, particleI->getFZ() + dz*act);
//	particleJ->setF(particleJ->getFX() - dx*act, particleJ->getFY() - dy*act, particleJ->getFZ() - dz*act);
//
//};
//
//void Simulation::ForcesCalculationRestrainedParticles(int i, int j){
//
//
//	unsigned int typeParticle_i = pParticles[i]->getParticleType();
//	unsigned int typeParticle_j = pParticles[j]->getParticleType();
//
//	double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
//	double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
//	double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());
//
//	double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
//	double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
//	double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//		qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
//		qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
//		qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	double dxOld = qijxOld / rOld;
//	double dyOld = qijyOld / rOld;
//	double dzOld = qijzOld / rOld;
//
//	double act = 0.0;
//	double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j) - InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	pParticles[i]->setF(pParticles[i]->getFX() - dxOld*actOld + dx*act, pParticles[i]->getFY() - dyOld*actOld + dy*act, pParticles[i]->getFZ() - dzOld*actOld + dz*act);
//	//	pParticles[j]->setF(pParticles[j]->getFX() - (-dxOld*actOld + dx*act), pParticles[j]->getFY() - (-dyOld*actOld + dy*act));
//
//	//return 1.0;
//};
//
//
//
//void Simulation::ForcesCalculationRestrainedParticlesNewton(int i, int j){
//
//
//	unsigned int typeParticle_i = pParticles[i]->getParticleType();
//	unsigned int typeParticle_j = pParticles[j]->getParticleType();
//
//	double qijx = (pParticles[i]->getPositionX() - pParticles[j]->getPositionX());
//	double qijy = (pParticles[i]->getPositionY() - pParticles[j]->getPositionY());
//	double qijz = (pParticles[i]->getPositionZ() - pParticles[j]->getPositionZ());
//
//	double qijxOld = (pParticles[i]->getPositionXOld() - pParticles[j]->getPositionXOld());
//	double qijyOld = (pParticles[i]->getPositionYOld() - pParticles[j]->getPositionYOld());
//	double qijzOld = (pParticles[i]->getPositionZOld() - pParticles[j]->getPositionZOld());
//
//
//	if (PBC == 1){
//		qijx = qijx - BoxLength*roundMy(qijx / BoxLength);
//		qijy = qijy - BoxLength*roundMy(qijy / BoxLength);
//		qijz = qijz - BoxLength*roundMy(qijz / BoxLength);
//
//		qijxOld = qijxOld - BoxLength*roundMy(qijxOld / BoxLength);
//		qijyOld = qijyOld - BoxLength*roundMy(qijyOld / BoxLength);
//		qijzOld = qijzOld - BoxLength*roundMy(qijzOld / BoxLength);
//	}
//
//	double r = sqrt(qijx*qijx + qijy*qijy + qijz*qijz);
//	double rOld = sqrt(qijxOld*qijxOld + qijyOld*qijyOld + qijzOld*qijzOld);
//
//	double dx = qijx / r;
//	double dy = qijy / r;
//	double dz = qijz / r;
//
//	double dxOld = qijxOld / rOld;
//	double dyOld = qijyOld / rOld;
//	double dzOld = qijzOld / rOld;
//
//	double act = 0.0;
//	double actOld = 0.0;
//
//
//	act = InteractionForce(r, typeParticle_i, typeParticle_j);
//	actOld = InteractionForce(rOld, typeParticle_i, typeParticle_j);
//
//	if (typeParticle_i + typeParticle_j == 2) 			 {
//
//		DimerDistance = DimerDistance + r;
//		DimerPotential = DimerPotential + InteractionPotential(r, typeParticle_i, typeParticle_j);
//	}
//
//	potentialEnergy = potentialEnergy + InteractionPotential(r, typeParticle_i, typeParticle_j) - InteractionPotential(rOld, typeParticle_i, typeParticle_j);
//
//	pParticles[i]->setF(pParticles[i]->getFX() - dxOld*actOld + dx*act, pParticles[i]->getFY() - dyOld*actOld + dy*act, pParticles[i]->getFZ() - dzOld*actOld + dz*act);
//	pParticles[j]->setF(pParticles[j]->getFX() - (-dxOld*actOld + dx*act), pParticles[j]->getFY() - (-dyOld*actOld + dy*act), pParticles[j]->getFZ() - (-dzOld*actOld + dz*act));
//
//	//return 1.0;
//};


double Simulation::InteractionForce(double r, unsigned int typeParticle_i, unsigned int typeParticle_j){

	//for (int i = 0; i < 10000; i++);

	unsigned int caseType = (typeParticle_j + typeParticle_i);

	//if (NumberOfDimerParticles == 0 && (typeParticle_i == 1 || typeParticle_j == 1)) std::cout << "Error in particles types..." << std::endl;

	switch (caseType){

	case 0:{ //solvent/solvent interactions


#ifdef SOLVENT_LJ

			   if (r > cut_off) return 0.0;
			   else return fWCA(r);
			   break;
#endif


#ifdef SOLVENT_COS

			   if (r > cut_off) return 0.0;
			   else return -sin(r*PI/cut_off);
			   break;
#endif

	}

	case 1: {//dimer/solvent interactions


				if (r > cut_off) return 0.0;
				else return fWCA(r);
				break;
	}

	case 2://dimer/dimer interactions
	{


			   return fDW(r);//4*(r*r-1);//
			   break;
	}
	}


};

double Simulation::InteractionPotential(double r, unsigned int typeParticle_i, unsigned int typeParticle_j){


	unsigned int caseType = (typeParticle_j + typeParticle_i);

	switch (caseType){

	case 0:{ //solvent/solvent interactions

#ifdef SOLVENT_LJ

			   if (r > cut_off) return 0.0;
			   else return WCA(r);
			   break;

#endif


#ifdef SOLVENT_COS

			   if (r > cut_off) return 0.0;
			   else return cos(r*PI/cut_off)*cut_off/PI;
			   break;

#endif

	}

	case 1: {//dimer/solvent interactions


				if (r > cut_off) return 0.0;
				else return WCA(r);
				break;
	}

	case 2://dimer/dimer interactions
	{


			   return DW(r);//4*(r*r-1);//
			   break;
	}
	}

};




//------------ WCA potential energy ----------
double Simulation::WCA(double x)
{
#ifdef NO_FORCES

	return 0.0;

#endif

	double pot;
	if (x > cut_off)
		pot = 0;
	else pot = epsilon_LJ + 4 * epsilon_LJ*(pow(sigma_LJ / x, 12.) - pow(sigma_LJ / x, 6.));
	return  pot;


}

//------------ WCA force -------------
double Simulation::fWCA(double x)
{
#ifdef NO_FORCES

	return 0.0;

#endif
	double ff;
	if (x > cut_off)
		ff = 0;
	else ff = -24.*epsilon_LJ*(2 * pow(sigma_LJ, 12.)*pow(1. / x, 13.) - pow(sigma_LJ, 6)*pow(1. / x, 7.));
	return  ff;
}

//----------- Dimer potential -------------
double Simulation::DW(double x)
{
#ifdef NO_FORCES

	return 0.0;

#endif
	double pot;
	double dd = (x - cut_off - wDimer*sigma_LJ) / (wDimer*sigma_LJ);
	pot = hDimer*pow(1 - pow(dd, 2.), 2.);//pow(dd, 2.0);//
	return pot;
}

//----------- Dimer force -----------------
double Simulation::fDW(double x)
{
#ifdef NO_FORCES

	return 0.0;

#endif
	double ff;
	double dd = (x - cut_off - wDimer*sigma_LJ) / (wDimer*sigma_LJ);
	ff = -4 * hDimer*(1 - pow(dd, 2.))*dd / (wDimer*sigma_LJ);//2 * dd / (wDimer*sigma_LJ);//
	return ff;
}


void Simulation::WritePositionForBlockAveraging(){

#if BA_STYLE == 1
	if(current_n >= (NumberOfTimeSteps-MAX_ARRAY_LENGTH)) {

		positionBA[BA_count]=ComputePotentialEnergy();
		//positionBA[current_n-(NumberOfTimeSteps-MAX_ARRAY_LENGTH)]=ComputePotentialEnergy();
		BA_count++;
	}

#endif

#if BA_STYLE == 2


	if (((current_n %  WritingPeriodBlockaveraging == 0)) && BA_count < MAX_ARRAY_LENGTH) {

		//std::cout<<BA_count<<std::endl;
		positionBA_scaled[BA_count] = Observable1();
		positionBA_scaled2[BA_count] = Observable2();

		BA_count++;



	}


#endif


};

void Simulation::WritePositionsInFile(){

	std::ofstream FilePositions("CurrentPositionsParticles.m");


		FilePositions << "L=" << BoxLength << ";\n NrDimer=" << NumberOfDimerParticles << ";\n" << "q=[\n";


		for (int i = 0; i < NumberOfParticles; i++)
		{
			double qx = pParticles[i]->getPositionX();
			double qy = pParticles[i]->getPositionY();
			double qz = pParticles[i]->getPositionZ();
			unsigned int isRestrained = pParticles[i]->getActiveStatus();

			FilePositions << qx << "\t" << qy << "\t" << qz << "\t" << isRestrained << "\n";

		}
		FilePositions << "];\n" << "\n";




};


void Simulation::WriteOnScreen(){

	if ((current_n != 0 && (current_n % WritingPeriodOnScreen == 0)) || errorIndicator == true) {

		//grid->print();

		std::cout << "dt= " << dt << "\n epsR/epsF " << pParticles[NumberOfDimerParticles + 1]->getEpsr() << "/" << pParticles[NumberOfDimerParticles + 1]->getEpsf() << " % restrained part = " << RestrainedParticles.getAverage() << "%" << std::endl;
		std::cout << "Temperature: " << Temperature.getAverage() << "\t";//ComputeTemperature()<<"\t";//
	//	std::cout << "Config temperature: " << ConfigTemperature.getAverage() << std::endl;//ComputeTemperature()<<"\t";//

		std::cout << "Potential: " << Potential.getAverage() << "\t";//ComputeTemperature()<<"\t";//
		std::cout << "Position: " << Position.getAverage() << std::endl;//ComputeTemperature()<<"\t";//

		std::cout << "SS A_1(q): " << VarianceOne.getAverage() << "\t";
		std::cout << "D A_2(q): " << VarianceTwo.getAverage() << "\t";
		std::cout << "DS A_3(q): " << VarianceThree.getAverage() << std::endl;
		std::cout << "Nr Interactions: " << NumberOfInteractions.getAverage() << std::endl;
		std::cout << "Nr Neighbors: " << AverageNumberOfNeighbors.getAverage() << std::endl;
		std::cout << "Dimer transition time: " << TransitionTime.getAverage() << std::endl;


		//	std::cout << "Check of computation of potential energy: " << potentialEnergy - ComputePotentialEnergy() << std::endl;

			//potentialEnergy - SolventSolventPotential - DimerPotential - DimerSolventPotential << std::endl;//ComputePotentialEnergy()- SolventSolventPotential - DimerPotential - DimerSolventPotential << std::endl;


		/*
		for (int i = 0; i < NumberOfParticles; i++){

			if (pParticles[i]->getMomentumX() > 10)
			std::cout << "particle "<<i<< "has momentum x "<< pParticles[i]->getMomentumX() << std::endl;

			if (pParticles[i]->getMomentumY() > 10)
				std::cout << "particle "<<i<< "has momentum y "<< pParticles[i]->getMomentumY() << std::endl;

			if (pParticles[i]->getMomentumZ() > 10)
				std::cout << "particle "<<i<< "has momentum z "<< pParticles[i]->getMomentumZ() << std::endl;

			if (pParticles[i]->getPositionX() > BoxLength)
				std::cout << "particle "<<i<< "has position x "<< pParticles[i]->getPositionX() << std::endl;

			if (pParticles[i]->getPositionY() > BoxLength)
				std::cout << "particle " << i << "has position y " << pParticles[i]->getPositionY() << std::endl;


			if (pParticles[i]->getPositionZ() > BoxLength)
				std::cout << "particle " << i << "has position z " << pParticles[i]->getPositionZ() << std::endl;

		}
		*/

#ifdef KINETIC_ENERGY_DOMAIN_COUNTERS
		//counters for domain identification - in which part of the kinetic energy I am the most of the time

		std::cout << "Full: " << FullDynamicsDomain.getAverage() << "\t";
		std::cout << "Restr: " << RestrainedDynamicsDomain.getAverage() << "\t";
		std::cout << "Spline: " << SplineDomain.getAverage() << std::endl;




#endif



#ifdef NVE


		double eCurrent = ComputeEnergy();// ComputeKineticEnergy();


		double errorEnergyPercentage = 100 * abs((eCurrent - intialEnergyNVE)) / abs(intialEnergyNVE);


		std::cout << "Average rel. error: " << EnergyErrorNVE.getAverage() << "%" << std::endl;
		std::cout << "Rel. error: " << errorEnergyPercentage << "%" << std::endl;
#endif

#ifdef METROPOLIZATION
		std::cout << "Rejection rate Determ " << RejectionRate.getAverage() << std::endl;
		std::cout << "Rejection rate FD " << RejectionRateFD.getAverage() << std::endl;

#endif

		std::cout << std::endl;

	}



};

#ifdef VARIANCE_BY_AUTOCORRELATION_FUNCTION

void Simulation::WriteAutocorrelationAveragedInReplicas(){

	int maxLengthOfData = MAX_LENGTH_OF_VARIANCE_ARRAY; // for writing convergence of variance in number of replicas

	int modTimeInterval1 = (current_n % (int)NumberOfStepsInTimeInterval1);
	int modTimeInterval2 = (current_n % (int)NumberOfStepsInTimeInterval2);
	int modTimeInterval3 = (current_n % (int)NumberOfStepsInTimeInterval3);

	int arrayVarianceCounter1 = ((int)countReplicasAutocorrelation1 % maxLengthOfData); //move in the array of variance onver number of replicas
	int arrayVarianceCounter2 = ((int)countReplicasAutocorrelation2 % maxLengthOfData);
	int arrayVarianceCounter3 = ((int)countReplicasAutocorrelation3 % maxLengthOfData);

	if (arrayVarianceCounter1 > maxLengthOfData || arrayVarianceCounter2 > maxLengthOfData || arrayVarianceCounter3 > maxLengthOfData) std::cout << "error in WriteAutocorrelationAveragedInReplicas" << std::endl;

	double writingToFile = (current_n % WRITING_AUTOCORRELATION_TO_FILE);


#ifdef BACKUP_VARIANCE_EVOLUTION

	if (arrayVarianceCounter1 == 0 && current_n != 0) {

		//std::cout << countReplicasAutocorrelation << std::endl;
		//while (1);

		std::ofstream fileAR1("VARIANCE_IN_NUMBER_OF_REPLICAS_1_tmp.m");
		fileAR1 << "NumberOfReplicas=" << countReplicasAutocorrelation1 << ";\n";
		fileAR1 << "Variance1_tmp=[";

		for (int i = 0; i < maxLengthOfData; i++){

			double tmp = ConvergenceEstimatedVarianceInNumberOfReplicas1[i];
			fileAR1 << tmp<< " ";

		}
		fileAR1 << "];\n" << std::endl;

		//for (int i = 0; i < maxLengthOfData; i++)ConvergenceEstimatedVarianceInNumberOfReplicas[i] = 0.0;
	}
	if (arrayVarianceCounter2 == 0 && current_n != 0) {

		std::ofstream fileAR22("VARIANCE_IN_NUMBER_OF_REPLICAS_2_tmp.m");

		//std::cout << countReplicasAutocorrelation_2 << std::endl;

		fileAR22 << "NumberOfReplicas2=" << countReplicasAutocorrelation2 << ";\n";
		fileAR22 << "Variance2_tmp=[";

		for (int i = 0; i < maxLengthOfData; i++){

			double tmp = ConvergenceEstimatedVarianceInNumberOfReplicas2[i];
			fileAR22 <<  tmp<< " ";

		}

		fileAR22 << "];\n" << std::endl;


		//for (int i = 0; i < maxLengthOfData; i++)	ConvergenceEstimatedVarianceInNumberOfReplicas2[i] = 0.0;

	}
	if (arrayVarianceCounter3 == 0 && current_n != 0) {

		//std::cout << countReplicasAutocorrelation << std::endl;
		//while (3);

		std::ofstream fileAR3("VARIANCE_IN_NUMBER_OF_REPLICAS_3_tmp.m");
		fileAR3 << "NumberOfReplicas=" << countReplicasAutocorrelation3 << ";\n";
		fileAR3 << "Variance3_tmp=[";

		for (int i = 0; i < maxLengthOfData; i++){

			double tmp = ConvergenceEstimatedVarianceInNumberOfReplicas3[i];
			fileAR3 << tmp << " ";

		}
		fileAR3 << "];\n" << std::endl;

		//for (int i = 0; i < maxLengthOfData; i++)ConvergenceEstimatedVarianceInNumberOfReplicas[i] = 0.0;
	}

#endif

	// define new initial condition

	if (modTimeInterval1 == 0)		InitialValueAutocorrelation1 = Observable1();
	if (modTimeInterval2 == 0)		InitialValueAutocorrelation2 = Observable2();
	if (modTimeInterval3 == 0)		InitialValueAutocorrelation3 = Observable3();


	// add sample to the average for given time step A0*An

	AutocorrelationArrayReplica1[(int)modTimeInterval1].addSample(InitialValueAutocorrelation1 *Observable1());
	AutocorrelationArrayReplica2[(int)modTimeInterval2].addSample(InitialValueAutocorrelation2 *Observable2());
	AutocorrelationArrayReplica3[(int)modTimeInterval3].addSample(InitialValueAutocorrelation3 *Observable3());


	AveragesArrayReplica1[(int)modTimeInterval1].addSample(Observable1());
	AveragesArrayReplica2[(int)modTimeInterval2].addSample(Observable2());
	AveragesArrayReplica3[(int)modTimeInterval3].addSample(Observable3());



#ifdef ACF_TRAPEZ
	// ***compute the variance from autocorrelation function for first observable

	if (modTimeInterval1 == 0) {


		for (int i = 0; i < NumberOfStepsInTimeInterval1; i++) {

			if (i == 0 || i == NumberOfStepsInTimeInterval1 - 1)
				IntegralAutocorrelation1 = IntegralAutocorrelation1 + 0.5*(AutocorrelationArrayReplica1[i].getAverage() - VarianceOne.getAverage()*VarianceOne.getAverage());
			else
				IntegralAutocorrelation1 = IntegralAutocorrelation1 + AutocorrelationArrayReplica1[i].getAverage() - VarianceOne.getAverage()*VarianceOne.getAverage();
		}

		IntegralAutocorrelation1 = 2.0*IntegralAutocorrelation1*dt;

		//save to variance vector (to see the convergence, i.e. the most precise value is the last one)

		ConvergenceEstimatedVarianceInNumberOfReplicas1[(int)arrayVarianceCounter1] = IntegralAutocorrelation1;

		//add replicas count

		countReplicasAutocorrelation1++;


	}



	//***compute the variance from autocorrelation function for second observable

	if (modTimeInterval2 == 0) {


		for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {

			if (i == 0 || i == NumberOfStepsInTimeInterval2 - 1)
				IntegralAutocorrelation2 = IntegralAutocorrelation2 + 0.5*(AutocorrelationArrayReplica2[i].getAverage() - VarianceTwo.getAverage()*VarianceTwo.getAverage());
			else
				IntegralAutocorrelation2 = IntegralAutocorrelation2 + AutocorrelationArrayReplica2[i].getAverage() - VarianceTwo.getAverage()*VarianceTwo.getAverage();
		}
		IntegralAutocorrelation2 = 2.0*IntegralAutocorrelation2*dt;

		ConvergenceEstimatedVarianceInNumberOfReplicas2[(int)arrayVarianceCounter2] = IntegralAutocorrelation2;

		countReplicasAutocorrelation2++;


	}
	if (modTimeInterval3 == 0) {


		IntegralAutocorrelation3 = 0;
		for (int i = 0; i < NumberOfStepsInTimeInterval3; i++) {

			if (i == 0 || i == NumberOfStepsInTimeInterval3 - 1)
				IntegralAutocorrelation3 = IntegralAutocorrelation3 + 0.5*(AutocorrelationArrayReplica3[i].getAverage() - VarianceThree.getAverage()*VarianceThree.getAverage());
			else
				IntegralAutocorrelation3 = IntegralAutocorrelation3 + AutocorrelationArrayReplica3[i].getAverage() - VarianceThree.getAverage()*VarianceThree.getAverage();
		}
		IntegralAutocorrelation3 = 2.0*IntegralAutocorrelation3*dt;

		//save to variance vector (to see the convergence, i.e. the most precise value is the last one)

		ConvergenceEstimatedVarianceInNumberOfReplicas3[(int)arrayVarianceCounter3] = IntegralAutocorrelation3;

		//add replicas count

		countReplicasAutocorrelation3++;


	}

#endif




//******************************************************************
	/*AutocorrelationArrayReplica1[(int)modTimeInterval1].addSample((InitialValueAutocorrelation1 - VarianceOne.getAverage())*(Observable1() - VarianceOne.getAverage()));
	AutocorrelationArrayReplica2[(int)modTimeInterval2].addSample((InitialValueAutocorrelation2 - VarianceTwo.getAverage())*(Observable2() - VarianceTwo.getAverage()));
	AutocorrelationArrayReplica3[(int)modTimeInterval3].addSample((InitialValueAutocorrelation3 - VarianceThree.getAverage())*(Observable3() - VarianceThree.getAverage()));
	*/

//
//#ifdef ACF_TRAPEZ
//	//***compute the variance from autocorrelation function for first observable
//
//	if (modTimeInterval1 == 0) {
//
//
//		for (int i = 0; i < NumberOfStepsInTimeInterval1; i++) {
//
//			if (i == 0 || i == NumberOfStepsInTimeInterval1-1)
//				IntegralAutocorrelation1 = IntegralAutocorrelation1 + 0.5*AutocorrelationArrayReplica1[i].getAverage();
//			else
//				IntegralAutocorrelation1 = IntegralAutocorrelation1 + AutocorrelationArrayReplica1[i].getAverage();
//		}
//		IntegralAutocorrelation1 = 2.0*IntegralAutocorrelation1*dt;
//
//		//save to variance vector (to see the convergence, i.e. the most precise value is the last one)
//
//		ConvergenceEstimatedVarianceInNumberOfReplicas1[(int)arrayVarianceCounter1] = IntegralAutocorrelation1;
//
//		//add replicas count
//
//		countReplicasAutocorrelation1++;
//
//
//	}
//
//	//***compute the variance from autocorrelation function for second observable
//
//	if (modTimeInterval2 == 0) {
//
//
//		for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {
//
//			if (i == 0 || i == NumberOfStepsInTimeInterval2 - 1)
//				IntegralAutocorrelation2 = IntegralAutocorrelation2 + 0.5*AutocorrelationArrayReplica2[i].getAverage();
//			else
//				IntegralAutocorrelation2 = IntegralAutocorrelation2 + AutocorrelationArrayReplica2[i].getAverage();
//		}
//		IntegralAutocorrelation2 = 2.0*IntegralAutocorrelation2*dt;
//
//		ConvergenceEstimatedVarianceInNumberOfReplicas2[(int)arrayVarianceCounter2] = IntegralAutocorrelation2;
//
//		countReplicasAutocorrelation2++;
//
//
//	}
//	if (modTimeInterval3 == 0) {
//
//
//		IntegralAutocorrelation3 = 0;
//		for (int i = 0; i < NumberOfStepsInTimeInterval3; i++) {
//
//			if (i == 0 || i == NumberOfStepsInTimeInterval3 - 1)
//				IntegralAutocorrelation3 = IntegralAutocorrelation3 + 0.5*AutocorrelationArrayReplica3[i].getAverage();
//			else
//				IntegralAutocorrelation3 = IntegralAutocorrelation3 + AutocorrelationArrayReplica3[i].getAverage();
//		}
//		IntegralAutocorrelation3 = 2.0*IntegralAutocorrelation3*dt;
//
//		//save to variance vector (to see the convergence, i.e. the most precise value is the last one)
//
//		ConvergenceEstimatedVarianceInNumberOfReplicas3[(int)arrayVarianceCounter3] = IntegralAutocorrelation3;
//
//		//add replicas count
//
//		countReplicasAutocorrelation3++;
//
//
//	}
//
//#endif


	//*** simpson rule:

	//***compute the variance from autocorrelation function for first observable

#ifdef ACF_SIMPSON
	//simpson's rule

	if (modTimeInterval1 == 0) {

		IntegralAutocorrelation1 = AutocorrelationArrayReplica1[0].getAverage() + AutocorrelationArrayReplica1[(int)NumberOfStepsInTimeInterval1 - 1].getAverage();

		for (int i = 2; i < (int)(NumberOfStepsInTimeInterval1)-1 - 2; i = i + 2) {

			IntegralAutocorrelation1 = IntegralAutocorrelation1 + 2.0*AutocorrelationArrayReplica1[i].getAverage();

		}

		for (int i = 1; i < (int)(NumberOfStepsInTimeInterval1)-1; i = i + 2) {

			IntegralAutocorrelation1 = IntegralAutocorrelation1 + 4.0*AutocorrelationArrayReplica1[i].getAverage();

		}

		IntegralAutocorrelation1 = IntegralAutocorrelation1 / 3.0;

		IntegralAutocorrelation1 = 2.0*IntegralAutocorrelation1*dt;

		ConvergenceEstimatedVarianceInNumberOfReplicas1[(int)arrayVarianceCounter1] = IntegralAutocorrelation1;

		countReplicasAutocorrelation1++;



	}



	//simpson's rule

	if (modTimeInterval2 == 0) {

		IntegralAutocorrelation2 = AutocorrelationArrayReplica2[0].getAverage() + AutocorrelationArrayReplica2[(int)NumberOfStepsInTimeInterval2 - 1].getAverage();

		for (int i = 2; i < (int)(NumberOfStepsInTimeInterval2)-1 - 2; i = i + 2) {

			IntegralAutocorrelation2 = IntegralAutocorrelation2 + 2.0*AutocorrelationArrayReplica2[i].getAverage();

		}

		for (int i = 1; i < (int)(NumberOfStepsInTimeInterval2)-1; i = i + 2) {

			IntegralAutocorrelation2 = IntegralAutocorrelation2 + 4.0*AutocorrelationArrayReplica2[i].getAverage();

		}

		IntegralAutocorrelation2 = IntegralAutocorrelation2 / 3.0;

		IntegralAutocorrelation2 = 2.0*IntegralAutocorrelation2*dt;

		ConvergenceEstimatedVarianceInNumberOfReplicas2[(int)arrayVarianceCounter2] = IntegralAutocorrelation2;

		countReplicasAutocorrelation2++;



	}

	//simpson's rule

	if (modTimeInterval3 == 0) {

		IntegralAutocorrelation3 = AutocorrelationArrayReplica3[0].getAverage() + AutocorrelationArrayReplica3[(int)NumberOfStepsInTimeInterval3 - 1].getAverage();

		for (int i = 2; i < (int)(NumberOfStepsInTimeInterval3)-1 - 2; i = i + 2) {

			IntegralAutocorrelation3 = IntegralAutocorrelation3 + 2.0*AutocorrelationArrayReplica3[i].getAverage();

		}

		for (int i = 1; i < (int)(NumberOfStepsInTimeInterval3)-1; i = i + 2) {

			IntegralAutocorrelation3 = IntegralAutocorrelation3 + 4.0*AutocorrelationArrayReplica3[i].getAverage();

		}

		IntegralAutocorrelation3 = IntegralAutocorrelation3 / 3.0;

		IntegralAutocorrelation3 = 2.0*IntegralAutocorrelation3*dt;

		ConvergenceEstimatedVarianceInNumberOfReplicas3[(int)arrayVarianceCounter3] = IntegralAutocorrelation3;

		countReplicasAutocorrelation3++;



	}

#endif



	if (writingToFile == 0){

		// *******  write autocorrelation function into file

		std::ofstream file("AUTOCORRELATION_FUNCTION_1.m");

		file << "dt= " << dt << ";" << std::endl;
		file << "AAverEnd1= " << VarianceOne.getAverage() << ";" << std::endl;
		file << "A1=[";

		for (int i = 0; i < NumberOfStepsInTimeInterval1; i++) {

			//the writing interval is reduced to avoid big data

			//if (i % 100 == 0)
			file << AutocorrelationArrayReplica1[i].getAverage() << " ";
		}
		file << "\n];" << std::endl;

		file << "AAver1=[";
		for (int i = 0; i < NumberOfStepsInTimeInterval1; i++) {

			file << AveragesArrayReplica1[i].getAverage() << " ";
		}
		file << "\n];" << std::endl;

		file << "Variance1 = " << IntegralAutocorrelation1 << ";" << std::endl;


		std::ofstream file2("AUTOCORRELATION_FUNCTION_2.m");

		file2 << "dt= " << dt << ";" << std::endl;
		file2 << "AAverEnd2= " << VarianceTwo.getAverage() << ";" << std::endl;
		file2 << "A2=[";

		for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {
			//if (i % 1 == 0)
			file2 << AutocorrelationArrayReplica2[i].getAverage() << " ";
		}
		file2 << "\n];" << std::endl;


		file2 << "AAver2=[";
		for (int i = 0; i < NumberOfStepsInTimeInterval2; i++) {

			file2 << AveragesArrayReplica2[i].getAverage() << " ";
		}
		file2 << "\n];" << std::endl;

		file2 << "Variance2 = " << IntegralAutocorrelation2 << ";" << std::endl;

		std::ofstream file3("AUTOCORRELATION_FUNCTION_3.m");

		file3 << "dt= " << dt << ";" << std::endl;
		file3 << "AAverEnd3= " << VarianceThree.getAverage() << ";" << std::endl;
		file3 << "A3=[";

		for (int i = 0; i < NumberOfStepsInTimeInterval3; i++) {
			//if (i % 1 == 0)
			file3 << AutocorrelationArrayReplica3[i].getAverage() << " ";
		}
		file3 << "\n];" << std::endl;


		file3 << "AAver3=[";
		for (int i = 0; i < NumberOfStepsInTimeInterval3; i++) {

			file3 << AveragesArrayReplica3[i].getAverage() << " ";
		}
		file3 << "\n];" << std::endl;

		file3 << "Variance3 = " << IntegralAutocorrelation3 << ";" << std::endl;

		// *******  write estimated variance into file

		std::ofstream fileAR("VARIANCE_IN_NUMBER_OF_REPLICAS_1.m");
		fileAR << "NumberOfReplicas=" << countReplicasAutocorrelation1 << ";\n";

		fileAR << "Variance1=" << IntegralAutocorrelation1 << ";" << std::endl;
		//fileAR << "VarianceVect1=[";

		//for (int i = 0; i < arrayVarianceCounter1; i++)	fileAR << ConvergenceEstimatedVarianceInNumberOfReplicas1[i] << " ";
	//	fileAR << "\n];" << std::endl;

		std::ofstream fileAR2("VARIANCE_IN_NUMBER_OF_REPLICAS_2.m");

		fileAR2 << "NumberOfReplicas2=" << countReplicasAutocorrelation2 << ";\n";
		fileAR2 << "Variance2=" << IntegralAutocorrelation2 << ";" << std::endl;
		//fileAR2 << "VarianceVect2=[";

		//for (int i = 0; i < arrayVarianceCounter2; i++)	fileAR2 << ConvergenceEstimatedVarianceInNumberOfReplicas2[i] << " ";

		//fileAR2 << "\n];" << std::endl;

		std::ofstream fileAR3("VARIANCE_IN_NUMBER_OF_REPLICAS_3.m");

		fileAR3 << "NumberOfReplicas3=" << countReplicasAutocorrelation3 << ";\n";
		fileAR3 << "Variance3=" << IntegralAutocorrelation3 << ";" << std::endl;
		//fileAR3 << "VarianceVect3=[";

		//for (int i = 0; i < arrayVarianceCounter3; i++)	fileAR3 << ConvergenceEstimatedVarianceInNumberOfReplicas3[i] << " ";

		//fileAR3 << "\n];" << std::endl;


	}



	if ((current_n != 0 && (current_n % WritingPeriodOnScreen == 0)) || errorIndicator == true) {

		std::cout << "Variance: ";
		std::cout << "SS : " << IntegralAutocorrelation1 << "\t";
		std::cout << "D : " << IntegralAutocorrelation2 << "\t";
		std::cout << "DS : " << IntegralAutocorrelation3 << std::endl;
	}


};

#endif


void Simulation::SampleVarianceIR(){


	ReplicasAverageThroughPathN.addSample(Observable2());

	//if (replicasCounter < (double) NUMBER_OF_REPLICAS){
	if (current_n % (int)LengthOfPartialTrajectoryN == 0){

		ReplicasAverage.addSample(ReplicasAverageThroughPathN.getAverage());

		//if (ReplicasAverage.getCurrentIndex() > 1000)
		VarianceIR.addSample((ReplicasAverageThroughPathN.getAverage() - ReplicasAverage.getAverage())*(ReplicasAverageThroughPathN.getAverage() - ReplicasAverage.getAverage()));

		ReplicasAverageThroughPathN.clean();



	}

	VarianceIndependentReplicas = ((double)LengthOfPartialTrajectoryN) * VarianceIR.getAverage()*dt;

	if (WritingPeriodFile == 0){

		std::ofstream filevARiR("VARIANCE_IR.m");

		filevARiR << VarianceIndependentReplicas << std::endl;
	}


	if (current_n % WritingPeriodOnScreen == 0){
		//Observable2
		std::cout << "\t Variance IR D = " << VarianceIndependentReplicas << std::endl;
	}
	//}





};


void Simulation::SampleVarianceIR2(){


	ReplicasAverageThroughPathN2.addSample(Observable2());

	//if (replicasCounter < (double) NUMBER_OF_REPLICAS){
	if (current_n % (int)LengthOfPartialTrajectoryN2 == 0){

		ReplicasAverage2.addSample(ReplicasAverageThroughPathN2.getAverage());

		//if (ReplicasAverage.getCurrentIndex() > 1000)
		VarianceIR2.addSample((ReplicasAverageThroughPathN2.getAverage() - ReplicasAverage.getAverage())*(ReplicasAverageThroughPathN2.getAverage() - ReplicasAverage.getAverage()));

		ReplicasAverageThroughPathN2.clean();



	}

	VarianceIndependentReplicas2 = ((double)LengthOfPartialTrajectoryN2) * VarianceIR2.getAverage()*dt;

	if (WritingPeriodFile == 0){

		std::ofstream filevARiR("VARIANCE_IR2.m");

		filevARiR << VarianceIndependentReplicas2 << std::endl;
	}


	if (current_n % WritingPeriodOnScreen == 0){
		//Observable2
		std::cout << "\t Variance IR 2 D = " << VarianceIndependentReplicas2 << std::endl;
	}
	//}





};
