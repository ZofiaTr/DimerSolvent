// Zofia Trstanova, September 2014
// last modification February 2015
// Linux modification of original cpp code written for windows - main modification in file writing functions (linux didnt like file objects.. solution - kick them out :-D)

// Langevin with ARPS

#include "Simulation.h"
//#include "SimulationAdaptive.h"

#define ADAPTIVE 0 // 1 adaptive simulation, 0 full forces update
#define ERROR_ANALYSIS_INDEPENDENT_REPLICAS

//#define INPUT_TIME_STEP_SIZE
//#define INPUT_ARPS_PARAMETERS

//#define INPUT_ARPS_PARAMETERS_AND_DT_FROM_MATRIX


#define METROPOLIS 0

int main() {




	// write on the screen

#ifdef  METROPOLIZATION
	std::cout << "METROPOLIZATION" << std::endl;
#endif

#ifdef  FIRST_ORDER_SPLITTING
	std::cout << "FIRST_ORDER_SPLITTING" << std::endl;
#endif

#ifdef  SECOND_ORDER_SPLITTING
	std::cout << "SECOND_ORDER_SPLITTING" << std::endl;
#endif

#ifdef  SECOND_ORDER_FD_CN
	std::cout << "SECOND_ORDER_FD_CN trapez" << std::endl;
#endif

#ifdef  SECOND_ORDER_FD_MIDPOINT
	std::cout << "SECOND_ORDER_FD_MIDPOINT" << std::endl;
#endif
#ifdef ONLY_FLUCTUATION_DISSIPATION_PART
	std::cout << "ONLY_FLUCTUATION_DISSIPATION_PART" << std::endl;
#endif

	bool errorParameters = false;

#ifdef ERROR_ANALYSIS_INDEPENDENT_REPLICAS

#ifdef ONE_PARTICLE_ONLY
	double NrParticles=1;//100;
	double NrDimer = 0;


#endif
#ifndef ONE_PARTICLE_ONLY
	double NrParticles =  27;//100;
	double NrDimer=0;
#endif

	double ratio = 0.0;

	double epsfParticle = 6.0;// 100.0;
	double epsrParticle = 0.0;//ratio*epsfParticle;
	 // eps_f= ratio * eps_r

#ifndef INPUT_TIME_STEP_SIZE

	double timeStepSize = 0.01;//055;

#endif


#ifdef INPUT_TIME_STEP_SIZE

	// read dt value from input file
	std::ifstream INPUT_dt("INPUT_dt");


	double timeStepSize = 0.0;
	INPUT_dt >> timeStepSize >> std::ws;

#endif



#ifdef INPUT_ARPS_PARAMETERS


	std::ifstream INPUT_eps("INPUT_eps");
	//std::ifstream INPUT_epsF("INPUT_epsF");

	std::ifstream INPUT_delta("INPUT_delta");
//
//#ifndef INPUT_TIME_STEP_SIZE
//	double timeStepSize = 0.005;
//#endif


	INPUT_eps >> epsfParticle >> std::ws;
	//INPUT_eps >> epsrParticle >> std::ws;
	//INPUT_epsF >> epsfParticle >> std::ws;
	INPUT_delta >> epsrParticle >> std::ws;

	//epsrParticle = ratio;//*epsfParticle;

	// if (epsrParticle > epsfParticle) {
	//
	// 	std::cout << "Error in Parameters: eps_r< eps_f!" << std::endl;
	// 	std::cout << "Epsr= " << epsrParticle << "\tEpsf=" << epsfParticle << std::endl;
	// 	errorParameters = true;
	// }

#ifdef ONE_PARTICLE_ONLY

	std::cout <<  std::endl;
	std::cout << "-------------------------------------------------------" << std::endl;
	std::cout << "\n2D SIMULATION WITH ONE PARTICLE IN COSINUS POTENTIAL\n" << std::endl;
	std::cout << "-------------------------------------------------------" << std::endl;

#endif


#endif


#ifdef INPUT_ARPS_PARAMETERS_AND_DT_FROM_MATRIX

	// read from INPUT_NR_EPS and INPUT_NR_DT number various time steps corresponding to eps coeff-> each row is for one epsr and contains corresponding dt
	//in file MATRIX are arps parameters with corresponding dt s
	//MATRIX(NrEps, NrDt)
	// (INPUT_NR_EPS , INPUT_NR_DT) size of MATRIX -> read it all to array_2d
	// (INDEX_EPS, INDEX_DT) choose dt at this position from MATRIX


	int numberEpsR = 0;
	int numberDt = 0;
	int indexDt = 0;
	int indexEps = 0;

	std::cout << "reading matrix from file.." << std::endl;
	std::ifstream MATRIX("MATRIX");

	if (!MATRIX){
		std::cout << "Error: cannot find file MATRIX!" << std::endl;
		while (1);
	}

	std::ifstream INPUT_NR_EPS("INPUT_NR_EPS");
	INPUT_NR_EPS >> numberEpsR >> std::ws;

	std::ifstream INPUT_NR_DT("INPUT_NR_DT");
	INPUT_NR_DT >> numberDt >> std::ws;

	std::ifstream INDEX_EPS("INDEX_EPS");
	INDEX_EPS >> indexEps >> std::ws;

	std::ifstream INDEX_DT("INDEX_DT");
	INDEX_DT >> indexDt >> std::ws;

	std::cout << "Number of eps "<<numberEpsR << "\t Number of dt " << numberDt << std::endl;


	double ** array_2d = new double*[numberEpsR];
	for (int i = 0; i < numberEpsR; i++) {
		array_2d[i] = new double[numberDt];

		for (int j = 0; j < numberDt; j++) {
			MATRIX >> array_2d[i][j];
			std::cout << array_2d[i][j] << "\t";
		}
		std::cout <<std::endl;
	}


	timeStepSize =array_2d[indexEps - 1][indexDt - 1];// array_2d[numberEpsR - 1][numberDt - 1];//

	std::cout << "dt ="  << timeStepSize << std::endl;

	std::cout << "from MATRIX: dt=" << timeStepSize << std::endl;
	std::cout << " epsR=" << epsrParticle << std::endl;

	//while (1);

#endif


	double initialDensity = 1.0;

	// physical time interval (0, time]
	double time = pow(10.0, 9);
	//only initialization, then numberOfTimeSteps=(time/dt);
	double NumberOfTimeSteps=1.0;

	int numberOfReplicas=1;//100000;
	double numberOfEquilibrationSteps =  10000;

	//minimal and maximal time step size
	double dt_min = timeStepSize;
	double dt_max = timeStepSize;
	double deltat = timeStepSize;

	//double NrParticles=numberOfReplicas;

#if ADAPTIVE == 1
	SimulationAdaptive* NVT;
	NVT=new SimulationAdaptive(NrParticles);
#endif

#if ADAPTIVE == 0
	Simulation* NVT;
	NVT=new Simulation(NrParticles);
#endif



	//how many different dt..
	int NrOfDifferent_dt=0;
	for(double iter=dt_min; iter<=dt_max; iter=iter+deltat)  NrOfDifferent_dt=NrOfDifferent_dt+1;
	double dt_iter_count=0;

	std::cout<<"Number of different dt "<<NrOfDifferent_dt<<std::endl;


	// FILE init.
	std:: ofstream fileSamples ("data/FINAL");

#ifdef WRITE_AVER_TIME_PER_TIME_STEP
	std::ofstream fileTimeMeasurings("data/TIME_MEASURMENT_ADAPTIVELY_NEWTON");
#endif

//#ifdef NO_FORCES
//	std::ofstream fileTimeMeasurings("data/TIME_MEASURMENT_NO_FORCES");
//#endif

	//MATLAB format with % comment.. to desribe matrix entries -> columns- observables, rows- different dt

	/*
	//Writing to data/FINAL

	fileSamples << "%samples[0] = (double)(current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN); " << std::endl;
	fileSamples << "%samples[1] = Temperature.getAverage();" << std::endl;
	fileSamples << "%samples[2] = Pressure.getAverage();" << std::endl;
	fileSamples << "%samples[3] = Position.getAverage();" << std::endl;
	fileSamples << "%samples[4] = Potential.getAverage();" << std::endl;
	fileSamples << "%samples[5] = RestrainedParticles.getAverage();" << std::endl;
	fileSamples << "%samples[6] = VarianceOne.getAverage();" << std::endl;
	fileSamples << "%samples[7] = VarianceTwo.getAverage();" << std::endl;
	fileSamples << "%samples[8] = VarianceThree.getAverage();" << std::endl;
	fileSamples << "%samples[9] = ForceAppliedOnDimer.getAverage();" << std::endl;
	fileSamples << "%samples[10] = NumberOfInteractions.getAverage();" << std::endl;
	fileSamples << std::endl;
	*/

	fileTimeMeasurings << "%Time per time step depending on the percentage of restrained particles \n";

	double *varianceArray=new double[NrOfDifferent_dt];
	double *varianceArray2 = new double[NrOfDifferent_dt];
	double *dtArray=new double[NrOfDifferent_dt];


	//	std::cout<<"epsParticle="<<epsrParticle<<std::endl;

	// time variables
	double t0=0.;
	double time_for_trajectory=0.;

	int numberOfSamples = NUMBER_OF_AVERAGES;
	double*samplesVector = new double[numberOfSamples];


	/*samples[0] = (double)(current_n + (double)numberOfCurrent_n*(double)MaxLengthOfCurrentN);
	samples[1] = Temperature.getAverage();
	samples[2] = Pressure.getAverage();
	samples[3] = Position.getAverage();
	samples[4] = Potential.getAverage();
	samples[5] = RestrainedParticles.getAverage();
	samples[6] = VarianceOne.getAverage();
	samples[7] = VarianceTwo.getAverage();
	samples[8] = VarianceThree.getAverage();
	samples[9] = ForceAppliedOnDimer.getAverage();
	samples[10] = NumberOfInteractions.getAverage();
	*/

		NVT->setNumberStepsFORloop(time, dt_min, NrOfDifferent_dt, deltat);

	//-------------------------------------------------
	// TIME STEP SIZE LOOP

	int countdt=0;
	double *averageValueReplica=new double [numberOfReplicas];
	double *averageValueReplicaTwo = new double[numberOfReplicas];
	double *averageValueReplicaThree = new double[numberOfReplicas];


	std::cout << "Physical time: " << time << " Time step sizes from dtMin = " << dt_min << " to dtMax = " << dt_max << " by delta = " << deltat << std::endl;
	NumberOfTimeSteps = (time / dt_min);


#ifdef  OPTIMAL_PARAMETERS

	double** threePoints = new double*[3];
	for (int i = 0; i < 3; i++) threePoints[i] = new double[2];

	//three points for (epsR, epsF): (r1,f1),(r2,f2),(r3,f3) with r1=r2 and f1=f3!
	double cStab = 0.8;

	threePoints[0][1] = 10;//3.2; //f1
	threePoints[0][0] = cStab*threePoints[0][1];//0.8*3.2; //r1

	threePoints[1][1] = 5;//1.0; //r2
	threePoints[1][0] = cStab*threePoints[1][1];//3;

	threePoints[2][1] = 3;//2.0;
	threePoints[2][0] = cStab*threePoints[2][1];//3.2; //f3

	//if (abs(threePoints[0][0] - threePoints[1][0]) > 0) std::cout << "Error: r1 noneq r2\n" << std::endl;
	//if (abs(threePoints[0][1] - threePoints[2][1]) > 0) std::cout << "Error: f1 noneq f3\n" << std::endl;

	if (errorParameters == false){

		std::cout << "SEARCH OPTIMAL PARAMETERS " << std::endl;

		std::ofstream fileOptimalParameters("OPTIMAL_PARAMETERS.m");

		//run reference simulation with standard
#ifdef RUN_REFERENCE_SIMULATION_WITH_STD
		NVT->setSimulationParameters(NumberOfTimeSteps, dt_min, 0, 0, NrDimer, numberOfReplicas, initialDensity, numberOfEquilibrationSteps);
		NVT->InitialCondition();
		NVT->PerformTrajectory();

		double referenceStandardDeviation = NVT->getStandardDeviation();

		fileOptimalParameters << "VarianceSTD = " << referenceStandardDeviation*referenceStandardDeviation << ";\n" << std::endl;


		NVT->CleanSamples();

		NVT->~Simulation();

		double varStd = referenceStandardDeviation*referenceStandardDeviation;
#endif
		//	bool weCouldDoBetter = true;

		//while (weCouldDoBetter == true){

		fileOptimalParameters << "M=[" << std::endl;

		double *convergenceSpeedUpVector = new double [3];
		double *varianceVector = new double[3];
		double *percRestrVector = new double[3];


		for (int i = 0; i < 3; i++){

			Simulation* NVT;
			NVT = new Simulation(NrParticles);

			epsrParticle = threePoints[i][0];
			epsfParticle = threePoints[i][1];


			NVT->setSimulationParameters(NumberOfTimeSteps, dt_min, epsrParticle, epsfParticle, NrDimer, numberOfReplicas, initialDensity, numberOfEquilibrationSteps);
			NVT->InitialCondition();
			NVT->PerformTrajectory();

		double standardDeviationARPS = NVT->getStandardDeviation();

		// linear interpolation of the variance

		double varArps = standardDeviationARPS*standardDeviationARPS;
		varianceVector[i] = varArps;
		//double y = varStd / (x*(varArps - varStd) / epsR + varStd);

		double percentageOfRestrained = NVT->GetPercentageOfRestrained();


		//if (convergenceSpeedUp > 1){
#ifdef RUN_REFERENCE_SIMULATION_WITH_STD
		double convergenceSpeedUp = 0.0;
		if (percentageOfRestrained > 0.7)
			convergenceSpeedUp = sqrt(0.25 / (1 - 0.01*percentageOfRestrained))*referenceStandardDeviation / standardDeviationARPS;
		else
			convergenceSpeedUp = 1 * referenceStandardDeviation / standardDeviationARPS;

		convergenceSpeedUpVector[i] = convergenceSpeedUp;
		percRestrVector[i] = percentageOfRestrained;

		std::cout << "Convergence speed up is " << convergenceSpeedUp << std::endl;
		std::cout << "referenceStandardDeviation is " << referenceStandardDeviation << std::endl;
#endif

		std::cout << "percentageOfRestrained is " << percentageOfRestrained << std::endl;
		std::cout << "standardDeviationARPS is " << standardDeviationARPS << std::endl;
		//	}



		fileOptimalParameters << epsrParticle << "\t" << epsfParticle << "\t" << percentageOfRestrained << "\t" << standardDeviationARPS*standardDeviationARPS << "\n" << std::endl;

		NVT->CleanSamples();
		NVT->~Simulation();


	}

		std::cout << "Variance and convergence speed up three simulations: \n" << std::endl;
		std::cout << "1) Var=" << varianceVector[0]<<"\tTotal speed up= " << convergenceSpeedUpVector[0] <<"\t % restr"<<percRestrVector[0]<< std::endl;
		std::cout << "2) Var=" << varianceVector[1] << "\tTotal speed up=" << convergenceSpeedUpVector[1] << "\t % restr" << percRestrVector[1] << std::endl;
		std::cout << "3) Var=" << varianceVector[2] << "\tTotal speed up=" << convergenceSpeedUpVector[2] << "\t % restr" << percRestrVector[2] << std::endl;

	fileOptimalParameters << "];\n" << std::endl;

				//while (1);
	//	}
	}


#endif


#ifndef OPTIMAL_PARAMETERS

	if (errorParameters == false){

		for (double dt_iter = dt_min; dt_iter <= dt_max; dt_iter = dt_iter + deltat){

			std::cout << "Physical time: " << time << " Time step sizes from dtMin = " << dt_min << " to dtMax = " << dt_max << " by delta = " << deltat << std::endl;

			NumberOfTimeSteps = (time / dt_iter);

			NVT->setSimulationParameters(NumberOfTimeSteps, dt_iter, epsrParticle, epsfParticle, NrDimer, numberOfReplicas, initialDensity, numberOfEquilibrationSteps);


			for (int i_replica = 0; i_replica < numberOfReplicas; i_replica++){

				std::cout << "***********************\n starting trajectory for replica number " << i_replica + 1 << std::endl;



				NVT->InitialCondition();
			//	NVT->WritePositionsInFile();

#ifndef CHANGE_PARAMETERS_ON_THE_FLY
				NVT->PerformTrajectory();
#endif

#ifdef CHANGE_PARAMETERS_ON_THE_FLY
				NVT->PerformTrajectoryChangingParameters();
#endif
				averageValueReplica[i_replica] = NVT->getObservableForVariance(1);
				averageValueReplicaTwo[i_replica] = NVT->getObservableForVariance(2);
				averageValueReplicaThree[i_replica] = NVT->getObservableForVariance(3);

				//	std::cout<<"here"<<std::endl;


				samplesVector = NVT->getSamples();
				for (int i = 0; i < numberOfSamples; i++)	fileSamples << samplesVector[i] << "\t";
				fileSamples << NVT->getVar_mu() << "\t";
				fileSamples << dt_iter << "\t\n";

				double* averTimePerTimeStep = new double[100];
				averTimePerTimeStep = NVT->getArrayTimePerTimeStep();
				for (int i = 0; i < 100; i++) fileTimeMeasurings << averTimePerTimeStep[i] << "\t";
				fileTimeMeasurings << std::endl;


				//if(i_replica !=0)
				NVT->CleanSamples();

			}


			// -------------- compute variance


			double irepSum = 0;
			double irepSum2 = 0;
			double irepSum3 = 0;

			for (int irep = 0; irep < numberOfReplicas; irep++){


				irepSum = irepSum + averageValueReplica[irep];
				irepSum2 = irepSum2 + averageValueReplicaTwo[irep];
				irepSum3 = irepSum3 + averageValueReplicaThree[irep];

				//Observable.addSample(ObservableVariance(pParticles[irep]->getPosition()));

			}

			double variance = 0;
			double variance2 = 0;
			double variance3 = 0;

			for (int irep = 0; irep < numberOfReplicas; irep++){

				variance = variance + ((irepSum / (double)numberOfReplicas - averageValueReplica[irep])*(irepSum / (double)numberOfReplicas - averageValueReplica[irep]));
				variance2 = variance2 + ((irepSum2 / (double)numberOfReplicas - averageValueReplicaTwo[irep])*(irepSum2 / (double)numberOfReplicas - averageValueReplicaTwo[irep]));
				variance3 = variance3 + ((irepSum3 / (double)numberOfReplicas - averageValueReplicaThree[irep])*(irepSum3 / (double)numberOfReplicas - averageValueReplicaThree[irep]));
			}

			double var = NumberOfTimeSteps* variance / (double)numberOfReplicas;
			double var2 = NumberOfTimeSteps* variance2 / (double)numberOfReplicas;
			double var3 = NumberOfTimeSteps* variance3 / (double)numberOfReplicas;


			double standardDeviation = sqrt(var);
			double standardDeviation2 = sqrt(var2);
			double standardDeviation3 = sqrt(var3);

			std::cout << "\n ************** end replicas *************\n" << std::endl;



			if (numberOfReplicas > 1){
				std::cout << "dt=" << dt_iter << std::endl;
				std::cout << "Standard deviation of observable1 from IR is " << standardDeviation << std::endl;
				std::cout << "Standard deviation of observable2 from IR is " << standardDeviation2 << std::endl;
				std::cout << "Standard deviation of observable3 from IR is " << standardDeviation3 << std::endl;

				std::ofstream fileVariance("data/Variance");
				fileVariance << "stdDev=[" << standardDeviation << "\t" << standardDeviation2 << "\t" << standardDeviation3 << "];" << std::endl;
				fileVariance << "var=[" << var << "\t" << var2 << "\t" << var3 << "];" << std::endl;

			}



			dtArray[countdt] = dt_iter;
			varianceArray[countdt] = standardDeviation;
			varianceArray2[countdt] = standardDeviation2;
			//varianceArray3[countdt] = standardDeviation3;
			countdt++;

			//double tmpVar=varianceArray[countdt];

			dt_iter_count++;




			NVT->CleanSamples(); // clean averages objects - new dt is coming...

		}

	}
	std::cout << "\n ************* end dt loop *********" << std::endl;

	if (numberOfReplicas > 1){
		std::cout << "\n ***** Standard deviation for DimerDistance " << std::endl;
		for (int i = 0; i < NrOfDifferent_dt; i++) std::cout << "\nContin.Variance (IR) is " << varianceArray[i] * varianceArray[i] * dtArray[i] << "\t for dt=" << dtArray[i] << std::endl;
		double meanVariance = 0;
		for (int i = 0; i < NrOfDifferent_dt; i++) meanVariance = meanVariance + varianceArray[i] * varianceArray[i] * dtArray[i];
		double varianceIndt = meanVariance / NrOfDifferent_dt;
		std::cout << "Mean variance (IR) estimated from all dt " << varianceIndt << std::endl;
		std::cout << "Variance of DimerDistance from IR averaged in dt's is " << varianceIndt << std::endl;

		std::cout << "\n ***** Standard deviation for DimerPotential " << std::endl;
		for (int i = 0; i < NrOfDifferent_dt; i++) std::cout << "\n Contin.Variance (IR) is " << varianceArray2[i] * varianceArray2[i] * dtArray[i] << "\t for dt=" << dtArray[i] << std::endl;
		double meanVariance2 = 0;
		for (int i = 0; i < NrOfDifferent_dt; i++) meanVariance2 = meanVariance2 + varianceArray2[i] * varianceArray2[i] * dtArray[i];
		double varianceIndt2 = meanVariance2 / NrOfDifferent_dt;
		std::cout << "Mean variance (IR) estimated from all dt " << varianceIndt2 << std::endl;
		std::cout << "Variance of DimerPotential from IR averaged in dt's is " << varianceIndt << std::endl;

		fileSamples << "\n Variance_replicas_averaged_in_dt_DimerDistance = " << varianceIndt << "\n";
		fileSamples << "\n Variance_replicas_averaged_in_dt_DimerPotential = " << varianceIndt2 << "\n";


	}


	delete varianceArray;
	delete varianceArray2;
	delete dtArray;
	delete samplesVector;
	delete averageValueReplica;
	delete averageValueReplicaTwo;
	delete averageValueReplicaThree;

#endif
	//endif optimal_parameters


#endif



	std::cout<< "------------------- THE END ---------------------------"<<std::endl;


	//while(1);
}
