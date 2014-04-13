	//! String and io management modules
#include <iostream>
#include <fstream>
#include <regex>

//! Exception modules
#include <stdexcept>
#include <assert.h>
#include <signal.h>

//! Time module
#include <unistd.h>
#include "CH_Timer.H"

//! Math modules
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>

//! cGromCorr modules
#include "grom2atomFMO.h"

static bool keepRunning=true;

void intHandler(int dummy=0){
  std::cerr << "Received calcellation, exiting after current frame..."<<std::endl;
	keepRunning=false;
}


int main_Cr(std::string filename, std::string energyOutFileName, int n_samples, std::vector<std::string>& atomNameTable, std::vector<float> atomMassTable, std::vector<float> atomChargeTable){
	CH_TIMERS("MAIN");
	CH_TIMER("dE_compute",dE);
	CH_TIMER("FileRead",read);
	CH_TIMER("Correlate",corr);
	bool posIsDefined = false;
	std::streampos currentPos;
	bool resign = false;

	//bool corr_mtx_built = false;
	int dEsize = 0;
	int numGroups = 0;
	int numChromo = 0;
	int numTot = 0;

	// Delete the file that used to be at energyOutFile.time
	std::ofstream outputStream;
	outputStream.open(energyOutFileName+".time", std::ofstream::out);
	outputStream.close();

#ifdef DEBUG_ON
  std::cerr << "Number of samples: " << n_samples << "\tInternal keepRunning variable:"<< keepRunning <<std::endl;
	std::cout << "\t=================MergeNames present in the atom lookup table:======================" << std::endl;
	for (unsigned int i = 0 ; i < atomNameTable.size() ; i++){
		std::cout << atomNameTable[i] << " ";
	}
	std::cout << "\n\t=================End MergeName lookup table========================================" << std::endl;
	std::cout << std::endl;
#endif

	signal(SIGINT,intHandler);

	int timeCount = 0;
	std::vector< std::vector<float > > dE_sum_k_i(7);
	std::vector< std::vector<float > > dEsq_sum_k_ij(7); // row-major matrix
	bool corr_allocated = false;
	for (timeCount = 0 ; timeCount < n_samples && keepRunning ; timeCount++)
		{
		    std::stringstream this_time_stdout;
		    this_time_stdout << "--------------------------------------------------------------" << std::endl;
				this_time_stdout << "timeCount:" << timeCount << ", Press Ctrl-C to exit safely after this time." << std::endl;
				
		// --------------------------------- READ FILE -------------------------------------------
				
				std::vector<float > Xi;
				std::vector<std::string > atomTypei;
				std::vector<int > atomGroupi;
				float systemSize_nm = 0;
			
				std::cout << "Reading file..." <<std::endl;
		
				bool success = false;
				CH_START(read);
				while (!success){
					try{
						std::streampos outputPos;
						std::streampos* ptrPos = posIsDefined? &currentPos : 0;
						ReadOneTimeGrom2Atom(ptrPos,
															 	outputPos,
															 	Xi, 
															 	atomTypei,
															 	atomGroupi,
															 	systemSize_nm,
															 	filename,
                                true);
						currentPos = outputPos;
						posIsDefined = true;
						success = true;
					}
					catch (ReadError& e){
						if (resign){
							throw ReadError("Error reading file " + filename +", Fatal Error.");
						}
						std::cout << "Caught ReadError: " << e.what() << "\n 	Trying ReadOneTimeGrom2Atom again on file "<< filename <<std::endl;
						resign = true;
					}
					catch (InputFileMalformed& e){
						std::cout << "Caught InputFileMalformed: " << e.what() << "\n  Closing program. I should really tell you where, but I can't do that super well. The best you get is timeCount = "<<timeCount <<std::endl;
						break;
					}
				}
				CH_STOP(read);
		
				this_time_stdout << "Note: Press Ctrl-C at any time to exit safely after this frame."<<std::endl;
				this_time_stdout << "Looking up atomtypes..." << std::endl;
		
				// numGroups is the max element of the group numbers, and thus the number of distinct groups, which are assumed to start at 1 and go to numGroups.
				numGroups = *std::max_element(atomGroupi.begin(), atomGroupi.end());
				numChromo = std::abs(*std::min_element(atomGroupi.begin(), atomGroupi.end()));
				numTot = numGroups + numChromo;
				if (!corr_allocated)
					{
						for (int k = 0 ; k < 7 ; k++)
							{
								dE_sum_k_i[k].assign(numTot, 0);
								dEsq_sum_k_ij[k].assign(numTot*numTot, 0);
							}
						corr_allocated = true;
					}
			
				std::vector<float > cdc_kcal;
				this_time_stdout << "Computing site";
				for (int site = 1 ; site <= 7 ; site++)
					{
						CH_START(dE);
						this_time_stdout << " " << site;
		// --------------------------- CONVERT ATOMTYPE -------------------------------------------
		
						std::vector<float > atomMassi(  dEsize, 0);
						std::vector<float > atomSizei(  dEsize, 0);
						std::vector<float > atomChargei(dEsize, 0);
						std::vector<float > ChromoGroundCharges(0,0);

            // Matched arrays for *Table enter, and using information from site, atom[Mass/Size/Charge] and ChromoGroundCharges are generated.
            // Indexes for atom* match indexes for atomTypei and Xi.
            AtomDataLookup_v2(atomTypei, atomMassi, atomSizei, atomChargei, atomGroupi, site, ChromoGroundCharges, atomNameTable, atomMassTable, atomChargeTable);
		
		// ---------------------------- COMPUTE DENSITY -------------------------------------------
    
						ComputeCDC_v1(site, cdc_kcal, Xi, atomChargei, atomGroupi, ChromoGroundCharges);
            auto cm_kcal = [](float x){return x * 349.7;};
            std::vector<float > cdc_cm;
            cdc_cm.resize(cdc_kcal.size());
            std::transform(cdc_kcal.begin(), cdc_kcal.end(), cdc_cm.begin(), cm_kcal);
            
						double current_dEsum = std::accumulate(cdc_cm.begin(), cdc_cm.end(), 0.);
						this_time_stdout << "="<< current_dEsum << " cm-1";
						outputStream.open(energyOutFileName+".time", std::ofstream::out | std::ofstream::app);
						for (unsigned int i = 0 ; i < cdc_kcal.size() ; i++)
							{
								if (i != 0) outputStream << ",";
								outputStream << cdc_cm[i];
							}
						outputStream << std::endl;
						outputStream.close();
						CH_STOP(dE);

						CH_START(corr);
						
						int k = site - 1;	
						for (int i = 0 ; i < numTot ; i++)
							{
								dE_sum_k_i[k][i] += cdc_kcal[i];
								for (int j = 0 ; j < numTot ; j++){
									dEsq_sum_k_ij[k][(numTot * i) + j] += cdc_kcal[i] * cdc_kcal[j];
								}
							}
						CH_STOP(corr);
				}

				std::cout << this_time_stdout.str() << std::endl;

// --------------------------------- CORRELATE --------------------------------------
	}



	CH_START(corr);
	for (int k = 0 ; k < 7 ; k++){
		for (int i = 0 ; i < numTot ; i++){
			dE_sum_k_i[k][i] /= timeCount; 
			for (int j = 0 ; j < numTot ; j++){
				dEsq_sum_k_ij[k][(dEsize * i) + j] /= timeCount;
			}
		}
	}

	float corr_mtx[7][numTot][numTot];
	for (int k = 0 ; k < 7 ; k++){
		for (int i = 0 ; i < numTot ; i++){
			for (int j = 0 ; j < numTot ; j++){
				corr_mtx[k][i][j] = dEsq_sum_k_ij[k][(dEsize*i) + j] 
																				- dE_sum_k_i[k][i] * dE_sum_k_i[k][j];
			}
		}
	}

  std::cout << "Saving to file " << energyOutFileName << ".mtx" <<std::endl;
	std::string filename200(energyOutFileName+".mtx");
	std::ofstream f200;
	f200.open(filename200, std::ofstream::out);
	for (int k = 0 ; k < 7 ; k++){
		for (int i = 0 ; i < numTot ; i++){
			for (int j = 0 ; j < numTot ; j++){
				if (j != 0) f200 << ",";
				f200 << corr_mtx[k][i][j];
			}
			f200 << "\n";
		}
	}
	f200.close();

	CH_STOP(corr);
/*

// --------------------------------- RESCALE AND CORRECT AVERAGES --------------------------------------
	
	for (int boxCount = 0 ; boxCount < numGroups ; boxCount++){
		density_av[boxCount] /= n_samples;
	}
	
	float meanRMS = 0;
	float meanRMS_v2 = 0;
	for (int i = 0 ; i < numGroups ; i++){
		for (int j = 0 ; j < numGroups ; j++){
			corr_av[(i * numGroups) + j] /= n_samples; 
			prod_sum[(i * numGroups) + j] /= n_samples;
			prod_sum[(i * numGroups) + j] -= (density_av[i] * density_av[j]);
		}
		meanRMS 	 += corr_av [(i* numGroups) + i];
		meanRMS_v2 += prod_sum[(i* numGroups) + i];
	}
	meanRMS /= numGroups;


	for (int i = 0 ; i < numGroups ; i++){
		for (int j = 0 ; j < numGroups ; j++){
			corr_av[(i * numGroups) + j]  /= meanRMS;
			prod_sum[(i * numGroups) + j] /= meanRMS_v2;
		}
	}

*/
// --------------------------------- OUTPUT FILE -------------------------------------------

/*
	std::ofstream outputStream;
	std::string outputName;
	
	outputName = "CsvOut/DensityCorrAvLin.csv";
	outputStream.open(outputName, std::ofstream::out);
	for (int i = 0 ; i < density_av.size() ; i++){
		if (i!=0){
			outputStream << ",";
		}
		outputStream << density_av[i];
	}
	outputStream << std::endl;
	for (int i = 0 ; i < numGroups;  i++){
		for (int j = 0 ; j < numGroups ; j++){
			if (j!=0){
				outputStream << ",";
			}
			outputStream << corr_av[(i * numGroups) + j];
		}
		outputStream << std::endl;
	}
	outputStream << std::endl;
	for (int i = 0 ; i < numGroups;  i++){
		for (int j = 0 ; j < numGroups ; j++){
			if (j!=0){
				outputStream << ",";
			}
			outputStream << prod_sum[(i * numGroups) + j];
		}
		outputStream << std::endl;
	}
	outputStream << std::endl;
	//std::cout << std::endl;
	outputStream.close();
	std::cout << "Wrote the linear stream, density and correlation" << std::endl;
*/



	return 0;
}






int main(int argc, char** argv){
  ConfigData cfg;

#ifdef OLDPARSER
	std::string cfgFileName;
	if (argc == 1){
		cfgFileName = "configFMO.cfg";
		std::cout << "Using default config file " << cfgFileName << std::endl;
	}
	else if (argc == 2){
	 	cfgFileName = argv[1];
	}
	else{
		throw std::invalid_argument("Usage: ./parse configFile");
	}
	int rtn = ReadFMOConfig(cfgFileName, cfg)
	if (rtn == EXIT_FAILURE){
		return(EXIT_FAILURE);
	}

#else
  if (argc != 6){
    std::cout << "Six inputs required for this program; have you considered using a python script to control this program??" << std::endl;
    std::cout << "Order of args: fmo_traj_path, fmo_top_path, bcx_itp_path, output_path, num_frames" << std::endl;
    return(EXIT_FAILURE);
  }
  else{
    cfg.fmo_traj_path = argv[1];
    cfg.fmo_top_path  = argv[2];
    cfg.bcx_itp_path  = argv[3];
    cfg.output_path  = argv[4];
    cfg.num_frames    = std::atoi(argv[5]);
  }
#endif


	std::vector<std::string> atomNames;
	std::vector<float> atomMasses;
	std::vector<float> atomCharges;
		

/*! 
 *  Uses grom2atomFMO.cpp: ReadFMOConfig to work through the config file.
 * 		Parses the topology files and puts them into paired list [atomNames, atomMass, atomCharge]
 * 			ParseTopology gets called twice:
 * 				> Once to get the difference between ground and excited charges from the bchlFile
 * 				> Again to get the topology file masses and charges
 * 			atomNames are stored as "molName,atomName" (e.g. BCL,CHA) (stores atomName, not atomType)
 **/

    
  try{
	  ParseTopology(cfg.bcx_itp_path, atomNames, atomMasses, atomCharges, Excited::YES);
  }
  catch (InputFileMinorMalformed &mal){
    std::cout << "Minor issues found with input file'" << cfg.bcx_itp_path << "'" << std::endl;
  }
  try{
	ParseTopology(cfg.fmo_top_path, atomNames, atomMasses, atomCharges);
  }
  catch (InputFileMinorMalformed &mal){
    std::cout << "Minor issues found with input file '" << cfg.fmo_top_path << "'" << std::endl;
  }


/*! 
 *  Uses an executable from this module to call, from grom2atomFMO.cpp:
 *  	> ReadOneTimeGrom2Atom to load a single frame starting at location "streampos"
 *  			. if streampos input is NULL, then it starts at the beginning
 *  	> AtomDataLookup to convert .gro file into std::vector of charges and other data from topology
 *  	> ComputeCDC to compute the energies for the frame
 **/
  std::cout << "Running main_Cr..." << std::endl;
	return main_Cr(cfg.fmo_traj_path, cfg.output_path, cfg.num_frames,  atomNames, atomMasses, atomCharges);
}
