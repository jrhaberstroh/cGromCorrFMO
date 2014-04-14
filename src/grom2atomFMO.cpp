// String and io management modules
#include <iostream>
#include <fstream>
#include <regex>
#include <boost/filesystem.hpp>
#include <libconfig.h++>

// Exception modules
#include <stdexcept>
#include <assert.h>

// Math modules
#include <cmath>
#include <algorithm>
#include <vector>

// module to load sleep(float seconds)
#include <unistd.h>
#include "grom2atomFMO.h"

//! \brief Class to assist in interpreation of strings during Topology and .gro file reading
class IsChars
{
	public:
//! \brief Ctor for IsChars class, to mark characters with a boolean if they are in charsToRemove 
	IsChars(const char* charsToRemove) : chars(charsToRemove) {};
//! \brief Asserts whether input character is in charsToRemove
	bool operator()(char c)
	{
		for(const char* testChar = chars; *testChar != 0; ++testChar)
		{
			if(*testChar == c) { return true; }
		}
		return false;
	}

	private:
		const char* chars;
};


const char* EXC_STRING = "EXC";






/////////////////////////////////////////////////
/// ParseTopology
/////////////////////////////////////////////////
/*! \brief 	Loads topology file into a lookup table
 *
 *  \param 		filename			 INPUT Name of topology file
 *
 *  \param  	mergeNameTable OUTPUT Table with "merge names", combination of MOL,ATM
 *  \param		massTable			 OUTPUT Table with atom masses
 *  \param		chargeTable		 OUTPUT Table with atom charges
 *  \param		xc						 OUTPUT Use excited charge lookup (i.e. input dq = qE - qO in the merge name (MOL,ATM + ',EXC') and q0 into MOL,ATM)
 *
 *  Adds data from topology to the table of coulomb forcefield values. Currently only used for
 *  charge computations because of the simplicity of their interactions and because of the
 *  relevance to computing molecular energy gaps.
 *
 *  Names are stored in a MOL,ATM format, where ATM is the atomName, not the atomType. When xc
 *  is on, charges diferences between ground and excitation are computed and entered into the table instead
 *  of ground state or excited state charges. An additional entry is also added, if the entry is accepted,
 *  which is of the form MOL,ATM,EXC. EXC here is literally the string "EXC". This allows lookup of
 *  excited state charges and ground state charges in the same computation.
 *
 *  Names are accepted or rejected based on whether another entry with the same mergeName exists. Thus,
 *  if either the MOL or the ATM are different, the entry is accepted.
 *
 *  Runs recursively on successive #includes within the topology files, up to arbitrary depth.
 */

// REQUIRES: empty massTable, empyt chargeTable
void ParseTopology(std::string const& filename,
									 std::vector<std::string>& mergeNameTable, 
									 std::vector<float>& massTable,
									 std::vector<float>& chargeTable,
									 Excited xc){
	
	//std::cout << "\tEntering function ParseTopology for " << filename << std::endl;

	std::ifstream atomFstream;
	atomFstream.open(filename);
	std::string line;
  bool returnWarning = false;
  std::stringstream warningMessage;
	bool atomFlagOn = false;
	bool firstAtomLine = false;
	if (!atomFstream.is_open()){
		throw InputFileMalformed("Requested topology file " + filename + " could not be opened.");
	}
	while (getline(atomFstream, line)){
		// Treat semicolons as comment delimiters
		line.assign(line.substr(0,line.find(";")));

		std::istringstream iss(line);
		std::string str1;
		std::string str2;
		std::string str3;

		iss >> str1 >> str2 >> str3;

		if (str1.compare("[") == 0 && str3.compare("]") == 0){
			if (str2.compare("atoms") == 0){
				//std::cout << "\t\tTurned atom flag ON  at line: " <<  line << std::endl;
				atomFlagOn = true;
				firstAtomLine = true;
				sleep(.1);
			}
			else if (atomFlagOn){
				//std::cout << "\t\tTurned atom flag OFF at line: " <<  line << std::endl;
				atomFlagOn = false;
				sleep(.1);
			}
		}

    if (str1.compare("#include") == 0 && xc == Excited::NO){
			str2.erase( std::remove(str2.begin(), str2.end(), '"'), str2.end() );
			std::string parentPath = boost::filesystem::path(filename).parent_path().string();
			ParseTopology(parentPath + "/" + str2, mergeNameTable, massTable, chargeTable);
			//std::cout << "\tReturning to ParseTopology on file "<<filename <<std::endl;
		}

    if (atomFlagOn && !firstAtomLine && line.compare("") != 0){
			std::istringstream iss(line);
			std::string atomName;
			std::string molName;
			std::string trash;
			float atomCharge;
			float atomMass;
			float atomChargeEX;
			float atomMassEX;
      int atomID;
	
			if (xc == Excited::YES){
				//std::cout << line << std::endl;
        std::vector<std::string> line_arr(0);
        do
        {
            std::string substr;
				    iss >> substr;
            line_arr.push_back(substr);
        } while (iss);

        if (line_arr.size() != 12 && line_arr.size() != 9){
            std::stringstream err_msg;
            err_msg << "Bad Topology file, found wrong number of line entries for FEP topology (expected either 8 or 11 and received <"<<line_arr.size()<<">";
            err_msg << "\nBad line = "<<line;
            throw InputFileMalformed(err_msg.str());
        }
        atomID  = std::atoi(line_arr[0].c_str());
        molName = line_arr[3];
        atomName = line_arr[4];
        atomCharge   = std::atof(line_arr[6 ].c_str());
        atomMass     = std::atof(line_arr[7 ].c_str());
        if (line_arr.size() == 12){
            atomChargeEX = std::atof(line_arr[9 ].c_str());
            atomMassEX   = std::atof(line_arr[10].c_str());
        }
        else{
            atomChargeEX = atomCharge;
            atomMassEX   = atomMass  ;
        }

        if (atomMassEX != atomMass){
            std::stringstream err_msg;
            err_msg << "Bad Topology file, atom in " << molName << " has bad masses M"<< atomMassEX<<"and Mex=" << atomMass;
            err_msg << "\nBad line: "<<line;
            throw InputFileMalformed(err_msg.str());
        }
				atomChargeEX = atomChargeEX - atomCharge;
				// Use only the relative charge to compute CDC
			}
			else{
				iss >> atomID >> trash >> trash >> molName >> atomName >> trash >> atomCharge >> atomMass;
			}

			std::stringstream mergeName;
      mergeName << molName << "," << atomName << atomID;
			std::stringstream mergeNameEX;
      mergeNameEX << molName << "," << atomName << atomID << "," << EXC_STRING;

			bool atomIncluded = false;
			for (int i = 0 ; i < mergeNameTable.size() ; i++){
				if (mergeNameTable[i].compare(mergeName.str()) == 0){
					if(massTable[i] != atomMass){
            returnWarning = true;
            warningMessage << "Mass of " << mergeNameTable[i] << " conflicts with prior entry" <<std::endl;
#if DEBUG_ON
						std::cerr << "Note: " << mergeNameTable[i] << ", mass " << massTable[i] << " conflicts with "<< mergeName.str() << ", mass " <<atomMass << std::endl;
#endif
					}
					if (chargeTable[i] != atomCharge){
            returnWarning = true;
            warningMessage << "Charge of " << mergeNameTable[i] << " conflicts with prior entry" <<std::endl;
#if DEBUG_ON
						std::cerr << "Note: " << mergeNameTable[i] << ", charge " << chargeTable[i] << " conflicts with "<< mergeName.str() << ", charge" <<atomCharge << std::endl;
#endif
					}
					atomIncluded = true;
					break;
				}
			}

			if (!atomIncluded && atomName.compare("") != 0){
				mergeNameTable.push_back(mergeName.str());
				massTable.push_back(atomMass);
				chargeTable.push_back(atomCharge);

				if (xc == Excited::YES) {
					mergeNameTable.push_back(mergeNameEX.str());
					massTable.push_back(atomMassEX);
					chargeTable.push_back(atomChargeEX);
				}
			}
#if DEBUG_ON
      else if (atomIncluded){
        std::cerr << "Note: found duplicate entries on "<<mergeName.str()<<std::endl;
      }
#endif
		}

		if (firstAtomLine){
			firstAtomLine =false;
		}
	}
  if (returnWarning){
    throw InputFileMinorMalformed(warningMessage.str());
  }
}





#ifdef OLDPARSER
int ReadFMOConfig(std::string const& cfgFileName, Config){

//-------------------------- LOADING CFG FILE  -------------------------------
	libconfig::Config cfg;
	try
  {
    cfg.readFile(cfgFileName.c_str());
  }
  catch(const libconfig::FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const libconfig::ParseException &pex)
  {
    std::cerr << "Warning: error while reading file: " << cfgFileName.c_str() << std::endl;
    std::cerr << pex.what() << std::endl;
    return(EXIT_FAILURE);
  }

//-------------------------- .gro TRAJ FILE  -------------------------------
  
	try
  {
		std::string name = cfg.lookup("fmoTraj");
		trajFileName = name;
		std::cout << "Traj filename: " << name << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex)
  {
    std::cerr << "No 'fmoTraj' setting in configuration file." << std::endl;
  }
  
//-------------------------- NUM SAMPLES -------------------------------
//
	try
  {
		int numSamples = cfg.lookup("numSamples");
		n_entries = numSamples;
		std::cout << "Number of samples to use: " << numSamples << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex)
  {
    std::cerr << "No 'numSamples' setting in configuration file." << std::endl;
  }

//-------------------------- TOPOLOGY -------------------------------

	std::string topologyFile = "";
	try
  {
		std::string topologyFile2 = cfg.lookup("topology");
		topologyFile = topologyFile2;
		std::cout << "Looking up topology: " << topologyFile << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex)
  {
    std::cerr << "No 'topology' setting in configuration file." << std::endl;
  }

//-------------------------- BCHL TOPOLOGY -------------------------------
//
	std::string bchlFile = "";
	try
  {
		std::string bchlFile2 = cfg.lookup("bchlFile");
		bchlFile = bchlFile2;
		std::cout << "Looking up bchlFile: " << bchlFile << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex)
  {
    std::cerr << "No 'bchlFile' setting in configuration file." << std::endl;
  }


//-------------------------- OUTPUT PATH ---------------------------------
//
	std::string energyFile = "";
	try
  {
		std::string energyFile2 = cfg.lookup("energyFileOut");
		energyFile = energyFile2;
		std::cout << "Looking up energyFileOut: " << energyFile << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex)
  {
    std::cerr << "No 'energyFileOut' setting in configuration file." << std::endl;
  }
	outFileName = energyFile;
	
}
#endif

/////////////////////////////////////////////////
// ReadOneTimeGrom2Atom
/////////////////////////////////////////////////

void ReadOneTimeGrom2Atom(std::streampos const * const  currentPos, std::streampos & outPos, std::vector<float > & Xi, std::vector<std::string > & atomTypei, std::vector<int>& atomGroupi, float & systemSize_nm, std::string const& fileName, bool should_print){

	std::ifstream gromFile;
	gromFile.open(fileName);
	// Check for nullpointer
	if (currentPos != 0){
		gromFile.seekg(*currentPos);
	}
	
	std::string line;
	int atomTotal = NAN;
	int atomCount = 0;
	bool encounteredHeader = false;
	int lastMoleculeNum = 0;
	int groupNumberMinus = 0;
	int groupNumberPlus = 0;
  int atomID_prot = 1;
  int atomID_bcl = 1;

	if (gromFile.is_open()){
		while (getline(gromFile,line)){
			
			std::regex time_regex(".*t=.*");
			std::regex atom_regex(".*SOL.*");
	// Description of flow of function, steps 1-4:
	// 1 : Read the first two lines. The first tells you the time, the second sets atomTotal to the number of atoms
	// 4 : Assert that Xi filled up appropriately, and that atomCount has reached atomTotal, then break.
	
			if (std::regex_match(line, time_regex)){
				if (encounteredHeader){
					assert (std::none_of(Xi.begin(), Xi.end(), [](float x){return std::isnan(x);}));
					assert (atomCount == atomTotal);
					break;
				}
				// When you encounter the time, take the following line as the number of atoms in the timeslice
				encounteredHeader = true;
				if (should_print) std::cout << line << std::endl;
				if (getline(gromFile,line)){
					try{
						atomTotal = Convert<int>(line);
					}
					catch(std::invalid_argument & e){
						throw InputFileMalformed("Encountered non-integer entry when number-of-atoms was excpected after time-header line.");
					}
				}
				else{
					throw InputFileMalformed("Encountered EOF when number-of-atoms was expected after time-header line.");
				}

				Xi.assign(atomTotal*3,nan(""));
				atomTypei.assign(atomTotal,"");
				atomGroupi.assign(atomTotal,NAN);
			}
	
	// 2 : Add atom types to the atom array and positions to the position array
			
			
			else if (atomCount < atomTotal && encounteredHeader){
				std::istringstream iss(line);

				int loop_count = 0;
				std::string molName;
				std::string atomName;
        int atomID;
				int groupNumber;

				do{
					std::string sub;
					iss >> sub;
					if (loop_count == 0){
            // For the first column, go to the first character that is non-numeric to grab the molecule name
						int i = 0;
						std::string digits("1234567890");
						while (std::any_of(digits.begin(), digits.end(), [i, sub](char j){ return j == sub[i]; } )){
							i++;
						}
						int thisMoleculeNum = std::atoi(std::string(sub.begin(), sub.begin() + i).c_str());
						std::string thisMolName(sub.begin()+i, sub.end());
						
						if (thisMolName.compare("BCL") == 0
#ifdef BCX_VALID
                || thisMolName.compare("BCX") == 0
#endif
                ){
              thisMolName = "BCX";
							if (thisMoleculeNum != lastMoleculeNum){
								groupNumberMinus--;
                atomID_bcl = 1;
								if (should_print) std::cout << thisMolName << ":" << groupNumberMinus << "\t";
							}
							groupNumber = groupNumberMinus;
							lastMoleculeNum = thisMoleculeNum;
              atomID = atomID_bcl;
              atomID_bcl++;
						}
						else{
							if (thisMoleculeNum != lastMoleculeNum){
								groupNumberPlus++;
								if (should_print) std::cout << thisMolName << ":" << groupNumberPlus << "\t";
							}
							groupNumber = groupNumberPlus;
							lastMoleculeNum = thisMoleculeNum;
              atomID = atomID_prot;
              atomID_prot++;
						}
						
						molName = thisMolName;
					}
					if (loop_count == 1) atomName = sub;
					//std::cout << "Substring: " << sub << std::endl;
					if (loop_count == 3)  Xi[atomCount*3 + 0] = std::atof(sub.c_str());
					if (loop_count == 4)  Xi[atomCount*3 + 1] = std::atof(sub.c_str());
					if (loop_count == 5)  Xi[atomCount*3 + 2] = std::atof(sub.c_str());
					//if (loop_count == 6)  Vi[atomCount][0] = sub;
					//if (loop_count == 7)  Vi[atomCount][0] = sub;
					//if (loop_count == 8)  Vi[atomCount][0] = sub;
					loop_count++;
				} while(iss);

        assert(atomID != -1);
        std::stringstream this_atom_type;

        this_atom_type << molName<<","<<atomName << atomID;

				atomTypei[atomCount] = this_atom_type.str();
				atomGroupi[atomCount] = groupNumber;

				//std::cout << Xi[atomCount][0] << std::endl;
				atomCount++;
			}

	// 3 : Read the sytem size off of the bottom of a time slice
	//
	// WARNING: This section works by chance, because the input format is very spectific. Search here if bugs arise.
			else if (encounteredHeader){
				if (should_print) std::cout << line << std::endl;

				// Extract system size
				std::istringstream iss(line);
				std::string sub;
				iss >> sub;
				systemSize_nm = std::atof(sub.c_str());

				// Set stream pos after this line
				std::streampos outPos_try;
				outPos_try = gromFile.tellg(); //Can throw error
				outPos = outPos_try;
			}
			else {
				throw ReadError("Mysterious line appeared!");
			}
	// End switch 
		}
		gromFile.close();
		if (should_print) std::cout<<std::endl;
	}
	//If file opening failed:
	else{
		throw ReadError("Open failed on file passed in.");
	}
}



/////////////////////////////////////////////////
/// AtomDataLookup
/////////////////////////////////////////////////
/*! \brief 	Uses lookup tables from the topology file (ParseTopology) to convert .gro data (ReadOneTimeGrom2Atom) into charges and atomtypes
 */


void AtomDataLookup(std::vector<std::string > const& atomNamei, 
		std::vector<float > & atomMassi_amu, 
		std::vector<float > & atomSizei_nm,
		std::vector<float > & atomChargei_e,
		std::vector<int > const & atomGroupi,
		int excitationSite,
	 	std::vector<std::string > const& atomTypeTable, 
		std::vector<float > const & atomMassTable,
		std::vector<float > const & atomChargeTable){

	atomMassi_amu.assign(atomNamei.size(), 0);
	atomSizei_nm.assign(atomNamei.size(), 0);
	atomChargei_e.assign(atomNamei.size(), 0);
	// Loop atoms i, then loop the lookup table
	for (int i = 0 ; i < atomNamei.size() ; i++){
		std::string thisAtomName = atomNamei[i];
		if (atomGroupi[i] == -excitationSite){
			thisAtomName = thisAtomName +","+EXC_STRING;
		}
		bool found = false;
		for (int j = 0 ; j < atomTypeTable.size() ; j++){
			//std::cout << "Looking for regex match with " << thisAtomName << " to " << atomTypeTable[j]<< "..." << std::endl;
			if (std::regex_match(thisAtomName, std::regex(atomTypeTable[j]))){
				//std::cout << "\tMatch found on " << thisAtomName << std::endl;
				atomMassi_amu[i] = atomMassTable[j];
				atomSizei_nm[i] = .3;
				atomChargei_e[i] = atomChargeTable[j];
				found = true;
				break;
			}
		}
		if (!found){
			throw InputFileMalformed("Atom name " + thisAtomName + " from input file not found in lookup table.");
		}
	}
}

#include <map>
typedef std::map<std::string, int> IndexMap;

void AtomDataLookup_v2(
		std::vector<std::string > const& atomNamei, 
		std::vector<float > & atomMassi_amu, 
		std::vector<float > & atomSizei_nm,
		std::vector<float > & atomChargei_e,
		std::vector<int > const & atomGroupi,
		int excitationSite,
	 	std::vector<float > & atomGroundCharges_excitedGroup,
	 	std::vector<std::string > const& atomTypeTable, 
		std::vector<float > const & atomMassTable,
		std::vector<float > const & atomChargeTable){

	atomMassi_amu.assign(atomNamei.size(), 0);
	atomSizei_nm.assign(atomNamei.size(), 0);
	atomChargei_e.assign(atomNamei.size(), 0);

	IndexMap indexMap;
	for (int j = 0 ; j < atomTypeTable.size() ;j++){
		indexMap.insert(std::pair<std::string,int> (atomTypeTable[j], j));
	}

	atomGroundCharges_excitedGroup.clear();
	bool excited;
	std::string excitedAtom_groundName;
	// Loop atoms i, then loop the lookup table
	for (int i = 0 ; i < atomNamei.size() ; i++){
		std::string thisAtomName = atomNamei[i];

    // "excited" will activate when the atomGroup integer (1..7) is the excitation site -(1..7)
    // When it matches, it will append EXC_STRING to find the excited-state topology values
		excited = (atomGroupi[i] == -excitationSite);
		if (excited){
			excitedAtom_groundName = thisAtomName;
			thisAtomName = thisAtomName +","+EXC_STRING;
		}

		IndexMap::iterator it;

		it = indexMap.find(thisAtomName);
		if (it != indexMap.end()){
			int j = it->second;
			atomMassi_amu[i] = atomMassTable[j];
			atomSizei_nm[i] = .3;
			atomChargei_e[i] = atomChargeTable[j];
		}
		else {
			throw InputFileMalformed("Atom name " + thisAtomName + " from input file not found in lookup table.");
		}
		if (excited)
		{
			it = indexMap.find(excitedAtom_groundName);
			if (it != indexMap.end()){
				int j = it->second;
				atomGroundCharges_excitedGroup.push_back(atomChargeTable[j]);
				//std::cout << "Pushing " << atomChargeTable[j] << " for " << atomTypeTable[j] << " ground charge." << std::endl;
			}
			else {
				throw InputFileMalformed("Atom name " + excitedAtom_groundName + " from input file not found in lookup table.");
			}
		}
	}
}

// ComputeCDC_v1
/* \brief 	Uses input from topology and .gro (AtomDataLookup) to compute the energy vector (partitioned by molecule) with a particular site excited
 * 	
 * 	\input		chromoSite	Site to compute dE for
 *						Xi					Position of atom i [nm]
 *						qi					Charge of atom i (dQ for chromo, Q for others) [e]
 *						atomGroups	
 *
 *	\output		cdc_kcal		Output of cdc computation
 *
 */


void ComputeCDC_v1( int chromoSite, 
										std::vector<float > & cdc_kcal,
										std::vector<float >const& Xi,
										std::vector<float >const& qi,
										std::vector<int > const & atomGroups,
										std::vector<float > const& excitedAtom_groundCharges){
	// Assert that site is a valid entry of atomGroups; that there are entries in atomGroup of -chromoSite
	assert(std::any_of(atomGroups.begin(), atomGroups.end(), [chromoSite](int i){return i == -chromoSite;}));
	int nGroups = *std::max_element(atomGroups.begin(), atomGroups.end());
	int nChromo = std::abs(*std::min_element(atomGroups.begin(), atomGroups.end()));


	cdc_kcal.assign(nGroups + nChromo,0);
	std::vector<int > chromoInd(0,0);
	for (unsigned int i = 0 ; i < atomGroups.size() ; i++){
		if (atomGroups[i] == -chromoSite){
			chromoInd.push_back(i);
		}
	}
	
	bool printExAtoms = false;
	if (printExAtoms)
		{
			for (unsigned int i = 0 ; i < chromoInd.size() ; i++)
				{
					int ind = chromoInd[i];
					float x = Xi[3 * ind + 0];
					float y = Xi[3 * ind + 1];
					float z = Xi[3 * ind + 2];
					float q0 = excitedAtom_groundCharges[i];
					float dq = qi[ind];
					float qE = q0 + dq;
					std::cout << x << "\t"
										<< y << "\t"
										<< z << "\t"
										<< q0<< "\t"
									 	<< qE<< std::endl;
				}
		}

	float esConst_kCalnm_e2 = 33.2;
	assert(excitedAtom_groundCharges.size() == chromoInd.size() || excitedAtom_groundCharges.size() == 0);
	// Natural counting (one-based indexing) for counting atomGroups
	for (int i = 1 ; i <= nGroups ; i++){
		//bool printed_atm = false;
		for (int atm = 0 ; atm < atomGroups.size() ; atm++)
		{
      // Include all atomgroups except the current excited-site
			if (abs(atomGroups[atm]) == i && -atomGroups[atm] != chromoSite)
			{
				for (int chromoAtm = 0 ; chromoAtm < chromoInd.size() ; chromoAtm ++)
				{
					float dist_temp = 0.;
					for (int dir = 0 ; dir < 3 ; dir++)
					{
						float dx = Xi[3 * atm + dir] - Xi[3 * chromoInd[chromoAtm] + dir];
						dist_temp += dx * dx;
					}
					dist_temp = std::sqrt(dist_temp);
					//if (dist_temp < .2){ std::cout << dist_temp <<std::endl; }
					// Value of q0_i * qE_j - q0_i * q0_j (only includes one chromo atom)
          float coupling = esConst_kCalnm_e2 / dist_temp * (qi[atm] * qi[chromoInd[chromoAtm]])  ;
          if (atomGroups[atm] < 0){
					  cdc_kcal[i -1] +=  coupling;
          }
          else if (atomGroups[atm] > 0){
					  cdc_kcal[nChromo + i -1] += coupling;
          }
				}
			}
#ifdef INTRACHROMO_COUPLE
			if (-atomGroups[atm] == i && i == chromoSite)
			{
				//if (!printed_atm) { std::cout << "Computing for site " << i << std::endl; printed_atm = true;}
				if (excitedAtom_groundCharges.size() != 0)
				{
					float q0 = excitedAtom_groundCharges[excited_atom_counter];
					assert(atm == chromoInd[excited_atom_counter]);
					for (int chromoAtm = 0 ; chromoAtm < chromoInd.size() ; chromoAtm ++)
					{
						if (chromoAtm != excited_atom_counter)
						{
							// atm is the list index of excited_atom_counter, analogous to chromoInd[chromoAtm]
							//std::cout << atm << ":"<<excited_atom_counter <<"; " << chromoInd[chromoAtm] << ":"<<chromoAtm << std::endl;
							int ind_c1 = atm;
							int ind_c2 = chromoInd[chromoAtm];
							float dist_temp = 0.;
							for (int dir = 0 ; dir < 3 ; dir++){
								float dx = Xi[(3 * ind_c1) + dir] - Xi[(3 * ind_c2) + dir];
								dist_temp += dx * dx;
							}
							assert(dist_temp != 0);
							dist_temp = std::sqrt(dist_temp);
							//if (dist_temp < .11) {std::cout << dist_temp <<std::endl;}
							// There are cross terms from  qE_i * qE_j - q0_i*q0_j in the substitution for qE_i = q0_i + dq_j because we have two chromophore charges in this computation.
							// Because we iterate over all of the sites, include a factor of 1/2 to avoid double counting, and include only half of the overall contribution between the sites.
							float q0_1  = excitedAtom_groundCharges[excited_atom_counter];
							float dq_1  = qi[ind_c1];
							float qE_1  = q0_1 + dq_1;
							float q0_2  = excitedAtom_groundCharges[chromoAtm];
							float dq_2  = qi[ind_c2];
							float qE_2  = q0_2 + dq_2;
							//cdc_kcal[nChromo - i] += cm_kCal * esConst_kCalnm_e2 / dist_temp * (dq_2 * q0 + .5 * dq_1 * dq_2) ;
							cdc_kcal[nChromo - i] += esConst_kCalnm_e2 / dist_temp * .5 * (qE_1 * qE_2 - q0_1 * q0_2) ;
						}
					}
					excited_atom_counter++;
				}
			}
#endif
		}
	}
	
}




float DensityContributionErf(float dist_nm[3], float systemSize_nm, float halfBin_nm, float atomSize_nm, float CUTOFFsq_nm2){
	//Compute nearest periodic neighbors
	float halfBox_nm = systemSize_nm / 2;
	float dsq_nm2 = 0;

	for (int hat = 0 ; hat < 3 ; hat++){
		dist_nm[hat] = std::fmod(dist_nm[hat], systemSize_nm);

		if (dist_nm[hat] < -halfBox_nm){
			dist_nm[hat] += systemSize_nm;
		}
		if (dist_nm[hat] > halfBox_nm){
			dist_nm[hat] -= systemSize_nm; 
		}

		dsq_nm2 = dist_nm[hat] * dist_nm[hat];
		if (dsq_nm2 >= CUTOFFsq_nm2){
			return 0;
		}
	}
	// Compute the contribution to the density field
	float densityContribution = 1.;
	float halfBin_red = halfBin_nm / atomSize_nm;

	for (int hat = 0 ; hat < 3 ; hat++){

		float absDist_red = std::abs(dist_nm[hat])/atomSize_nm;

		densityContribution *= - std::erf(-halfBin_red - absDist_red ) 
												 	 + std::erf(+halfBin_red - absDist_red );
		
		//if (atomCount < 10 && boxCount < 10){
		//	std::cout << "absDist_red["	<< hat << "]: "<<  absDist_red << ", " << densityContribution << "\t";
		//}
		
	}
	return densityContribution / 8.;

}

void ComputeErfDensity_v1(std::vector<float > const & Xi_nm,
													std::vector<float > const & atomMassi_amu,
													std::vector<float > const & atomSizei_nm,
													float systemSize_nm,
													int Nb,
													std::vector<float > & density_amu){

	density_amu.assign(Nb*Nb*Nb, 0);
	double CUTOFF_nm = 1;
	double CUTOFFsq_nm2= CUTOFF_nm * CUTOFF_nm;

	float halfBox_nm = systemSize_nm / 2.;
	float halfBin_nm = halfBox_nm / Nb;


	std::cout << "SystemSize (nm), HalfSystem (nm):"
						<< systemSize_nm << "," << halfBox_nm << std::endl ;


	for (int atomCount = 0 ; atomCount < Xi_nm.size()/3 ; atomCount++){
		for (int boxCount = 0 ; boxCount < density_amu.size() ; boxCount++){
			float dist_nm[3] = {};

			dist_nm[0] = (float(boxCount % Nb)      * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 0];
			dist_nm[1] = (float(boxCount / Nb % Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 1];
			dist_nm[2] = (float(boxCount / Nb / Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 2];
			
			float densityContribution = DensityContributionErf(dist_nm, systemSize_nm, halfBin_nm, atomSizei_nm[atomCount],CUTOFFsq_nm2);
			density_amu[boxCount] +=  densityContribution * atomMassi_amu[atomCount];
		}
	}
}





void ComputeErfMomentum_v1(std::vector<float > const & Xi_nm,
													 std::vector<float > const & Vi_nm_ps,
													 std::vector<float > const & atomMassi_amu,
													 std::vector<float > const & atomSizei_nm,
													 float systemSize_nm,
													 int Nb,
													 std::vector<float > & momentum_amunm_ps){

	momentum_amunm_ps.assign(Nb*Nb*Nb*3, 0);
	double CUTOFF_nm = 1.;
	double CUTOFFsq_nm2= CUTOFF_nm * CUTOFF_nm;

	float halfBox_nm = systemSize_nm / 2.;
	float halfBin_nm = halfBox_nm / Nb;


	std::cout << "SystemSize (nm), HalfSystem (nm):"
						<< systemSize_nm << "," << halfBox_nm << std::endl ;


	for (int atomCount = 0 ; atomCount < Xi_nm.size()/3 ; atomCount++){
		for (int boxCount = 0 ; boxCount < Nb*Nb*Nb ; boxCount++){
			float dist_nm[3] = {};

			dist_nm[0] = (float(boxCount % Nb)      * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 0];
			dist_nm[1] = (float(boxCount / Nb % Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 1];
			dist_nm[2] = (float(boxCount / Nb / Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 2];
			
			float densityContribution = DensityContributionErf(dist_nm, systemSize_nm, halfBin_nm, atomSizei_nm[atomCount],CUTOFFsq_nm2);
			densityContribution *= atomMassi_amu[atomCount];

			for (int hat = 0 ; hat < 3 ; hat++){
				momentum_amunm_ps[boxCount*3 + hat] += densityContribution * Vi_nm_ps[atomCount*3+ hat];
			}

		}
	}
}




/*
void AtomLookup::AddAtom(std::string atomName, float atomMass_amu, float atomSize_nm){
	m_atomNames.push_back(atomName);
	m_atomMass_amu.push_back(atomMass_amu);
	m_atomSize_nm.push_back(atomSize_nm);
	assert(m_atomNames.size() == m_atomMass_amu.size() && m_atomNames.size() == m_atomSize_nm.size());
}

void AtomLookup::GetData(std::string const& atomNameIn, float & massOut_amu, float & sizeOut_nm){
	int i = 0;
	while (i < m_atomNames.size() && ! atomNames[i].compare(atomNameIn)){
		i++;
	}

	if (i == m_atomNames.size()){
		throw InputFlieMalformed("Atom name " << atomNameIn << " from input file not found in lookup table.")
	}
	massOut_amu = m_atomMass_amu[i];
	sizeOut_nm  = m_atomSize_nm[i];
}

std::vector<float > AtomLookup::ConvertNameToMass(std::vector<std::string > & atomNames){
	std::vector<float > masses(atomNames.size(), 0);
	return masses;
}

std::vector<float > AtomLookup::ConvertNameToSize(std::vector<std::string > & atomNames){
	std::vector<float > sizes(atomNames.size(), 0);
	return sizes;
}
*/

