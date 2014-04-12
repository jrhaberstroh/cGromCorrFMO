#ifndef __GROM2ATOMFMO_H_DEF__
#define __GROM2ATOMFMO_H_DEF__

/*! \file grom2atomFMO.h 
		\brief Module containing classes that convert between topology and .gro file types into processed data (energies, correlations, etc)
 		
		Big description right hurr.
*/

/*! \brief Thrown on failure to format and interpret data in any of the input files
 */
class InputFileMalformed : public std::logic_error
{ public: 
//! \brief Ctor for class, default constructor that allows for passing description of error	
	explicit InputFileMalformed(const std::string& what_arg) : std::logic_error(what_arg){} };

/*! \brief Thrown on failure to perform io on file
 */
class ReadError : public std::runtime_error
{ public: 
//! \brief Ctor for class, default constructor that allows for passing description of error	
	explicit ReadError(const std::string& what_arg) : std::runtime_error(what_arg){} };


/*! \brief Interprets string data as an alternative format using the stringstream class (e.g. "3.14" into 3.14f), 
 *  \throws std::invalid_argument String cannot be automatically formatted into the desired type
 */
template< typename T > inline T Convert(const std::string& str)
{
	std::istringstream iss(str);
	T obj;
			
	iss >> std::ws >> obj >> std::ws;

	if(!iss.eof()) throw std::invalid_argument("Convert was not given an appropriate string");

	return obj; 
}


struct ConfigData
{
  std::string fmo_traj_path;
  std::string fmo_top_path;
  std::string bcx_itp_path;
  std::string output_path;
  int num_frames;
};

enum Excited{ NO = 0, YES = 1};

/////////////////////////////////////////////////
//  ReadFMOConfig
/////////////////////////////////////////////////
/*! \brief 	Reads the configuration file and loads the topology from files
 * \param 	cfgFileName   INPUT Name of input configuration file
 * \param		trajFileName	OUTPUT Location of .gro input read from config
 * \param		outFileName		OUTPUT Location to save output read from config
 * \param		n_entries			OUTPUT Max number of entries to read
 * \param		atomNames			OUTPUT The name of an atom within a molecule, "MOL,ATM"
 * \param		atomMasses		OUTPUT The mass for corresponding atomName atom
 * \param		atomCharges		OUTPUT The charge for corresponding atomName atom. dQ between grd and exc for chromophores, and Q for all other atoms
 */

int ReadFMOConfig(std::string const& cfgFileName, 
									std::string& trajFileName, 
									std::string& outFileName, 
									int& n_entries, 
									std::vector<std::string>& atomNames, 
									std::vector<float>& atomMasses,
									std::vector<float>& atomCharges);



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
									 Excited xc = Excited::NO);


/////////////////////////////////////////////////
// ReadOneTimeGrom2Atom
/////////////////////////////////////////////////
/*! \brief Loads a single time frame from .gro file
	\param currentPos 		INPUT Position in file to start from. 0 if beginning of file
	\param outPos 				OUTPUT Reference to output variable, outputs the non-atom line before the next time-header
	\param Xi							OUTPUT Position vector to be populated, all data is overwritten
	\param atomTypei			OUTPUT Atom type vector to be populated, all data is overwritten
	\param atomGroupi 		OUTPUT Molecule type vector to be populated, all data is overwritten
	\param systemSize_nm	OUTPUT The output system size for the current frame
	\param fileName 			INPUT  Name of the file to read from
	\param should_print		INPUT  Print the diagnostic print statements from the function

	\throws ReadError
	\throws	InputFileMalformed

	\return void
 */


void ReadOneTimeGrom2Atom(std::streampos const* const currentPos, 
													std::streampos & outPos,
												 	std::vector<float > & Xi,
												 	std::vector<std::string > & atomTypei,
													std::vector<int>& atomGroupi,
													float &systemSize_nm,
												 	std::string const& fileName = "md_traj.gro",
													bool should_print = false);


/////////////////////////////////////////////////
/// AtomDataLookup
/////////////////////////////////////////////////
/* \brief 	Uses lookup tables from the topology file (ParseTopology) to convert .gro data (ReadOneTimeGrom2Atom) into charges and atomtypes
 *
 * \param 	atomTypei 
 * \param		atomMassi_amu : REQ size must be equal to atomTypei size
 * \param		atomSize_nm : REQ size must be equal to atomTypei size
 *
 * \throws  InputFileMalformed
 *		basic exception guarantee
 *
 */
void AtomDataLookup(std::vector<std::string > const& atomTypei,
												 std::vector<float > & atomMassi_amu,
												 std::vector<float > & atomSizei_nm , 
												 std::vector<float > & atomChargei_nm , 
												 std::vector<int > const & atomGroupTablei,
												 int excitationSite,
												 std::vector<std::string > const& atomTypeTable, 
												 std::vector<float > const & atomMassTable,
												 std::vector<float > const & atomChargeTable);

void AtomDataLookup_v2(std::vector<std::string > const& atomTypei,
												 std::vector<float > & atomMassi_amu,
												 std::vector<float > & atomSizei_nm , 
												 std::vector<float > & atomChargei_nm , 
												 std::vector<int > const & atomGroupTablei,
												 int excitationSite,
												 std::vector<float > & atomGroundCharges_excitedGroup,
												 std::vector<std::string > const& atomTypeTable, 
												 std::vector<float > const & atomMassTable,
												 std::vector<float > const & atomChargeTable);


/////////////////////////////////////////////////
// ComputeCDC_v1
/////////////////////////////////////////////////
/*! \brief 	Uses input from topology and .gro (AtomDataLookup) to compute the energy vector (partitioned by molecule) with a particular site excited
 * \param		chromoSite	INPUT Site to compute dE for
 * \param		Xi					INPUT Position of atom i [nm]
 * \param		qi					INPUT Charge of atom i (dQ for chromo, Q for others) [e]
 * \param		atomGroups	INPUT 
 * \param 		cdc_kcal		OUTPUT of cdc computation
 * \param 	excitedAtom_groundCharges		INPUT/OPTIONAL The ground charges of the excited molecule, in the order that they appear in qi. If excluded, only intermolecular charges will be computed. Asserted to be equal in size to the number of entries of atomGroups with group chromoSite, and asserts that no sites are at the same location except the corresponding indices, which provides a fairly through check.
 *
 * Computes the CDC estimate for site energy, using the atomic charges. The results are clustered by atomGroup to form a coarse grained picture of contributions to the energetics.
 * The indexing for the outputs of cdc_kcal are by atomGroup number.
 * The algorithm contained is brutish, but sufficient:
 * First, the indexes for the current chromophore are collected into a vector and the maximum group number $nGroups is computed 
 * Then a loop is constructed to go between 1 and $nGroups.
 * In each loop, all of the atoms are iterated through. A match is found when the abs() of the atomGroup for the atom is equal to the loop index but not equal to the chromophore index. For each match, the contribution between each chromophore atom and the current atom is added to cdc_kcal[i].
 * The bigO is thus the larger of (nGroups * nAtoms) [for index looping] or (nAtoms * nAtomInChromophore) [for computation].
 *
 * Values in atomGroups go from 1 - nGroups for protein and (-1) - (-7) for chromophores.
 * Negative values of atomGroups are treated as the first indexes, in the order of their negative index values (e.g. -3, -2, -1, 1, 2, 3).
 * Size of cdc_kcal is nGroups + nChromo. 
 * Index of group g (one-based indexing) is [nChromo + g - 1].
 * Index of chromophore i is [nChromo + i] when i is a one-based indexing negative number.
 * 
 * 
 * Asserts
 *
 * \warning The order of the ground state charges corresponding to the excited molecule's atoms must be in the order that they exist in the atom list in atomGroups. Shuffling things will break everything.
 */

void ComputeCDC_v1( int chromoSite, 
										std::vector<float > & cdc_kcal,
										std::vector<float >const& Xi,
										std::vector<float >const& qi,
										std::vector<int > const & atomGroups,
										std::vector<float > const& excitedAtom_groundCharges = std::vector<float>(0,0));


/////////////////////////////////////////////////
// ComputeErfDensity_v1 
/////////////////////////////////////////////////
/*! \brief first implementation of a coarsegraining procedure 
 * \param Xi_nm            INPUT positions of atom i 
 * \param atomMassi_amu    INPUT mass of atom i
 * \param atomSizei_nm     INPUT size of atom i (to be used as a width for the gaussian spread)
 * \param systemSize_nm    INPUT size of the system in one dimension, assuming cubic.
 * \param Nb               INPUT Number of boxes to use along one length
 * \param density_amu      OUTPUT density field out
 * */

void ComputeErfDensity_v1(std::vector<float > const & Xi_nm,
													std::vector<float > const & atomMassi_amu,
													std::vector<float > const & atomSizei_nm,
													float systemSize_nm,
													int Nb,
													std::vector<float > & density_amu);

/*! \brief Computes Erf-coarsegrained momentum correlation field
 *  \param 	Xi_nm						The vector of atom positions, in nm
 *  \param  Vi_nm_ps			  The vector of atom positions, in nm per ps
 * 	\param 	atomMassi_amu 	The vector of atom masses, in amu
 * 	\param 	atomSizei_nm		The vector of atom sizes, in nm
 * 	\param 	systemSize_nm		The length of the system (assumes cubic)
 * 	\param 	Nb							The number of boxes to use for a single length
 * 	\param 	momentum_amunm_ps			The output density field for this configuration, all data is overwritten
 */
void ComputeErfMomentum_v1(std::vector<float > const & Xi_nm,
													 std::vector<float > const & Vi_nm_ps,
													 std::vector<float > const & atomMassi_amu,
													 std::vector<float > const & atomSizei_nm,
													 float systemSize_nm,
													 int Nb,
													 std::vector<float > & momentum_amunm_ps);

													

























/*
struct AtomParams{
	std::vector<float > Xi_nm;
	std::vector<float > Vi_nm_ps;
	std::vector<std::string > atomTypei;
	float time_ps;
	float boxLength_nm;
	int atomTotal;
};

class AtomLookup{
	public:
		AtomLookup(): m_atomNames(), m_atomMass_amu(), m_atomSize_nm(){}
		void AddAtom(std::string const& atomName, float atomMass_amu, float atomSize_nm);
		 *! GetData 	inputs the mass and the size into the last two arguments by dictionary-style lookup.
		 *	
		 *	\throw 		InputFileMalformed
		 *		strong throw guarantee
		 *
		void GetData(std::string const& atomName, float & massOut_amu, float & sizeOut_nm);
		std::vector<float > ConvertNameToMass(std::vector<std::string > & atomNames);
		std::vector<float > ConvertNameToSize(std::vector<std::string > & atomNames);

	private:
		std::vector<std::string > m_atomNames;
		std::vector<float > m_atomMass_amu;
		std::vector<float > m_atomSize_nm;

}
*/

#endif //__GROM2ATOM_H_DEF__
