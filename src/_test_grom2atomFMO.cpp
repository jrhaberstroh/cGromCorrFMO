#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "grom2atomFMO.h"


TEST_CASE("","[ParseTopology]"){
  std::string filename = "./test/4BCL_pp.top"; 
  std::vector<std::string > mergeNameTable;
  std::vector<float > massTable;
  std::vector<float > chargeTable;

  SECTION("Try parsing topology"){
    ParseTopology(filename, mergeNameTable, massTable, chargeTable);
  }

}


TEST_CASE("Looking up atom data with constructed arrays","[AtomDataLookup_v2]"){
  std::vector<std::string> atomName;
  std::vector<float > atomMass;
  std::vector<float > atomSize;
  std::vector<float > atomCharge;
  std::vector<int > atomGroup;
  std::vector<float > excitedAtomGroundCharges;

  atomName.push_back("BCL,MG");
  atomGroup.push_back(-1); // Chromophores are given negative indices

  std::vector<std::string > atomTypeTable;
  std::vector<float > atomMassTable;
  std::vector<float > atomChargeTable;
 
  // Ground state data 
  atomTypeTable.push_back("BCL,MG");
  atomMassTable.push_back(33.0);
  atomChargeTable.push_back(.20);

  // Excited state data
  atomTypeTable.push_back("BCL,MG,EXC");
  atomMassTable.push_back(0.2);
  atomChargeTable.push_back(.70);

  SECTION("Ground-state lookup"){
    AtomDataLookup_v2(atomName, atomMass, atomSize, atomCharge, atomGroup, 0, excitedAtomGroundCharges, atomTypeTable, atomMassTable, atomChargeTable);

    REQUIRE(atomMass[0] == Approx(33.0));
    REQUIRE(atomCharge[0] == Approx(.20));
  }

  SECTION("Excited-state lookup"){
    AtomDataLookup_v2(atomName, atomMass, atomSize, atomCharge, atomGroup, 1, excitedAtomGroundCharges, atomTypeTable, atomMassTable, atomChargeTable);

    REQUIRE(atomCharge[0] == Approx(.70));
    REQUIRE(atomMass[0] == Approx(.2));
    REQUIRE(excitedAtomGroundCharges[0] == Approx(.20));
  }
}
