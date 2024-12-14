// The main file which performs the model
package main

import (
	"os"
	"strconv"
)

func main() {
	// Parse input values
	filepath := os.Args[1]
	iteration1, _ := strconv.Atoi(os.Args[2])     // steps for energy minimization
	iteration2, _ := strconv.Atoi(os.Args[3])     // steps for simulation
	protein, err := readProteinFromFile(filepath) //"../data/calmodulin_noCA.pdb"
	Check(err)

	// Commented for testing
	// filepath := "../data/calmodulin_noCA.pdb"
	// iteration1 := 2 // steps for energy minimization
	// iteration2 := 2 // steps for simulation
	// protein, err := readProteinFromFile(filepath)
	// Check(err)

	// Each step is 1 fs
	time := 1.0
	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/OPLS_atom_charge.rtp")
	Check(err)

	// Assign charges to the protein's atoms
	(&protein).AssignChargesToProtein(chargeData)

	// Parse related parameter files
	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	// First perform energy minimization
	initialProtein := PerformEnergyMinimization(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter, iteration1)
	// Then simulate the protein
	timepoints := SimulateMD(*initialProtein, time, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter, iteration2)
	RMSD := CalculateRMSD(timepoints)
	TemporaryPlot(RMSD, time)                                              // Plot RMSD in GO
	writeRMSD(RMSD)                                                        // Write RMSD to csv
	WriteProteinToPDB(&timepoints[len(timepoints)-1], "result/output.pdb") // Write the final protein to pdb

}

func Check(err error) {
	if err != nil {
		panic(err)
	}

}
