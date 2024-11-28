package main

import "fmt"

func main() {
	/*
		secline, err := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
		Check(err)
		//fmt.Println(secline.atomPair)
		printParameterDatabase(&secline)
	*/
	protein, err := readProteinFromFile("../data/calmodulin_noCA.pdb")

	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/gromacs43_atom_charge.rtp")
	if err != nil {
		fmt.Printf("Error parsing charge file: %v\n", err)
		return
	}

	// Assign charges to the protein's atoms
	(&protein).assignChargesToProtein(chargeData)

	/*
		// Print out the charges to verify
		for _, residue := range protein.Residue {
			fmt.Printf("Residue: %s\n", residue.Name)
			for _, atom := range residue.Atoms {
				fmt.Printf("  Atom: %s, Charge: %f\n", atom.element, atom.charge)
			}
		}
	*/
	residueParameterValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
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
	initialProtein := PerformEnergyMinimization(&protein, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)
	time := 10.0
	timepoints := SimulateMD(*initialProtein, time, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)
	RSMD := CalculateRMSD(timepoints)
	TemporaryPlot(RSMD)
}

func Check(err error) {
	if err != nil {
		panic(err)
	}

}

func printParameterDatabase(db *parameterDatabase) {
	for _, pair := range db.atomPair {
		fmt.Println("Atom Names:")
		for _, name := range pair.atomName {
			fmt.Printf("  %s\n", name)
		}
		fmt.Printf("Function: %d\n", pair.Function)
		fmt.Println("Parameters:")
		for _, param := range pair.parameter {
			fmt.Printf("  %.2f\n", param)
		}
		fmt.Println()
	}
}
