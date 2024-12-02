package main

import (
	"fmt"
)

func main() {
	// filepath := os.Args[1]
	// time, _ := strconv.ParseFloat(os.Args[2], 64)
	// protein, err := readProteinFromFile(filepath) //"../data/calmodulin_noCA.pdb"
	// Check(err)

	filepath := "../data/calmodulin_noCA.pdb"
	time := 0.00001
	protein, err := readProteinFromFile(filepath)
	Check(err)

	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/OPLS_atom_charge.rtp")
	Check(err)

	// Assign charges to the protein's atoms
	(&protein).AssignChargesToProtein(chargeData)
	fmt.Println(protein.Residue[0])
	fmt.Println(protein.Residue[0].Atoms[0])
	fmt.Println(protein.Residue[0].Atoms[1])
	fmt.Println(protein.Residue[0].Atoms[2])
	// Check if the assigned charges are correct
	//CheckAssignedCharges(&protein, chargeData)

	residueParameterValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
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
	timepoints := SimulateMD(*initialProtein, time, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)
	RSMD := CalculateRMSD(timepoints)
	TemporaryPlot(RSMD, time)
	writeRMSD(RSMD)
	WriteProteinToPDB(&timepoints[len(timepoints)-1], "result/output.pdb")

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

func printProtein(protein *Protein) {
	fmt.Printf("Protein Name: %s\n", protein.Name)
	for _, residue := range protein.Residue {
		fmt.Printf("  Residue Name: %s, ID: %d, ChainID: %s\n", residue.Name, residue.ID, residue.ChainID)
		for _, atom := range residue.Atoms {
			fmt.Printf("    Atom Index: %d, Element: %s, Position: (%.2f, %.2f, %.2f)\n",
				atom.index, atom.element, atom.position.x, atom.position.y, atom.position.z)
		}
	}
}

func CheckAssignedCharges(protein *Protein, chargeData map[string]map[string]float64) {
	for _, residue := range protein.Residue {
		residueName := residue.Name

		// Get the charge data for this residue, if it exists
		residueChargeData, residueExists := chargeData[residueName]

		if !residueExists {
			fmt.Printf("Warning: No charge data found for residue %s\n", residueName)
			continue
		}

		for _, atom := range residue.Atoms {
			atomName := atom.element

			// Try to get the charge data for this atom
			expectedCharge, atomExists := residueChargeData[atomName]
			if !atomExists {
				fmt.Printf("Warning: No charge data found for atom %s in residue %s\n", atomName, residueName)
				continue
			}

			// Check if the assigned charge matches the expected charge
			if atom.charge != expectedCharge {
				fmt.Printf("Discrepancy found: Residue %s, Atom %s, Assigned Charge: %f, Expected Charge: %f\n",
					residueName, atomName, atom.charge, expectedCharge)
			}
		}
	}
}
