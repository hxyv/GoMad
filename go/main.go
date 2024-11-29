package main

import "fmt"

func main() {
	protein, err := readProteinFromFile("../data/calmodulin_noCA.pdb")
	Check(err)

	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/gromacs43_atom_charge.rtp")
	Check(err)

	// Assign charges to the protein's atoms
	(&protein).AssignChargesToProtein(chargeData)

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
	time := 0.00000001
	timepoints := SimulateMD(*initialProtein, time, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)
	RSMD := CalculateRMSD(timepoints)
	TemporaryPlot(RSMD, time)
	WriteProteinToPDB(&timepoints[len(timepoints)-1], "output.pdb")

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
