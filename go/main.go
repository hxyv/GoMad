// The main file which performs the model
package main

func main() {
	// filepath := os.Args[1]
	// time, _ := strconv.ParseFloat(os.Args[2], 64)
	// protein, err := readProteinFromFile(filepath) //"../data/calmodulin_noCA.pdb"
	// Check(err)

	filepath := "../data/hsucp1_af3.pdb"
	time := 1.0
	protein, err := readProteinFromFile(filepath)
	Check(err)

	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/OPLS_atom_charge.rtp")
	Check(err)

	// Assign charges to the protein's atoms
	(&protein).AssignChargesToProtein(chargeData)

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

	iteration1 := 50
	iteration2 := 100
	initialProtein := PerformEnergyMinimization(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter, iteration1)

	timepoints := SimulateMD(*initialProtein, time, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter, iteration2)
	RMSD := CalculateRMSD(timepoints)
	TemporaryPlot(RMSD, time)
	writeRMSD(RMSD)
	WriteProteinToPDB(&timepoints[len(timepoints)-1], "result/output.pdb")

}

func Check(err error) {
	if err != nil {
		panic(err)
	}

}
