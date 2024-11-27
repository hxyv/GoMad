package main

import "fmt"

func main() {
	//secline, err := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	//Check(err)
	//fmt.Println(secline.atomPair)
	//printParameterDatabase(&secline)

	protein, err := readProteinFromFile("../data/calmodulin_noCA.pdb")
	fmt.Println(len(protein.Residue))
	Check(err)

	// Parse the charge data file
	chargeData, err := parseChargeFile("../data/gromacs43_atom_charge.rtp")
	if err != nil {
		fmt.Printf("Error parsing charge file: %v\n", err)
		return
	}

	// Assign charges to the protein's atoms
	(&protein).AssignChargesToProtein(chargeData)

	// // Print out the charges to verify
	// for _, residue := range protein.Residue {
	// 	fmt.Printf("Residue: %s\n", residue.Name)
	// 	for _, atom := range residue.Atoms {
	// 		fmt.Printf("  Atom: %s, Charge: %f\n", atom.element, atom.charge)
	// 	}
	// }

	verletList := NewVerletList()
	verletList.BuildVerlet(&protein)
	// Print nearest neighbors in the Verlet list
	// for atom, neighbors := range verletList.Neighbors {
	// 	fmt.Printf("Atom %d (%s) neighbors:\n", atom.index, atom.element)
	// 	for _, neighbor := range neighbors {
	// 		fmt.Printf("  Neighbor Atom %d (%s) at (%.2f, %.2f, %.2f)\n", neighbor.index, neighbor.element, neighbor.position.x, neighbor.position.y, neighbor.position.z)
	// 	}
	// }
	// printProtein(&protein)

	// Calculate energy for each atom
	for i := 0; i < len(protein.Residue); i++ {
		energyMap := CalculateTotalUnbondEnergy(protein.Residue[i].Atoms, verletList)

		// Output the energy for each atom
		for index, energy := range energyMap {
			fmt.Printf("Atom %d: Energy = %f\n", index, energy)
		}
	}

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
