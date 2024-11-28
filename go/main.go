package main

import (
	"bufio"
	"fmt"
	"os"
)

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

	// Print out the charges to verify
	/*
		for _, residue := range protein.Residue {
			fmt.Printf("Residue: %s\n", residue.Name)
			for _, atom := range residue.Atoms {
				fmt.Printf("  Atom: %s, Charge: %f\n", atom.element, atom.charge)
			}
		}
	*/
	fmt.Println((&protein).Residue[0].ChainID)
	WriteProteinToPDB(&protein, "result/output.pdb")
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

// function WriteProteinToPDB takes a protein structure as input
// return a pdb file
func WriteProteinToPDB(protein *Protein, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// Write the title section
	if _, err := writer.WriteString(fmt.Sprintf("HEADER    %s\n", protein.Name)); err != nil {
		return err
	}

	atomIndex := 1
	for _, residue := range protein.Residue {
		for _, atom := range residue.Atoms {
			// Format atom data according to the PDB file format
			_, err := writer.WriteString(fmt.Sprintf(
				"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %-2s\n",
				atomIndex,                                         // Atom serial number
				atom.element,                                      // Atom name
				residue.Name,                                      // Residue name
				residue.ChainID,                                   // Chain identifier
				residue.ID,                                        // Residue sequence number
				atom.position.x, atom.position.y, atom.position.z, // Atom coordinates
				atom.element, // Element symbol
			))
			if err != nil {
				return err
			}
			atomIndex++
		}
	}

	// Write the termination line
	if _, err := writer.WriteString("END\n"); err != nil {
		return err
	}

	return writer.Flush()
}
