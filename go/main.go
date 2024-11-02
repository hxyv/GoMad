package main

import "fmt"

func main() {
	protein, _ := readProteinFromFile("../data/calmodulin_noCA.pdb")
	for _, residue := range protein.Residue {
		fmt.Printf("Residue: %s, ID: %d, Chain: %s\n", residue.Name, residue.ID, residue.ChainID)
		for _, atom := range residue.Atoms {
			fmt.Printf("  Atom: index: %d,Element=%s, Mass=%.2f, Position=(%.2f, %.2f, %.2f)\n",
				atom.index, atom.element, atom.mass, atom.position.x, atom.position.y, atom.position.y)
		}
	}

}

/*
for _, residue := range protein.Residue {
		fmt.Printf("Residue: %s, ID: %d, Chain: %s\n", residue.Name, residue.ID, residue.ChainID)
		for _, atom := range residue.Atoms {
			fmt.Printf("  Atom: Element=%s, Mass=%.2f, Position=(%.2f, %.2f, %.2f)\n",
				atom.element, atom.mass, atom.position.x, atom.position.y, atom.position.y)
		}
	}

// A function to print the contents of parameterDatabase
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
*/
