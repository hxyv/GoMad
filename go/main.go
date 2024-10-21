package main

import "fmt"

func main() {
	protien, _ := readProteinFromFile("../data/calmodulin_noCA.pdb")

	fmt.Println(protien.Residue[0])
	fmt.Println(protien.Residue[0].Atoms[0])
}
