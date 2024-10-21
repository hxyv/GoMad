package main

import "fmt"

func main() {
	protien, _ := readProteinFromFile("calmodulin_noCA.pdb")

	fmt.Println(protien.Residue[0].Atoms[0])
}
