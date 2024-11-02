package main

import (
	"fmt"
)

func main() {
	secline, _ := ReadAminoAcidsPara("../data/test.itp")
	fmt.Println(secline["ALA"].dihedrals[0])

}
