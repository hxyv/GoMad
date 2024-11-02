package main

import (
	"fmt"
)

func main() {
	secline, err := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(err)
	//fmt.Println(secline.atomPair)
	printParameterDatabase(&secline)

}

func Check(err error) {
	if err != nil {
		panic(err)
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
