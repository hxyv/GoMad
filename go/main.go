package main

import (
	"fmt"
)

func main() {
	secline, _ := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	fmt.Println(secline.atomPair)

}
