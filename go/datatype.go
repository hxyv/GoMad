package main

type Protein struct {
	Name    string
	Residue []*Residue
}

type Residue struct {
	Name    string
	ID      int
	ChainID string
	Atoms   []*Atom
}

type Atom struct {
	Position    triArray
	Velocity    triArray
	Force       triArray
	Accelerated triArray
	Mass        float64
	Element     string
}

type triArray struct {
	X, Y, Z float64
}

//ATOM      1       N        MET          A       1              26.457  24.555  27.324  1.00 20.00           N
//Atom   number atomName ResidueName  chainName Residueposition    X       Y       Z                       element number
