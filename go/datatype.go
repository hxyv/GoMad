package main

type TriTuple struct {
	x float64
	y float64
	z float64
}

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
	position    TriTuple
	velocity    TriTuple
	force       TriTuple
	accelerated TriTuple
	mass        float64
	element     string
}

//ATOM      1       N        MET          A       1              26.457  24.555  27.324  1.00 20.00           N
//Atom   number atomName ResidueName  chainName Residueposition    X       Y       Z                       element number
