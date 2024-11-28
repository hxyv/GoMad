package main

import (
	"math"
)

// a1: the atom 1
// a2: the atom 2
// epsilon: the vacuum dielectric permittivity
// r: distance between q1 and q2
func CalculateElectricPotentialEnergy(a1, a2 *Atom, r float64) float64 {
	return (a1.charge * a2.charge) / (4 * math.Pi * epsilon * r * 1e-10)
}

// A: coefficient 1
// B: coefficient 2
// r: distance between atom 1 and atom 2
func CalculateLJPotentialEnergy(B, A, r float64) float64 {
	r_6 := math.Pow(r, 6)
	r_12 := r_6 * r_6

	return (A / r_6) - (B / r_12)
}

// func CalculateLJForce take the result of CalculateLJPotentialEnergy
// return the LJForce between two atoms
func CalculateLJForce(A, B, r float64) float64 {
	r_6 := math.Pow(r, 6)
	r_12 := r_6 * r_6

	return (A / r_12) - (B / r_6)
}

// C: coefficient 1
// D: coefficient 2
// r: distance between atom 1 and atom 2
func CalculateHydrogenBondEnergy(C, D, r float64) float64 {
	r_10 := math.Pow(r, 10)
	r_12 := math.Pow(r, 12)

	return (C / r_10) - (D / r_12)
}

func NewVerletList() *VerletList {
	return &VerletList{
		Neighbors: make(map[*Atom][]*Atom),
		Cutoff:    verletCutOff,
		Buffer:    verletBuffer,
	}
}

func (v *VerletList) BuildVerlet(protein *Protein) {
	cutoffPlusBuffer := v.Cutoff + v.Buffer
	v.Neighbors = make(map[*Atom][]*Atom)

	for _, residue := range protein.Residue {
		for _, atom := range residue.Atoms {
			for _, targetResidue := range protein.Residue {
				for _, targetAtom := range targetResidue.Atoms {
					if atom != targetAtom {
						distance := Distance(atom.position, targetAtom.position)
						if distance <= cutoffPlusBuffer && (targetAtom.index < atom.index+3 || targetAtom.index > atom.index+3) {
							v.Neighbors[atom] = append(v.Neighbors[atom], targetAtom)
						}
					}
				}
			}
		}
	}
}

func (protein *Protein) AssignChargesToProtein(chargeData map[string]map[string]AtomChargeData) {
	for _, residue := range protein.Residue {
		residueName := residue.Name

		// Get the charge data for this residue, if it exists
		residueChargeData, residueExists := chargeData[residueName]

		for _, atom := range residue.Atoms {
			atomName := atom.element

			if residueExists {
				// Try to get the charge data for this atom
				atomChargeData, atomExists := residueChargeData[atomName]
				if atomExists {
					// Assign the charge from the charge data
					atom.charge = atomChargeData.AtomCharge
					continue
				}
			}

			// If there's no data for this atom or residue, assign a charge of 0
			atom.charge = 0.0
		}
	}
}

func CalculateTotalUnbondEnergy(atoms []*Atom, verletList *VerletList) map[int]float64 {
	energyMap := make(map[int]float64)

	for _, atom1 := range atoms {
		// Initialize energy for atom1
		energyMap[atom1.index] = 0.0

		// Access the Neighbors map using the dereferenced verletList
		neighbors, exists := verletList.Neighbors[atom1]
		if !exists {
			continue
		}

		for _, atom2 := range neighbors {
			// Compute the distance between atom1 and atom2
			r := Distance(atom1.position, atom2.position)

			// Calculate the electric potential energy between atom1 and atom2
			electricPotentialEnergy := CalculateElectricPotentialEnergy(atom1, atom2, r)

			// Update the energy map for atom1
			energyMap[atom1.index] += electricPotentialEnergy

			nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
			Check(error)

			parameterList := SearchParameter(2, nonbondedParameter, atom1, atom2)
			if len(parameterList) != 1 {
				// Calculate the electric potential energy between atom1 and atom2
				LJPotentialEnergy := CalculateLJPotentialEnergy(parameterList[0], parameterList[1], r)

				// Update the energy map for atom1
				energyMap[atom1.index] += LJPotentialEnergy

			}
		}
	}

	return energyMap
}
