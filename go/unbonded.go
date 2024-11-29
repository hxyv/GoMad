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

func CalculateTotalUnbondEnergy(p *Protein, nonbondedParameter parameterDatabase) map[int]float64 {
	energyMap := make(map[int]float64)
	verletList := NewVerletList()
	verletList.BuildVerlet(p)

	for _, residue := range p.Residue {
		for _, atom1 := range residue.Atoms {
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

				// Calculate the Lennard-Jones potential energy between atom1 and atom2
				parameterList := SearchParameter(2, nonbondedParameter, atom1, atom2)
				if len(parameterList) == 2 {
					LJPotentialEnergy := CalculateLJPotentialEnergy(parameterList[0], parameterList[1], r)

					// Update the energy map for atom1
					energyMap[atom1.index] += LJPotentialEnergy
				}
			}
		}
	}

	return energyMap
}

func CalculateTotalUnbondedEnergyForce(p *Protein, nonbondedParameter parameterDatabase) (float64, map[int]*TriTuple) {
	forceMap := make(map[int]*TriTuple)
	totalEnergy := 0.0
	verletList := NewVerletList()
	verletList.BuildVerlet(p)

	for _, residue := range p.Residue {
		for _, atom1 := range residue.Atoms {
			// Initialize force for atom1
			forceMap[atom1.index] = &TriTuple{0.0, 0.0, 0.0}

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
				totalEnergy += electricPotentialEnergy

				// Calculate the electric force between atom1 and atom2
				electricForce := CalculateElectricForce(atom1, atom2, r)

				// Update the force map for atom1
				forceMap[atom1.index].x += electricForce.x
				forceMap[atom1.index].y += electricForce.y
				forceMap[atom1.index].z += electricForce.z

				// Calculate the Lennard-Jones potential energy between atom1 and atom2
				parameterList := SearchParameter(2, nonbondedParameter, atom1, atom2)
				if len(parameterList) == 2 {
					LJPotentialEnergy := CalculateLJPotentialEnergy(parameterList[0], parameterList[1], r)
					totalEnergy += LJPotentialEnergy

					// Calculate the Lennard-Jones force between atom1 and atom2
					LJForce := CalculateLJForce(atom1, atom2, parameterList[0], parameterList[1], r)

					// Update the force map for atom1
					forceMap[atom1.index].x += LJForce.x
					forceMap[atom1.index].y += LJForce.y
					forceMap[atom1.index].z += LJForce.z
				}
			}
		}
	}

	return totalEnergy, forceMap
}

func CalculateElectricForce(a1, a2 *Atom, r float64) TriTuple {
	forceMagnitude := (a1.charge * a2.charge) / (4 * math.Pi * epsilon * r * r * 1e-10)
	unitVector := TriTuple{
		x: (a2.position.x - a1.position.x) / r,
		y: (a2.position.y - a1.position.y) / r,
		z: (a2.position.z - a1.position.z) / r,
	}
	return TriTuple{
		x: forceMagnitude * unitVector.x,
		y: forceMagnitude * unitVector.y,
		z: forceMagnitude * unitVector.z,
	}
}

func CalculateLJForce(a1, a2 *Atom, A, B, r float64) TriTuple {
	r_6 := math.Pow(r, 6)
	r_12 := r_6 * r_6

	forceMagnitude := (12*B/r_12 - 6*A/r_6)
	unitVector := TriTuple{
		x: (a2.position.x - a1.position.x) / r,
		y: (a2.position.y - a1.position.y) / r,
		z: (a2.position.z - a1.position.z) / r,
	}
	return TriTuple{
		x: forceMagnitude * unitVector.x,
		y: forceMagnitude * unitVector.y,
		z: forceMagnitude * unitVector.z,
	}
}
