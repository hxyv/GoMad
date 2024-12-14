package main

import (
	"math"
)

// CalculateLJPotentialEnergy calculates the Lennard-Jones potential energy
// between two atoms given coefficients A and B, and the distance r between them.
// Input: two Atom object and float64 r (distance)
// Output: a float64 of electrostatic potential energy
func CalculateElectricPotentialEnergy(a1, a2 *Atom, r float64) float64 {
	return f * a1.charge * a2.charge / (r * 0.1)
}

// CalculateLJPotentialEnergy calculates the Lennard-Jones potential energy
// between two atoms given coefficients A and B, and the distance r between them.
// Input: A, B, r float64
// Output: a float64 of LJ potential energy
func CalculateLJPotentialEnergy(B, A, r float64) float64 {
	r = r * 0.1
	r_6 := math.Pow(r, 6)
	r_12 := r_6 * r_6
	LJ := (A / r_12) - (B / r_6)
	return LJ
}

// NewNeighborList creates a new Verlet list for neighbor atoms.
// Input: None
// Output: A pointer to a newly created VerletList structure containing:
// - An empty map of atoms to their neighboring atoms.
// - A predefined cutoff distance for neighbor detection.
func NewNeighborList() *VerletList {
	return &VerletList{
		Neighbors: make(map[*Atom][]*Atom),
		Cutoff:    verletCutOff,
	}
}

// BuildNeighbor constructs the neighbor list for a given protein.
// Input: protein: A pointer to a Protein structure containing residues and atoms.
// Output: Updates the VerletList's Neighbors field.
func (v *VerletList) BuildNeighbor(protein *Protein) {
	v.Neighbors = make(map[*Atom][]*Atom)

	// Iterate through all residues and atoms in the protein.
	for _, residue := range protein.Residue {
		for _, atom := range residue.Atoms {
			// Initialize an empty neighbor list for the current atom.
			v.Neighbors[atom] = []*Atom{}

			// Check all other atoms in the protein for potential neighbors.
			for _, otherResidue := range protein.Residue {
				for _, otherAtom := range otherResidue.Atoms {
					// Skip if the atoms are the same.
					if atom == otherAtom {
						continue
					}

					// Skip if the other atom is within 3 bonds of the current atom.
					if otherAtom.index >= atom.index-4 && otherAtom.index <= atom.index+4 {
						continue
					}

					// Calculate the distance between the current atom and the other atom.
					distance := Distance(atom.position, otherAtom.position)

					// Add the other atom to the neighbor list if it satisfies the distance criteria.
					if distance <= v.Cutoff && distance > 2.5 {
						v.Neighbors[atom] = append(v.Neighbors[atom], otherAtom)
					}
				}
			}
		}
	}
}

func (protein *Protein) AssignChargesToProtein(chargeData map[string]map[string]float64) {
	for _, residue := range protein.Residue {
		residueName := residue.Name

		// Get the charge data for this residue, if it exists
		residueChargeData, residueExists := chargeData[residueName]

		for _, atom := range residue.Atoms {
			atomName := atom.element

			if residueExists {
				// Try to get the charge data for this atom
				atomCharge, atomExists := residueChargeData[atomName]
				if atomExists {
					// Assign the charge from the charge data
					atom.charge = atomCharge
					continue
				}
			}

			// If there's no data for this atom or residue, assign a charge of 0
			atom.charge = 0.0
		}
	}
}

// CalculateTotalUnbondedEnergyForce calculates the total unbonded energy and forces in a protein.
// Input:
// - p: A pointer to a Protein structure containing residues and atoms.
// - nonbondedParameter: A parameter database for nonbonded interactions.
// Output:
// - A float64 representing the total unbonded energy of the protein.
// - A map[int]*TriTuple containing the forces on each atom indexed by their unique IDs.
func CalculateTotalUnbondedEnergyForce(p *Protein, nonbondedParameter parameterDatabase) (float64, map[int]*TriTuple) {
	// Initialize the force map to store forces for each atom and the total energy to 0.
	forceMap := make(map[int]*TriTuple)
	totalEnergy := 0.0

	// Create and build a Verlet neighbor list for the protein.
	verletList := NewNeighborList()
	verletList.BuildNeighbor(p)

	// Iterate through each residue in the protein.
	for _, residue := range p.Residue {
		for _, atom1 := range residue.Atoms {

			// Initialize the force vector for the current atom.
			forceMap[atom1.index] = &TriTuple{0.0, 0.0, 0.0}

			// Retrieve the neighbors of the current atom using the Verlet list.
			neighbors, exists := verletList.Neighbors[atom1]
			if !exists {
				continue
			}

			// Iterate through the neighbors of the current atom.
			for _, atom2 := range neighbors {
				// Compute the distance between atom1 and atom2.
				r := Distance(atom1.position, atom2.position)

				// Retrieve Lennard-Jones parameters for the atom pair.
				parameterList := SearchParameter(2, nonbondedParameter, atom1, atom2)
				if len(parameterList) == 2 {
					// Calculate the Lennard-Jones potential energy and add it to the total energy.
					LJPotentialEnergy := CalculateLJPotentialEnergy(parameterList[0], parameterList[1], r)
					totalEnergy += LJPotentialEnergy

					// Calculate the Lennard-Jones force and update the force map for atom1.
					LJForce := CalculateLJForce(atom1, atom2, parameterList[0], parameterList[1], r)
					forceMap[atom1.index].x += LJForce.x
					forceMap[atom1.index].y += LJForce.y
					forceMap[atom1.index].z += LJForce.z
				}

				// Skip electric interactions if either atom has no charge.
				if atom1.charge == 0.0 || atom2.charge == 0.0 {
					continue
				}

				// Calculate the electric potential energy and add it to the total energy.
				electricPotentialEnergy := CalculateElectricPotentialEnergy(atom1, atom2, r)
				totalEnergy += electricPotentialEnergy

				// Calculate the electric force and update the force map for atom1.
				electricForce := CalculateElectricForce(atom1, atom2, r)
				forceMap[atom1.index].x += electricForce.x
				forceMap[atom1.index].y += electricForce.y
				forceMap[atom1.index].z += electricForce.z
			}
		}
	}

	// Return the total unbonded energy and the force map.
	return totalEnergy, forceMap
}

// CalculateElectricForce computes the electric force exerted between two atoms.
// Input:
// - a1, a2: Pointers to Atom structures representing the two interacting atoms.
// - r: The distance between the two atoms.
// Output:
// - A TriTuple representing the electric force vector acting on a1 due to a2.
func CalculateElectricForce(a1, a2 *Atom, r float64) TriTuple {
	// Compute the unit vector from a1 to a2.
	unitVector := TriTuple{
		x: (a2.position.x - a1.position.x) / r,
		y: (a2.position.y - a1.position.y) / r,
		z: (a2.position.z - a1.position.z) / r,
	}

	// Convert the distance from Ångströms to nanometers.
	r = r * 0.1

	// Calculate the magnitude of the charges.
	chargeMagnitude := a1.charge * a2.charge

	// Initialize the force magnitude.
	forceMagnitude := 0.0

	// Compute the force magnitude using Coulomb's law, applying the appropriate sign.
	if chargeMagnitude > 0.0 {
		forceMagnitude = f * chargeMagnitude / (r * r)
	} else {
		forceMagnitude = -f * chargeMagnitude / (r * r)
	}

	// Scale the force vector and apply unit vector components.
	force := TriTuple{
		x: forceMagnitude * unitVector.x * 1e-5,
		y: forceMagnitude * unitVector.y * 1e-5,
		z: forceMagnitude * unitVector.z * 1e-5,
	}

	// Return the electric force vector.
	return force
}

// CalculateLJForce computes the Lennard-Jones force exerted between two atoms.
// Input:
// - a1, a2: Pointers to Atom structures representing the two interacting atoms.
// - B, A: Lennard-Jones parameters (constants) for the atom pair.
// - r: The distance between the two atoms.
// Output:
// - A TriTuple representing the Lennard-Jones force vector acting on a1 due to a2.
func CalculateLJForce(a1, a2 *Atom, B, A, r float64) TriTuple {
	// Compute the unit vector from a1 to a2.
	unitVector := TriTuple{
		x: (a2.position.x - a1.position.x) / r,
		y: (a2.position.y - a1.position.y) / r,
		z: (a2.position.z - a1.position.z) / r,
	}

	// Convert the distance from Ångströms to nanometers.
	r = r * 0.1

	// Compute the sixth and twelfth powers of the distance.
	r_6 := math.Pow(r, 6)
	r_12 := r_6 * r_6

	// Calculate the force magnitude using the Lennard-Jones potential derivative.
	forceMagnitude := (-12*A/r_12 + 6*B/r_6) / r

	// Scale the force vector and apply unit vector components.
	force := TriTuple{
		x: forceMagnitude * unitVector.x * 1e-5,
		y: forceMagnitude * unitVector.y * 1e-5,
		z: forceMagnitude * unitVector.z * 1e-5,
	}

	// Return the Lennard-Jones force vector.
	return force
}
