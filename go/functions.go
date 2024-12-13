// This file contains all the function related to the calculation of bond potential energy and corresponeded forces
package main

import (
	"fmt"
	"math"
)

// Calculate the potential energy due to Bond stretching
// Input: bond force constant, the distance between two atoms and equilibrium bond distance
// Output: a float64 of bond stretch energy
func CalculateBondStretchEnergy(k, r, r_0 float64) float64 {
	return 0.5 * k * (r*0.1 - r_0) * (r*0.1 - r_0)
}

// Calculate the potential energy due to angle bending
// Input: angle force constant, the angle of the three atoms and equilibrium angle
// Output: a float64 of angle bending energy
func CalculateAnglePotentialEnergy(k, theta, theta_0 float64) float64 {
	return 0.5 * k * (theta - theta_0) / 180 * math.Pi * (theta - theta_0) / 180 * math.Pi
}

// Calculate the potential energy due to proper dihedral angle bending
// Input: dihedral angle force constant, the dihedral angle of the four atoms, periodicity and phase angle
// Output: a float64 of dihedral angle bending energy
func CalculateProperDihedralAngleEnergy(kd, phi, pn, phase float64) float64 {
	return 0.5 * kd * (1 + math.Cos(pn*phi-phase/180*math.Pi))
}

// Calculate the angle in the three atoms
// Input: the pointers of the three atoms
// Output: a float64 of angle
func CalculateAngle(atom1, atom2, atom3 *Atom) float64 {
	// calculate the vectors of (atom2, atom1) and (atom2, atom3)
	vector1 := CalculateVector(atom1, atom2)
	vector2 := CalculateVector(atom3, atom2)

	// calculate the product of the two vectors
	upperValue := vector1.dot(vector2)

	// calculate the product of the two distance
	lowerValue := Distance(atom1.position, atom2.position) * Distance(atom2.position, atom3.position)

	// calculate the angle according to arccos((r_12 \cdot r_32) / |r_12| |r_32|)
	value := upperValue / lowerValue

	// The range of the angle is from 0 to pi, positive
	return math.Abs(math.Acos(value) / math.Pi * 180)
}

// Calculate the distance between two positions
// Input: Two Triptuples of positions
// Output: a float64 of distance
func Distance(p1, p2 TriTuple) float64 {
	deltaX := p1.x - p2.x
	deltaY := p1.y - p2.y
	deltaZ := p1.z - p2.z
	return math.Sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ)

}

// Calculate the dihedral angle in the four atoms
// Input: the pointers of the four atoms
// Output: a float64 ofdihedral angle
func CalculateDihedralAngle(atom1, atom2, atom3, atom4 *Atom) float64 {
	// calculate the three vectors of (atom1, atom2), (atom2, atom3) and (atom3, atom4)
	vector1 := CalculateVector(atom2, atom1)
	vector2 := CalculateVector(atom3, atom2)
	vector3 := CalculateVector(atom4, atom3)

	// build the planar vector
	plane1 := BuildNormalVector(vector1, vector2)
	plane2 := BuildNormalVector(vector2, vector3)

	// calculate the dihedral angle according to two planar vectors
	x := plane1.dot(plane2)
	y := magnitude(plane1) * magnitude(plane2)
	angle := math.Acos(x / y)

	// convert the unit of degree into radius
	return angle * (180 / math.Pi)
}

// Build the planar vector
// Input: two vectors in TriTuple on the same planar
// Output: the planar vector in TriTuple
func BuildNormalVector(vector1, vector2 TriTuple) TriTuple {
	var normVector TriTuple
	normVector.x = vector1.y*vector2.z - vector1.z*vector2.y
	normVector.y = vector1.z*vector2.x - vector1.x*vector2.z
	normVector.z = vector1.x*vector2.y - vector1.y*vector2.x

	return normVector
}

// The inner product of two vectors
// Input: two vectors in TriTuple
// Output: a float64 of inner product of two vectors
func (vector1 TriTuple) dot(vector2 TriTuple) float64 {
	return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z
}

// Calculate the magnitude of a vector
// Input: a vector in TriTuple
// Output: a float64 of the magnitude of a vector
func magnitude(vector TriTuple) float64 {
	return math.Sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z)
}

// Calculate the total potential energy of bonded and unbonded types with forces on each atom
// Input: A pointer of a protein and parameter datasets of unbonded, bonded and whether atoms are bound
// Output: A float64 of total energy and forceMap
func CombineEnergyAndForce(p *Protein, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (float64, map[int]*TriTuple) {
	// Calculate total energy and forces of bonded interactions
	bondedEnergy, bondedForceMap := CalculateTotalEnergyForce(p, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
	// Calculate total energy and forces of unbonded interactions
	unbondedEnergy, unbondedForceMap := CalculateTotalUnbondedEnergyForce(p, nonbondParameter)
	fmt.Println("bondedEnergy is:", bondedEnergy)
	fmt.Println("unbondedEnergy is:", unbondedEnergy)
	// Combine energies
	totalEnergy := bondedEnergy //+ unbondedEnergy
	// Create a total force map
	totalForceMap := make(map[int]*TriTuple)

	// Add forces due to bonded interactions
	for index, force := range bondedForceMap {

		totalForceMap[index] = &TriTuple{
			x: force.x,
			y: force.y,
			z: force.z,
		}
	}

	// Add forces due to unbonded interactions
	for index, force := range unbondedForceMap {
		if _, exists := totalForceMap[index]; exists {
			totalForceMap[index].x += force.x
			totalForceMap[index].y += force.y
			totalForceMap[index].z += force.z
		} else {
			totalForceMap[index] = &TriTuple{
				x: force.x,
				y: force.y,
				z: force.z,
			}
		}
	}

	return totalEnergy, totalForceMap
}

// Perform Energy minimization
// Input: A pointer of a protein and parameter datasets of unbonded, bonded and whether atoms are bound and iteration
// Output: A pointer of updated Protein
func PerformEnergyMinimization(currentProtein *Protein, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase, iteration int) *Protein {

	// set maximum displacement
	h := 0.01

	// range over the iteration
	for i := 0; i < iteration; i++ {
		// Combine energies and forces
		totalEnergy, totalForceMap := CombineEnergyAndForce(currentProtein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
		fmt.Printf("Iteration %d: Total Energy = %f\n", i, totalEnergy)

		tempProtein := CopyProtein(currentProtein)

		// Perform SteepestDescent, update positions in protein
		SteepestDescent(tempProtein, h, totalForceMap)

		// Calculate total energy of updated protein
		newTotalEnergy, _ := CombineEnergyAndForce(tempProtein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
		fmt.Printf("Iteration %d: New Total Energy = %f\n", i, newTotalEnergy)

		// If total energy decreases, accept the changes of positions and increase maximum displacement h
		// Otherwise, reject the changes in positions and decrease maximum displacement h
		if newTotalEnergy < totalEnergy {
			currentProtein = tempProtein
			h *= 1.2
		} else {
			h *= 0.2 * h
		}
	}

	return currentProtein
}

// Calculate the energies and forces due to the bonded interactions
// Input: A pointer of a protein and parameter datasets of bonded interactions and whether atoms are bound
// Output: A float64 of bonded energy and forceMap
func CalculateTotalEnergyForce(p *Protein, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (float64, map[int]*TriTuple) {
	forceMap := make(map[int]*TriTuple)

	bondEnergy := 0.0
	angleEnergy := 0.0
	dihedralEnergy := 0.0

	// range over each residue in protein
	for w, residue := range p.Residue {

		// Calculate bondstretch energy
		// range over each bond pair of the residue and check whether there are some atoms in the pairs
		for _, bondPairs := range residueParameterBondValue[residue.Name].bonds {

			// search the related atoms in the residue
			for i := 0; i < len(residue.Atoms)-1; i++ {
				atom1 := residue.Atoms[i]

				if atom1.element == (*bondPairs).atoms[0] {
					for j := i + 1; j < len(residue.Atoms); j++ {
						if residue.Atoms[j].element == (*bondPairs).atoms[1] {
							// find the two atoms in the currentbond pair of the residue
							atom2 := residue.Atoms[j]

							// Calculate the distance
							r := Distance(atom1.position, atom2.position)

							// if distance is zero, skip
							if r == 0 {
								continue
							}

							// if the distance larger than 2, considered as no interaction
							if r > 2 {
								continue
							}

							// search over the bond parameters of the two atoms
							parameterList := SearchParameter(2, bondParameter, atom1, atom2)

							// check whether the parameters exist
							if len(parameterList) != 1 {
								// calculate the bond force
								force := CalculateBondForce(parameterList[1], r, parameterList[0], atom1, atom2)

								// if Nan, skip
								if math.IsNaN(force.x) {
									continue
								}

								// add bond energy
								bondEnergy += CalculateBondStretchEnergy(parameterList[1], r, parameterList[0])

								// add forces
								_, exist := forceMap[atom1.index]
								if exist {
									forceMap[atom1.index].x += force.x
									forceMap[atom1.index].y += force.y
									forceMap[atom1.index].z += force.z
								} else {
									forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
									forceMap[atom1.index].x = force.x
									forceMap[atom1.index].y = force.y
									forceMap[atom1.index].z = force.z
								}

								_, exist1 := forceMap[atom2.index]
								if exist1 {
									forceMap[atom2.index].x += -force.x
									forceMap[atom2.index].y += -force.y
									forceMap[atom2.index].z += -force.z
								} else {
									forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
									forceMap[atom2.index].x = -force.x
									forceMap[atom2.index].y = -force.y
									forceMap[atom2.index].z = -force.z
								}

							}

						}
					}

					// find the peptide bond which connects current residue and the new residue where N in the next residue
					if (*bondPairs).atoms[1] == "N" && w != len(p.Residue)-1 {
						// search the next residue
						for k := range p.Residue[w+1].Atoms {
							if p.Residue[w+1].Atoms[k].element == "N" {

								atom2 := p.Residue[w+1].Atoms[k]
								r := Distance(atom1.position, atom2.position)

								// skip if ditance is zero
								if r == 0 {
									continue
								}

								// if the distance larger than 2, considered as no interaction
								if r > 2 {
									continue
								}

								// search over the bond parameters of the two atoms
								parameterList := []float64{0.13830, 354803.2}
								if len(parameterList) != 1 {

									// Calculate the force
									force := CalculateBondForce(parameterList[1], r, parameterList[0], atom1, atom2)

									// If Nan, skip
									if math.IsNaN(force.x) {
										continue
									}

									// add bond energy
									bondEnergy += CalculateBondStretchEnergy(parameterList[1], r, parameterList[0])

									// add force
									_, exist := forceMap[atom1.index]
									if exist {
										forceMap[atom1.index].x += force.x
										forceMap[atom1.index].y += force.y
										forceMap[atom1.index].z += force.z
									} else {
										forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom1.index].x = force.x
										forceMap[atom1.index].y = force.y
										forceMap[atom1.index].z = force.z
									}

									_, exist1 := forceMap[atom2.index]
									if exist1 {
										forceMap[atom2.index].x += -force.x
										forceMap[atom2.index].y += -force.y
										forceMap[atom2.index].z += -force.z
									} else {
										forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom2.index].x = -force.x
										forceMap[atom2.index].y = -force.y
										forceMap[atom2.index].z = -force.z
									}

								}

							}
						}
					}
				}
			}

		}

		// Calculate angle bending energy
		// range over each angle combinations of the residue and check whether there are some atoms in the angle combinations
		for _, angleTris := range residueParameterOtherValue[residue.Name].angles {
			// find the special angle of C, N, H or C, N ,CA which are in the peptide bond in the connection of current residue and the last residue
			// C here is in the last residue
			// "-" represent the last residue
			if (*angleTris).atoms[0][0] == '-' {
				if w != 0 {
					// search the las residue
					for a := range p.Residue[w-1].Atoms {
						if p.Residue[w-1].Atoms[a].element == "C" {
							atom1 := p.Residue[w-1].Atoms[a]

							// search current residue
							for i := 0; i < len(residue.Atoms)-1; i++ {
								if residue.Atoms[i].element == (*angleTris).atoms[1] {
									atom2 := residue.Atoms[i]
									for j := i + 1; j < len(residue.Atoms); j++ {
										if residue.Atoms[j].element == (*angleTris).atoms[2] {
											// find the three atoms
											atom3 := residue.Atoms[j]

											theta := CalculateAngle(atom1, atom2, atom3)

											// if Nan, skip
											if math.IsNaN(theta) {
												continue
											}

											// if the distance greater than 2, considered as no interaction
											if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 {
												continue
											}

											// the specific parameters
											parameterList := []float64{119.200, 418.400}
											if len(parameterList) != 1 {

												// Calculate the force
												force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

												// this angle energy
												OneAngleEnergy := CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
												// If Nan, skip
												if math.IsNaN(OneAngleEnergy) {
													continue
												}
												// add angle energy
												angleEnergy += OneAngleEnergy

												// add forces
												_, exist := forceMap[atom1.index]
												if exist {
													forceMap[atom1.index].x += force_i.x
													forceMap[atom1.index].y += force_i.y
													forceMap[atom1.index].z += force_i.z
												} else {
													forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom1.index].x = force_i.x
													forceMap[atom1.index].y = force_i.y
													forceMap[atom1.index].z = force_i.z
												}
												_, exist1 := forceMap[atom2.index]
												if exist1 {
													forceMap[atom2.index].x += force_j.x
													forceMap[atom2.index].y += force_j.y
													forceMap[atom2.index].z += force_j.z
												} else {
													forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom2.index].x = force_j.x
													forceMap[atom2.index].y = force_j.y
													forceMap[atom2.index].z = force_j.z
												}

												_, exist2 := forceMap[atom3.index]
												if exist2 {
													forceMap[atom3.index].x += force_k.x
													forceMap[atom3.index].y += force_k.y
													forceMap[atom3.index].z += force_k.z
												} else {
													forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom3.index].x = force_k.x
													forceMap[atom3.index].y = force_k.y
													forceMap[atom3.index].z = force_k.z
												}

											}

										}
									}
								}
							}
						}
					}
				}
			}

			// range over each angle combinations of the current residue and check whether there are some atoms in the combination
			for i := 0; i < len(residue.Atoms)-2; i++ {
				atom1 := residue.Atoms[i]
				if atom1.element == (*angleTris).atoms[0] {
					for j := i + 1; j < len(residue.Atoms)-1; j++ {
						if residue.Atoms[j].element == (*angleTris).atoms[1] {
							atom2 := residue.Atoms[j]
							for k := j + 1; k < len(residue.Atoms); k++ {
								if residue.Atoms[k].element == (*angleTris).atoms[2] {
									// find the three atoms
									atom3 := residue.Atoms[k]
									theta := CalculateAngle(atom1, atom2, atom3)

									// If Nan, skip
									if math.IsNaN(theta) {
										continue
									}

									// If the distance larger than 2, considered as no interaction
									if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 {
										continue
									}

									// search the angle parameters
									parameterList := SearchParameter(3, angleParameter, atom1, atom2, atom3)
									if len(parameterList) != 1 {

										// calculate the force
										force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

										// this angle energy
										OneAngleEnergy := CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
										// If Nan, skip
										if math.IsNaN(OneAngleEnergy) {
											continue
										}
										// add angle energy
										angleEnergy += OneAngleEnergy

										// add forces
										_, exist := forceMap[atom1.index]
										if exist {
											forceMap[atom1.index].x += force_i.x
											forceMap[atom1.index].y += force_i.y
											forceMap[atom1.index].z += force_i.z
										} else {
											forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[atom1.index].x = force_i.x
											forceMap[atom1.index].y = force_i.y
											forceMap[atom1.index].z = force_i.z
										}
										_, exist1 := forceMap[atom2.index]
										if exist1 {
											forceMap[atom2.index].x += force_j.x
											forceMap[atom2.index].y += force_j.y
											forceMap[atom2.index].z += force_j.z
										} else {
											forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[atom2.index].x = force_j.x
											forceMap[atom2.index].y = force_j.y
											forceMap[atom2.index].z = force_j.z
										}

										_, exist2 := forceMap[atom3.index]
										if exist2 {
											forceMap[atom3.index].x += force_k.x
											forceMap[atom3.index].y += force_k.y
											forceMap[atom3.index].z += force_k.z
										} else {
											forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[atom3.index].x = force_k.x
											forceMap[atom3.index].y = force_k.y
											forceMap[atom3.index].z = force_k.z
										}
									}

								}
							}

							// find the special angle of O, C, N which are in the peptide bond in the connection of current residue and the next residue
							if (*angleTris).atoms[2][0] == '+' && w != len(p.Residue)-1 {
								// search the next residue
								for k := range p.Residue[w+1].Atoms {
									if p.Residue[w+1].Atoms[k].element == "N" {
										// find the three atoms
										atom3 := p.Residue[w+1].Atoms[k]
										theta := CalculateAngle(atom1, atom2, atom3)

										// if Nan, skip
										if math.IsNaN(theta) {
											continue
										}

										// if the distance larger than 2, considered as no interaction
										if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 {
											continue
										}

										// the specific parameters
										parameterList := []float64{120.900, 669.440}

										if len(parameterList) != 1 {

											// calculate the force
											force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

											// this angle energy
											OneAngleEnergy := CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
											// If Nan, skip
											if math.IsNaN(OneAngleEnergy) {
												continue
											}

											// add energy
											angleEnergy += OneAngleEnergy

											// add forces
											_, exist := forceMap[atom1.index]
											if exist {
												forceMap[atom1.index].x += force_i.x
												forceMap[atom1.index].y += force_i.y
												forceMap[atom1.index].z += force_i.z
											} else {
												forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom1.index].x = force_i.x
												forceMap[atom1.index].y = force_i.y
												forceMap[atom1.index].z = force_i.z
											}

											_, exist1 := forceMap[atom2.index]
											if exist1 {
												forceMap[atom2.index].x += force_j.x
												forceMap[atom2.index].y += force_j.y
												forceMap[atom2.index].z += force_j.z
											} else {
												forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom2.index].x = force_j.x
												forceMap[atom2.index].y = force_j.y
												forceMap[atom2.index].z = force_j.z
											}

											_, exist2 := forceMap[atom3.index]
											if exist2 {
												forceMap[atom3.index].x += force_k.x
												forceMap[atom3.index].y += force_k.y
												forceMap[atom3.index].z += force_k.z
											} else {
												forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom3.index].x = force_k.x
												forceMap[atom3.index].y = force_k.y
												forceMap[atom3.index].z = force_k.z
											}

										}
									}
								}

							}
						}
					}

				}
			}
		}

		// Calculate dihedral angle bending energy
		// range over each combination of the residue and check whether there are some atoms in the combination
		for _, dihedralValues := range residueParameterOtherValue[residue.Name].dihedrals {

			// find the special angle of -C, N,CA,C   which are in the peptide bond in the connection of current residue and the last residue
			// "-" represents in the last residue
			if (*dihedralValues).atoms[0] == "-CA" && w != 0 {
				// search the las residue
				for i := range p.Residue[w-1].Atoms {
					if p.Residue[w-1].Atoms[i].element == "CA" {
						atom1 := p.Residue[w-1].Atoms[i]
						atom2 := p.Residue[w-1].Atoms[i+1]

						for j := range p.Residue[w].Atoms {
							if p.Residue[w].Atoms[j].element == (*dihedralValues).atoms[2] {
								// find -C, N,CA,C
								atom3 := p.Residue[w].Atoms[j]
								atom4 := p.Residue[w].Atoms[j+1]
								phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)

								// if Nan, skip
								if math.IsNaN(phi) {
									continue
								}

								// if the distance larger than 2, considered as no interaction
								if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 || Distance(atom3.position, atom4.position) > 2 {
									continue
								}

								// special parameters
								parameterList := []float64{180.0, 6.06680, 2}
								if len(parameterList) != 1 {

									// calculate the forces
									force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

									oneDihedralEnergy := CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])

									// If Nan, skip
									if math.IsNaN(oneDihedralEnergy) {
										continue
									}

									// add energy
									dihedralEnergy += oneDihedralEnergy

									// add forces
									_, exist := forceMap[atom1.index]
									if exist {
										forceMap[atom1.index].x += force_i.x
										forceMap[atom1.index].y += force_i.y
										forceMap[atom1.index].z += force_i.z
									} else {
										forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom1.index].x = force_i.x
										forceMap[atom1.index].y = force_i.y
										forceMap[atom1.index].z = force_i.z
									}

									_, exist1 := forceMap[atom2.index]
									if exist1 {
										forceMap[atom2.index].x += force_j.x
										forceMap[atom2.index].y += force_j.y
										forceMap[atom2.index].z += force_j.z
									} else {
										forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom2.index].x = force_j.x
										forceMap[atom2.index].y = force_j.y
										forceMap[atom2.index].z = force_j.z
									}

									_, exist2 := forceMap[atom3.index]
									if exist2 {
										forceMap[atom3.index].x += force_k.x
										forceMap[atom3.index].y += force_k.y
										forceMap[atom3.index].z += force_k.z
									} else {
										forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom3.index].x = force_k.x
										forceMap[atom3.index].y = force_k.y
										forceMap[atom3.index].z = force_k.z
									}

									_, exist3 := forceMap[atom4.index]
									if exist3 {
										forceMap[atom4.index].x += force_l.x
										forceMap[atom4.index].y += force_l.y
										forceMap[atom4.index].z += force_l.z
									} else {
										forceMap[atom4.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[atom4.index].x = force_l.x
										forceMap[atom4.index].y = force_l.y
										forceMap[atom4.index].z = force_l.z
									}

								}

							}
						}
					}
				}
			}

			// find the special angle of -CA, -C, N, CA   which are in the peptide bond in the connection of current residue and the last residue
			// "-" represents in the last residue
			if (*dihedralValues).atoms[0] == "-C" && w != 0 {
				for i := range p.Residue[w-1].Atoms {
					if p.Residue[w-1].Atoms[i].element == "C" {
						// find -CA, -C, N, CA
						atom1 := p.Residue[w-1].Atoms[i]
						atom2 := residue.Atoms[0]
						atom3 := residue.Atoms[1]
						atom4 := residue.Atoms[2]
						phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)
						// if Nan, skip
						if math.IsNaN(phi) {
							continue
						}

						// If the distance larger than 2, considered as no interaction
						if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 || Distance(atom3.position, atom4.position) > 2 {
							continue
						}

						// search the parameters
						parameterList := SearchParameter(4, dihedralParameter, atom1, atom2, atom3, atom4)

						if len(parameterList) != 1 {

							// Calculate the forces
							force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

							oneDihedralEnergy := CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])

							// If Nan, skip
							if math.IsNaN(oneDihedralEnergy) {
								continue
							}

							// add energy
							dihedralEnergy += oneDihedralEnergy

							// add force
							_, exist := forceMap[atom1.index]
							if exist {
								forceMap[atom1.index].x += force_i.x
								forceMap[atom1.index].y += force_i.y
								forceMap[atom1.index].z += force_i.z
							} else {
								forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[atom1.index].x = force_i.x
								forceMap[atom1.index].y = force_i.y
								forceMap[atom1.index].z = force_i.z
							}

							_, exist1 := forceMap[atom2.index]
							if exist1 {
								forceMap[atom2.index].x += force_j.x
								forceMap[atom2.index].y += force_j.y
								forceMap[atom2.index].z += force_j.z
							} else {
								forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[atom2.index].x = force_j.x
								forceMap[atom2.index].y = force_j.y
								forceMap[atom2.index].z = force_j.z
							}

							_, exist2 := forceMap[atom3.index]
							if exist2 {
								forceMap[atom3.index].x += force_k.x
								forceMap[atom3.index].y += force_k.y
								forceMap[atom3.index].z += force_k.z
							} else {
								forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[atom3.index].x = force_k.x
								forceMap[atom3.index].y = force_k.y
								forceMap[atom3.index].z = force_k.z
							}

							_, exist3 := forceMap[atom4.index]
							if exist3 {
								forceMap[atom4.index].x += force_l.x
								forceMap[atom4.index].y += force_l.y
								forceMap[atom4.index].z += force_l.z
							} else {
								forceMap[atom4.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[atom4.index].x = force_l.x
								forceMap[atom4.index].y = force_l.y
								forceMap[atom4.index].z = force_l.z
							}

						}
					}
				}
			}

			// range over each dihedral angle combinations of the residue and check whether there are some atoms in the dihedral angle combinations
			for i := 0; i < len(residue.Atoms)-3; i++ {
				atom1 := residue.Atoms[i]
				if atom1.element == (*dihedralValues).atoms[0] {
					for j := i + 1; j < len(residue.Atoms)-2; j++ {
						if residue.Atoms[j].element == (*dihedralValues).atoms[1] {
							atom2 := residue.Atoms[j]
							for k := j + 1; k < len(residue.Atoms)-1; k++ {
								if residue.Atoms[k].element == (*dihedralValues).atoms[2] {
									atom3 := residue.Atoms[k]
									for l := k + 1; l < len(residue.Atoms); l++ {
										if residue.Atoms[l].element == (*dihedralValues).atoms[3] {
											// find the four atoms
											atom4 := residue.Atoms[l]
											phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)

											// if Nan, skip
											if math.IsNaN(phi) {
												continue
											}

											// if the distance larger than 2, considered as no interaction
											if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 || Distance(atom3.position, atom4.position) > 2 {
												continue
											}

											// search parameters
											parameterList := SearchParameter(4, dihedralParameter, atom1, atom2, atom3, atom4)

											// if parameters exist
											if len(parameterList) != 1 {

												// calculate force
												force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

												oneDihedralEnergy := CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])

												// If Nan, skip
												if math.IsNaN(oneDihedralEnergy) {
													continue
												}

												// add energy
												dihedralEnergy += oneDihedralEnergy

												// add forces
												_, exist := forceMap[atom1.index]
												if exist {
													forceMap[atom1.index].x += force_i.x
													forceMap[atom1.index].y += force_i.y
													forceMap[atom1.index].z += force_i.z
												} else {
													forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom1.index].x = force_i.x
													forceMap[atom1.index].y = force_i.y
													forceMap[atom1.index].z = force_i.z
												}

												_, exist1 := forceMap[atom2.index]
												if exist1 {
													forceMap[atom2.index].x += force_j.x
													forceMap[atom2.index].y += force_j.y
													forceMap[atom2.index].z += force_j.z
												} else {
													forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom2.index].x = force_j.x
													forceMap[atom2.index].y = force_j.y
													forceMap[atom2.index].z = force_j.z
												}

												_, exist2 := forceMap[atom3.index]
												if exist2 {
													forceMap[atom3.index].x += force_k.x
													forceMap[atom3.index].y += force_k.y
													forceMap[atom3.index].z += force_k.z
												} else {
													forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom3.index].x = force_k.x
													forceMap[atom3.index].y = force_k.y
													forceMap[atom3.index].z = force_k.z
												}

												_, exist3 := forceMap[atom4.index]
												if exist3 {
													forceMap[atom4.index].x += force_l.x
													forceMap[atom4.index].y += force_l.y
													forceMap[atom4.index].z += force_l.z
												} else {
													forceMap[atom4.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[atom4.index].x = force_l.x
													forceMap[atom4.index].y = force_l.y
													forceMap[atom4.index].z = force_l.z
												}

											}
										}
									}

									// find the special angle of  N, CA, C, +N   which are in the peptide bond in the connection of current residue and the next residue
									// "+" represents in the next residue
									if (*dihedralValues).atoms[3] == "+N" && w != len(p.Residue)-1 {
										atom4 := p.Residue[w+1].Atoms[0]
										phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)

										// if Nan, skip
										if math.IsNaN(phi) {
											continue
										}

										// if distance larger than 2, considered as no interaction
										if Distance(atom3.position, atom2.position) > 2 || Distance(atom2.position, atom1.position) > 2 || Distance(atom3.position, atom4.position) > 2 {
											continue
										}
										// the specific parameters
										parameterList := []float64{180.0, 15.167, 2}
										if len(parameterList) != 1 {
											// calculate force
											force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

											oneDihedralEnergy := CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])

											// If Nan, skip
											if math.IsNaN(oneDihedralEnergy) {
												continue
											}

											// add energy
											dihedralEnergy += oneDihedralEnergy

											// add forces
											_, exist := forceMap[atom1.index]
											if exist {
												forceMap[atom1.index].x += force_i.x
												forceMap[atom1.index].y += force_i.y
												forceMap[atom1.index].z += force_i.z
											} else {
												forceMap[atom1.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom1.index].x = force_i.x
												forceMap[atom1.index].y = force_i.y
												forceMap[atom1.index].z = force_i.z
											}

											_, exist1 := forceMap[atom2.index]
											if exist1 {
												forceMap[atom2.index].x += force_j.x
												forceMap[atom2.index].y += force_j.y
												forceMap[atom2.index].z += force_j.z
											} else {
												forceMap[atom2.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom2.index].x = force_j.x
												forceMap[atom2.index].y = force_j.y
												forceMap[atom2.index].z = force_j.z
											}

											_, exist2 := forceMap[atom3.index]
											if exist2 {
												forceMap[atom3.index].x += force_k.x
												forceMap[atom3.index].y += force_k.y
												forceMap[atom3.index].z += force_k.z
											} else {
												forceMap[atom3.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom3.index].x = force_k.x
												forceMap[atom3.index].y = force_k.y
												forceMap[atom3.index].z = force_k.z
											}

											_, exist3 := forceMap[atom4.index]
											if exist3 {
												forceMap[atom4.index].x += force_l.x
												forceMap[atom4.index].y += force_l.y
												forceMap[atom4.index].z += force_l.z
											} else {
												forceMap[atom4.index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
												forceMap[atom4.index].x = force_l.x
												forceMap[atom4.index].y = force_l.y
												forceMap[atom4.index].z = force_l.z
											}

										}

									}
								}
							}

						}
					}
				}
			}

		}

	}

	totalEnergy := bondEnergy + angleEnergy + dihedralEnergy
	return totalEnergy, forceMap
}

// search the parameters in the paramter dataset
// Input: value represent the number of atoms, related parameter dataset and pointers of atoms
func SearchParameter(value int, parameterData parameterDatabase, atoms ...*Atom) []float64 {
	// search parameter combinations in the Descending order
	for i := range parameterData.atomPair {
		sym := 0
		for j := range parameterData.atomPair[i].atomName {
			// 'X' can be any atom
			if parameterData.atomPair[i].atomName[j] == "X" {
				sym += 1
				continue
			}

			// consider H and ignore the specific name
			if parameterData.atomPair[i].atomName[j][0] == 'H' && atoms[j].element[0] == 'H' {
				sym += 1
				continue
			}

			// eliminate te influence of "*" in N which represent atom in the next residue
			if parameterData.atomPair[i].atomName[j] == "N*" && atoms[j].element == "N" {
				sym += 1
				continue
			}

			// eliminate te influence of "*" in C which represent atom in the last residue
			if parameterData.atomPair[i].atomName[j] == "C*" && atoms[j].element == "C" {
				sym += 1
				continue
			}

			// if the atom name mismatch, break the loop
			if atoms[j].element[0] != parameterData.atomPair[i].atomName[j][0] {
				break
			}
			// add 1 to indicate the match
			sym += 1
		}

		// if all the atom matched, return the paramters
		if sym == value {
			return parameterData.atomPair[i].parameter
		}

	}

	// search parameter combinations in the ascending order
	for i := range parameterData.atomPair {
		sym := 0
		for j := range parameterData.atomPair[i].atomName {
			// 'X' can be any atom
			if parameterData.atomPair[i].atomName[value-j-1] == "X" {
				sym += 1
				continue
			}

			// consider H and ignore the specific name
			if parameterData.atomPair[i].atomName[value-j-1][0] == 'H' && atoms[j].element[0] == 'H' {
				sym += 1
				continue
			}

			// eliminate te influence of "*" in N which represent atom in the next residue
			if parameterData.atomPair[i].atomName[value-j-1] == "N*" && atoms[j].element == "N" {
				sym += 1
				continue
			}

			// eliminate te influence of "*" in C which represent atom in the last residue
			if parameterData.atomPair[i].atomName[value-j-1] == "C*" && atoms[j].element == "C" {
				sym += 1
				continue
			}

			// if the atom name mismatch, break the loop
			if atoms[j].element[0] != parameterData.atomPair[i].atomName[value-j-1][0] {
				break
			}
			// add 1 to indicate the match
			sym += 1
		}

		// if all the atom matched, return the paramters
		if sym == value {
			return parameterData.atomPair[i].parameter
		}

	}

	// mismatch, only return one parameter which will be detected and eliminated
	return []float64{0.0}
}

// Perform position update in Steepest Descent
// Input: a pointer of protein, displacement and forceMap
// Output: a pointer of protein
func SteepestDescent(protein *Protein, h float64, forceMap map[int]*TriTuple) *Protein {
	// range over residues
	for i := range protein.Residue {
		// range over atoms
		for j := range protein.Residue[i].Atoms {
			// check the existence of force
			_, exist := forceMap[protein.Residue[i].Atoms[j].index+1]
			if exist {
				force := forceMap[protein.Residue[i].Atoms[j].index+1]

				// skip nan or zero force
				magn := magnitude(*force)
				if magn == 0 || math.IsNaN(magn) {
					continue
				}

				// update position according to p_new = p_pre + force/ magnitude of force * displacement
				protein.Residue[i].Atoms[j].position.x = protein.Residue[i].Atoms[j].position.x + (force.x*h)/magn
				protein.Residue[i].Atoms[j].position.y = protein.Residue[i].Atoms[j].position.y + (force.y*h)/magn
				protein.Residue[i].Atoms[j].position.z = protein.Residue[i].Atoms[j].position.z + (force.z*h)/magn
			}

		}
	}

	return protein

}

// Calculate the bond force
// Input: bond force constant, the distance between two atoms and equilibrium bond distance
// Output: a float64 of bond stretch force
func CalculateBondForce(k, r, r_0 float64, atom1, atom2 *Atom) TriTuple {
	// Calculate the distance
	bondLen := Distance(atom1.position, atom2.position)

	// calculate the unit vector
	unitVector := TriTuple{
		x: -(atom1.position.x - atom2.position.x) / bondLen,
		y: -(atom1.position.y - atom2.position.y) / bondLen,
		z: -(atom1.position.z - atom2.position.z) / bondLen,
	}

	// calculate the total force value
	fScale := k * (r - r_0*10)

	// calculate the force and convert the unit into g * A / (mol* fs^2)
	force := TriTuple{
		x: fScale * unitVector.x / math.Pow10(6),
		y: fScale * unitVector.y / math.Pow10(6),
		z: fScale * unitVector.z / math.Pow10(6),
	}

	return force
}

// Calculate the angle force
// Input: angle force constant, the angle of the three atoms and equilibrium angle
// Output: a float64 of angle force
func CalculateAngleForce(k, theta, theta_0 float64, atom1, atom2, atom3 *Atom) (TriTuple, TriTuple, TriTuple) {
	// the previous detrivatives between potential energy and theta, cos theta and theta
	der_that_cos := (-1) * (1 / math.Sin(theta/180*math.Pi))
	der_U_thate := k * (theta - theta_0) / 180 * math.Pi

	// if Nan, skip
	if math.IsNaN(der_that_cos) {
		return TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}
	}

	// the derivates of theta and position
	der_theta_x_12 := DerivateAnglePositionX(atom1, atom2, atom3, theta)
	der_theta_x_32 := DerivateAnglePositionX(atom3, atom2, atom1, theta)

	der_theta_y_12 := DerivateAnglePositionY(atom1, atom2, atom3, theta)
	der_theta_y_32 := DerivateAnglePositionY(atom3, atom2, atom1, theta)

	der_theta_z_12 := DerivateAnglePositionZ(atom1, atom2, atom3, theta)
	der_theta_z_32 := DerivateAnglePositionZ(atom3, atom2, atom1, theta)

	// calculate the forces and convert the unit into g*A/ (mol*fs^2)
	force_i := TriTuple{
		x: der_U_thate * der_that_cos * der_theta_x_12 / math.Pow10(4) * 2 * math.Pi,
		y: der_U_thate * der_that_cos * der_theta_y_12 / math.Pow10(4) * 2 * math.Pi,
		z: der_U_thate * der_that_cos * der_theta_z_12 / math.Pow10(4) * 2 * math.Pi,
	}

	force_k := TriTuple{
		x: der_U_thate * der_that_cos * der_theta_x_32 / math.Pow10(4) * 2 * math.Pi,
		y: der_U_thate * der_that_cos * der_theta_y_32 / math.Pow10(4) * 2 * math.Pi,
		z: der_U_thate * der_that_cos * der_theta_z_32 / math.Pow10(4) * 2 * math.Pi,
	}

	force_j := TriTuple{
		x: (-force_i.x - force_k.x) * 2 * math.Pi,
		y: (-force_i.y - force_k.y) * 2 * math.Pi,
		z: (-force_i.z - force_k.z) * 2 * math.Pi,
	}

	// if Nan skip
	if math.IsNaN(force_i.x) || math.IsNaN(force_j.x) || math.IsNaN(force_k.x) {
		return TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}
	}

	return force_i, force_j, force_k

}

// Calculate the angle derivatives between theta and X positions
// Input: pointers of three atoms and their theta
// Output: a float64 angle derivatives between theta and X positions
func DerivateAnglePositionX(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.x-atom2.position.x)/Distance(atom3.position, atom2.position) - (atom1.position.x-atom2.position.x)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

// Calculate the angle derivatives between theta and Y positions
// Input: pointers of three atoms and their theta
// Output: a float64 angle derivatives between theta and Y positions
func DerivateAnglePositionY(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.y-atom2.position.y)/Distance(atom3.position, atom2.position) - (atom1.position.y-atom2.position.y)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

// Calculate the angle derivatives between theta and Z positions
// Input: pointers of three atoms and their theta
// Output: a float64 angle derivatives between theta and Z positions
func DerivateAnglePositionZ(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.z-atom2.position.z)/Distance(atom3.position, atom2.position) - (atom1.position.z-atom2.position.z)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

// Calculate the proper dihedral force
// Input: dihedral angle force constant, the dihedral angle of the four atoms, periodicity and phase angle
// Output: a float64 of dihedral angle force
func CalculateProperDihedralsForce(kd, phi, pn, phase float64, atom1, atom2, atom3, atom4 *Atom) (TriTuple, TriTuple, TriTuple, TriTuple) {
	// the previous detrivatives between potential energy and phi, cos phi and phi
	der_U_phi := -0.5 * kd * pn * math.Sin((pn*phi - phase/180*math.Pi))
	der_phi_cos := -1 / math.Sin(phi)

	// if Nan, skip
	if math.IsNaN(der_phi_cos) || math.IsNaN(der_U_phi) {
		return TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}
	}

	// calculate the vector
	vector12 := CalculateVector(atom2, atom1)
	vector32 := CalculateVector(atom2, atom3)
	vector43 := CalculateVector(atom3, atom4)

	// calculate the planar vector
	v_t := Cross(vector12, vector32)
	v_u := Cross(vector43, vector32)

	// calculate the derivate of cos phi with planar vector
	der_cos_tx, der_cos_ty, der_cos_tz := CalculateDerivate(v_t, v_u, phi)
	der_cos_ux, der_cos_uy, der_cos_uz := CalculateDerivate(v_u, v_t, phi)

	// calculate the force, convert the unit in g*A/(mol*fs^2)
	force_i := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*(-vector32.z) + der_cos_tz*vector32.y) / math.Pow10(4) * 2 * math.Pi,
		y: der_U_phi * der_phi_cos * (der_cos_tz*(-vector32.x) + der_cos_tx*vector32.z) / math.Pow10(4) * 2 * math.Pi,
		z: der_U_phi * der_phi_cos * (der_cos_tx*(-vector32.y) + der_cos_ty*vector32.x) / math.Pow10(4) * 2 * math.Pi,
	}

	force_j := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*(-vector12.z+vector32.z) + der_cos_tz*(-vector32.y+vector12.y) + der_cos_uy*(-vector43.z) + der_cos_uz*(vector43.y)) / math.Pow10(4) * 2 * math.Pi,
		y: der_U_phi * der_phi_cos * (der_cos_tz*(-vector12.x+vector32.x) + der_cos_tx*(-vector32.z+vector12.z) + der_cos_uy*(vector43.z) + der_cos_uz*(-vector43.x)) / math.Pow10(4) * 2 * math.Pi,
		z: der_U_phi * der_phi_cos * (der_cos_tx*(-vector12.y+vector32.y) + der_cos_ty*(-vector32.x+vector12.x) + der_cos_uy*(-vector43.y) + der_cos_uz*(vector43.x)) / math.Pow10(4) * 2 * math.Pi,
	}

	force_k := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*vector12.z - der_cos_tz*vector12.y + der_cos_uy*(vector32.z+vector43.z) - der_cos_uz*(vector32.y+vector43.y)) / math.Pow10(4) * 2 * math.Pi,
		y: der_U_phi * der_phi_cos * (der_cos_tz*vector12.x - der_cos_tx*vector12.z + der_cos_uz*(vector32.x+vector43.x) - der_cos_ux*(vector32.z+vector43.z)) / math.Pow10(4) * 2 * math.Pi,
		z: der_U_phi * der_phi_cos * (der_cos_tx*vector12.y - der_cos_ty*vector12.z + der_cos_ux*(vector32.y+vector43.y) - der_cos_uy*(vector32.x+vector43.x)) / math.Pow10(4) * 2 * math.Pi,
	}

	force_l := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_uy*(-vector32.z) + der_cos_uz*vector32.y) / math.Pow10(4) * 2 * math.Pi,
		y: der_U_phi * der_phi_cos * (der_cos_uz*(-vector32.x) + der_cos_ux*vector32.z) / math.Pow10(4) * 2 * math.Pi,
		z: der_U_phi * der_phi_cos * (der_cos_ux*(-vector32.y) + der_cos_uy*vector32.x) / math.Pow10(4) * 2 * math.Pi,
	}

	return force_i, force_j, force_k, force_l

}

// calculate the derivates of cos phi to planar vector
// Input: two vectors in TriTuple and phi
// Output: three float64 of derivates
func CalculateDerivate(v1, v2 TriTuple, phi float64) (float64, float64, float64) {
	return (1 / magnitude(v1)) * (v2.x/magnitude(v2) - v1.x/magnitude(v1)*math.Cos(phi)),
		(1 / magnitude(v1)) * (v2.y/magnitude(v2) - v1.y/magnitude(v1)*math.Cos(phi)),
		(1 / magnitude(v1)) * (v2.z/magnitude(v2) - v1.z/magnitude(v1)*math.Cos(phi))
}

// Calculate the
func Cross(v1, v2 TriTuple) TriTuple {
	return TriTuple{
		x: v1.x * v2.x,
		y: v1.y * v2.y,
		z: v1.z * v2.z,
	}
}

// Copy protein in the pointer
// Input: a pointer of a protein
// Output: a new pointer of a new protein with the same value of the input
func CopyProtein(currentProtein *Protein) *Protein {
	var newProtein Protein
	newProtein.Name = currentProtein.Name

	newProtein.Residue = make([]*Residue, len(currentProtein.Residue))
	for i := range currentProtein.Residue {
		newProtein.Residue[i] = CopyResidue(currentProtein.Residue[i])
	}

	return &newProtein
}

// Copy residue
// Input: a pointer of a residue
// Output: a new pointer of a new residue with the same of the input
func CopyResidue(currRes *Residue) *Residue {
	var newRes Residue
	newRes.Name = currRes.Name
	newRes.ID = currRes.ID
	newRes.ChainID = currRes.ChainID

	newRes.Atoms = make([]*Atom, len(currRes.Atoms))
	for i := range currRes.Atoms {
		newRes.Atoms[i] = CopyAtom(currRes.Atoms[i])
	}

	return &newRes
}

// Copy atom
// Input: a pointer of an atom
// Output: a new pointer of a new atom with the same value of the input atom
func CopyAtom(currAtom *Atom) *Atom {
	var newAtom Atom
	newAtom.mass = currAtom.mass
	newAtom.force = CopyTriTuple(currAtom.force)
	newAtom.position = CopyTriTuple(currAtom.position)
	newAtom.velocity = CopyTriTuple(currAtom.velocity)
	newAtom.accelerated = CopyTriTuple(currAtom.accelerated)
	newAtom.element = currAtom.element
	newAtom.charge = currAtom.charge
	newAtom.index = currAtom.index
	return &newAtom
}

// Copy TriTuple
// Input: a TriTuple
// Output: a copied TriTuple
func CopyTriTuple(tri TriTuple) TriTuple {
	var newTri TriTuple
	newTri.x = tri.x
	newTri.y = tri.y
	newTri.z = tri.z

	return newTri
}

// build a vector between two atoms
// Input: pointers of two atoms
// Output: vector in TriTuple
func CalculateVector(atom1, atom2 *Atom) TriTuple {
	var vector TriTuple
	vector.x = (atom2.position.x - atom1.position.x)
	vector.y = (atom2.position.y - atom1.position.y)
	vector.z = (atom2.position.z - atom1.position.z)

	return vector
}
