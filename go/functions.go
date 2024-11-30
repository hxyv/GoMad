package main

import (
	"fmt"
	"math"
)

func CalculateBondStretchEnergy(k, r, r_0 float64) float64 {
	return 0.5 * k * (r*0.1 - r_0) * (r*0.1 - r_0)
}

func CalculateAnglePotentialEnergy(k, theta, theta_0 float64) float64 {
	return 0.5 * k * (theta - theta_0) / 360 * math.Pi * (theta - theta_0) / 360 * math.Pi
}

func CalculateProperDihedralAngleEnergy(kd, phi, pn, phase float64) float64 {
	return 0.5 * kd * (1 + math.Cos(pn*phi-phase/360*math.Pi))
}

func CalculateAngle(atom1, atom2, atom3 *Atom) float64 {
	vector1 := CalculateVector(atom2, atom1)
	vector2 := CalculateVector(atom3, atom1)

	upperValue := vector1.dot(vector2)
	lowerValue := Distance(atom1.position, atom2.position) * Distance(atom2.position, atom3.position)

	value := upperValue / lowerValue

	return math.Acos(value) / math.Pi * 360
}

func Distance(p1, p2 TriTuple) float64 {
	deltaX := p1.x - p2.x
	deltaY := p1.y - p2.y
	deltaZ := p1.z - p2.z
	return math.Sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ)

}

func CalculateDihedralAngle(atom1, atom2, atom3, atom4 *Atom) float64 {
	vector1 := CalculateVector(atom2, atom1)
	vector2 := CalculateVector(atom3, atom2)
	vector3 := CalculateVector(atom4, atom3)

	plane1 := BuildNormalVector(vector1, vector2)
	plane2 := BuildNormalVector(vector2, vector3)
	plane3 := BuildNormalVector(plane1, plane2)

	x := plane1.dot(plane2)
	y := plane3.dot(vector2) / magnitude(vector2)
	angle := math.Atan2(y, x)

	return angle * (180 / math.Pi)
}

func BuildNormalVector(vector1, vector2 TriTuple) TriTuple {
	var normVector TriTuple
	normVector.x = vector1.y*vector2.z - vector1.z*vector2.y
	normVector.y = vector1.z*vector2.x - vector1.x*vector2.z
	normVector.z = vector1.x*vector2.y - vector1.y*vector2.x

	return normVector
}

func (vector1 TriTuple) dot(vector2 TriTuple) float64 {
	return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z
}

func magnitude(vector TriTuple) float64 {
	return math.Sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z)
}

func CombineEnergyAndForce(p *Protein, residueParameterValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (float64, map[int]*TriTuple) {
	// Calculate total energy and forces of bonded interactions
	bondedEnergy, bondedForceMap := CalculateTotalEnergyForce(p, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)

	// Calculate total energy and forces of unbonded interactions
	unbondedEnergy, unbondedForceMap := CalculateTotalUnbondedEnergyForce(p, nonbondParameter)
	// Combine energies
	totalEnergy := bondedEnergy + unbondedEnergy

	// Create a total force map
	totalForceMap := make(map[int]*TriTuple)
	for index, force := range bondedForceMap {
		totalForceMap[index] = &TriTuple{
			x: force.x,
			y: force.y,
			z: force.z,
		}
	}

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

func PerformEnergyMinimization(currentProtein *Protein, residueParameterValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) *Protein {
	iteration := 5
	// set maximum displacement
	h := 0.01

	for i := 0; i < iteration; i++ {
		// Combine energies and forces
		totalEnergy, totalForceMap := CombineEnergyAndForce(currentProtein, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
		fmt.Printf("Iteration %d: Total Energy = %f\n", i, totalEnergy)

		tempProtein := CopyProtein(currentProtein)

		// Perform SteepestDescent, update positions in protein
		SteepestDescent(tempProtein, h, totalForceMap)

		// Calculate total energy of updated protein
		newTotalEnergy, _ := CombineEnergyAndForce(tempProtein, residueParameterValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
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

func CalculateTotalEnergyForce(p *Protein, residueParameterValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (float64, map[int]*TriTuple) {
	forceMap := make(map[int]*TriTuple)

	// range over each residue in protein

	bondEnergy := 0.0
	angleEnergy := 0.0
	dihedralEnergy := 0.0

	index := 0
	for _, residue := range p.Residue {

		// Calculate bondstretch energy
		// range over each bond

		for _, bondPairs := range residueParameterValue[residue.Name].bonds {
			for i := 0; i < len(residue.Atoms)-1; i++ {
				atom1 := residue.Atoms[i]
				if atom1.element[0] == (*bondPairs).atoms[0][0] {
					for j := i + 1; j < len(residue.Atoms); j++ {
						if residue.Atoms[j].element[0] == (*bondPairs).atoms[1][0] {

							atom2 := residue.Atoms[j]
							r := Distance(atom1.position, atom2.position)
							if r == 0 {
								continue
							}
							if r > 2 {
								continue
							}

							parameterList := SearchParameter(2, bondParameter, atom1, atom2)
							if len(parameterList) != 1 {

								force := CalculateBondForce(parameterList[1], r, parameterList[0], atom1, atom2)

								bondEnergy += CalculateBondStretchEnergy(parameterList[1], r, parameterList[0])
								_, exist := forceMap[i+index]
								if exist {
									forceMap[i+index].x += force.x
									forceMap[i+index].y += force.y
									forceMap[i+index].z += force.z
								} else {
									forceMap[i+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
									forceMap[i+index].x = force.x
									forceMap[i+index].y = force.y
									forceMap[i+index].z = force.z
								}

								_, exist1 := forceMap[j+index]
								if exist1 {
									forceMap[j+index].x += -force.x
									forceMap[j+index].y += -force.y
									forceMap[j+index].z += -force.z
								} else {
									forceMap[j+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
									forceMap[j+index].x = -force.x
									forceMap[j+index].y = -force.y
									forceMap[j+index].z = -force.z
								}

							}

						}
					}
				}
			}

		}

		for _, angleTris := range residueParameterValue[residue.Name].angles {
			for i := 0; i < len(residue.Atoms)-2; i++ {
				atom1 := residue.Atoms[i]
				if atom1.element[0] == (*angleTris).atoms[0][0] {
					for j := i + 1; j < len(residue.Atoms)-1; j++ {
						if residue.Atoms[j].element[0] == (*angleTris).atoms[1][0] {
							atom2 := residue.Atoms[j]
							for k := j + 1; k < len(residue.Atoms); k++ {
								if residue.Atoms[k].element[0] == (*angleTris).atoms[2][0] {
									atom3 := residue.Atoms[k]
									theta := CalculateAngle(atom1, atom2, atom3)

									if math.IsNaN(theta) {
										continue
									}
									parameterList := SearchParameter(3, angleParameter, atom1, atom2, atom3)
									if len(parameterList) != 1 {

										force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

										angleEnergy += CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
										_, exist := forceMap[i+index]
										if exist {
											forceMap[i+index].x += force_i.x
											forceMap[i+index].y += force_i.y
											forceMap[i+index].z += force_i.z
										} else {
											forceMap[i+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[i+index].x = force_i.x
											forceMap[i+index].y = force_i.y
											forceMap[i+index].z = force_i.z
										}

										_, exist1 := forceMap[j+index]
										if exist1 {
											forceMap[j+index].x += force_j.x
											forceMap[j+index].y += force_j.y
											forceMap[j+index].z += force_j.z
										} else {
											forceMap[j+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[j+index].x = force_j.x
											forceMap[j+index].y = force_j.y
											forceMap[j+index].z = force_j.z
										}

										_, exist2 := forceMap[k+index]
										if exist2 {
											forceMap[k+index].x += force_k.x
											forceMap[k+index].y += force_k.y
											forceMap[k+index].z += force_k.z
										} else {
											forceMap[k+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
											forceMap[k+index].x = force_k.x
											forceMap[k+index].y = force_k.y
											forceMap[k+index].z = force_k.z
										}

									}

								}
							}

						}
					}
				}
			}
		}

		for _, dihedralValues := range residueParameterValue[residue.Name].dihedrals {
			for i := 0; i < len(residue.Atoms)-3; i++ {
				atom1 := residue.Atoms[i]
				if atom1.element[0] == (*dihedralValues).atoms[0][0] {
					for j := i + 1; j < len(residue.Atoms)-2; j++ {
						if residue.Atoms[j].element[0] == (*dihedralValues).atoms[1][0] {
							atom2 := residue.Atoms[j]
							for k := j + 1; k < len(residue.Atoms)-1; k++ {
								if residue.Atoms[k].element[0] == (*dihedralValues).atoms[2][0] {
									atom3 := residue.Atoms[k]
									for l := k + 1; l < len(residue.Atoms); l++ {
										if residue.Atoms[l].element[0] == (*dihedralValues).atoms[3][0] {
											atom4 := residue.Atoms[l]
											phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)
											if math.IsNaN(phi) {
												continue
											}
											parameterList := SearchParameter(4, dihedralParameter, atom1, atom2, atom3, atom4)
											if len(parameterList) != 1 {

												force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

												dihedralEnergy += CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])
												_, exist := forceMap[i+index]
												if exist {
													forceMap[i+index].x += force_i.x
													forceMap[i+index].y += force_i.y
													forceMap[i+index].z += force_i.z
												} else {
													forceMap[i+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[i+index].x = force_i.x
													forceMap[i+index].y = force_i.y
													forceMap[i+index].z = force_i.z
												}

												_, exist1 := forceMap[j+index]
												if exist1 {
													forceMap[j+index].x += force_j.x
													forceMap[j+index].y += force_j.y
													forceMap[j+index].z += force_j.z
												} else {
													forceMap[j+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[j+index].x = force_j.x
													forceMap[j+index].y = force_j.y
													forceMap[j+index].z = force_j.z
												}

												_, exist2 := forceMap[k+index]
												if exist2 {
													forceMap[k+index].x += force_k.x
													forceMap[k+index].y += force_k.y
													forceMap[k+index].z += force_k.z
												} else {
													forceMap[k+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[k+index].x = force_k.x
													forceMap[k+index].y = force_k.y
													forceMap[k+index].z = force_k.z
												}

												_, exist3 := forceMap[l+index]
												if exist3 {
													forceMap[l+index].x += force_l.x
													forceMap[l+index].y += force_l.y
													forceMap[l+index].z += force_l.z
												} else {
													forceMap[l+index] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
													forceMap[l+index].x = force_l.x
													forceMap[l+index].y = force_l.y
													forceMap[l+index].z = force_l.z
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
		index += len(residue.Atoms)
	}

	// additional addtion in interactions between two AA
	// additional connecting bonds
	index1 := 0
	for m, aminoA := range p.Residue {
		if m == len(p.Residue)-1 {
			break
		}

		for n := range aminoA.Atoms {
			if aminoA.Atoms[n].element == "CB" {
				atom1 := aminoA.Atoms[n]
				for t := range p.Residue[m+1].Atoms {
					if p.Residue[m+1].Atoms[t].element == "N" {
						atom2 := p.Residue[m+1].Atoms[t]
						r := Distance(atom1.position, atom2.position)
						if r == 0 {
							continue
						}
						if r > 2 {
							continue
						}
						parameterList := SearchParameter(2, bondParameter, atom1, atom2)
						if len(parameterList) != 1 {
							force := CalculateBondForce(parameterList[1], r, parameterList[0], atom1, atom2)

							bondEnergy += CalculateBondStretchEnergy(parameterList[1], r, parameterList[0])

							_, exist := forceMap[index1+n]
							if exist {
								forceMap[index1+n].x += force.x
								forceMap[index1+n].y += force.y
								forceMap[index1+n].z += force.z
							} else {
								forceMap[index1+n] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[index1+n].x = force.x
								forceMap[index1+n].y = force.y
								forceMap[index1+n].z = force.z
							}

							_, exist1 := forceMap[index1+n+len(aminoA.Atoms)]
							if exist1 {
								forceMap[index1+t+len(aminoA.Atoms)].x += -force.x
								forceMap[index1+t+len(aminoA.Atoms)].y += -force.y
								forceMap[index1+t+len(aminoA.Atoms)].z += -force.z
							} else {
								forceMap[index1+t+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[index1+t+len(aminoA.Atoms)].x = -force.x
								forceMap[index1+t+len(aminoA.Atoms)].y = -force.y
								forceMap[index1+t+len(aminoA.Atoms)].z = -force.z
							}
						}
					}
				}
			}
			break
		}
		index1 += len(aminoA.Atoms)

	}

	// additonal for angle
	index2 := 0
	for m, aminoA := range p.Residue {
		if m == len(p.Residue)-1 {
			break
		}

		for n := range aminoA.Atoms {
			if aminoA.Atoms[n].element == "CB" {
				atom1 := aminoA.Atoms[n]

				for g := range p.Residue[m+1].Atoms {
					if p.Residue[m+1].Atoms[g].element == "N" {
						atom2 := p.Residue[m+1].Atoms[g]
						for h := g + 1; h < len(p.Residue[m+1].Atoms); h++ {
							if p.Residue[m+1].Atoms[h].element[0] == 'H' {
								atom3 := p.Residue[m+1].Atoms[h]
								theta := CalculateAngle(atom1, atom2, atom3)

								parameterList := SearchParameter(3, angleParameter, atom1, atom2, atom3)
								if len(parameterList) != 1 {

									force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

									angleEnergy += CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
									_, exist := forceMap[index2+n]
									if exist {
										forceMap[index2+n].x += force_i.x
										forceMap[index2+n].y += force_i.y
										forceMap[index2+n].z += force_i.z
									} else {
										forceMap[index2+n] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index2+n].x = force_i.x
										forceMap[index2+n].y = force_i.y
										forceMap[index2+n].z = force_i.z
									}

									_, exist1 := forceMap[index2+g+len(aminoA.Atoms)]
									if exist1 {
										forceMap[index2+g+len(aminoA.Atoms)].x += force_j.x
										forceMap[index2+g+len(aminoA.Atoms)].y += force_j.y
										forceMap[index2+g+len(aminoA.Atoms)].z += force_j.z
									} else {
										forceMap[index2+g+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index2+g+len(aminoA.Atoms)].x = force_j.x
										forceMap[index2+g+len(aminoA.Atoms)].y = force_j.y
										forceMap[index2+g+len(aminoA.Atoms)].z = force_j.z
									}

									_, exist2 := forceMap[index2+h+len(aminoA.Atoms)]
									if exist2 {
										forceMap[index2+h+len(aminoA.Atoms)].x += force_k.x
										forceMap[index2+h+len(aminoA.Atoms)].y += force_k.y
										forceMap[index2+h+len(aminoA.Atoms)].z += force_k.z
									} else {
										forceMap[index2+h+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index2+h+len(aminoA.Atoms)].x = force_k.x
										forceMap[index2+h+len(aminoA.Atoms)].y = force_k.y
										forceMap[index2+h+len(aminoA.Atoms)].z = force_k.z
									}

								}
							}
						}
					}
				}
			}

		}
		index2 += len(aminoA.Atoms)

	}

	index3 := 0
	for m, aminoA := range p.Residue {
		if m == len(p.Residue)-1 {
			break
		}

		for n := range aminoA.Atoms {
			if aminoA.Atoms[n].element == "O" {
				atom1 := aminoA.Atoms[n]
				atom2 := aminoA.Atoms[n+1]
				for g := range p.Residue[m+1].Atoms {
					if p.Residue[m+1].Atoms[g].element == "N" {
						atom3 := p.Residue[m+1].Atoms[g]
						theta := CalculateAngle(atom1, atom2, atom3)
						parameterList := SearchParameter(3, angleParameter, atom1, atom2, atom3)
						if len(parameterList) != 1 {

							force_i, force_j, force_k := CalculateAngleForce(parameterList[1], theta, parameterList[0], atom1, atom2, atom3)

							angleEnergy += CalculateAnglePotentialEnergy(parameterList[1], theta, parameterList[0])
							_, exist := forceMap[index3+n]
							if exist {
								forceMap[index3+n].x += force_i.x
								forceMap[index3+n].y += force_i.y
								forceMap[index3+n].z += force_i.z
							} else {
								forceMap[index3+n] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[index3+n].x = force_i.x
								forceMap[index3+n].y = force_i.y
								forceMap[index3+n].z = force_i.z
							}

							_, exist1 := forceMap[index3+n+1]
							if exist1 {
								forceMap[index3+n+1].x += force_j.x
								forceMap[index3+n+1].y += force_j.y
								forceMap[index3+n+1].z += force_j.z
							} else {
								forceMap[index3+n+1] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[index3+n+1].x = force_j.x
								forceMap[index3+n+1].y = force_j.y
								forceMap[index3+n+1].z = force_j.z
							}

							_, exist2 := forceMap[index3+g+len(aminoA.Atoms)]
							if exist2 {
								forceMap[index3+g+len(aminoA.Atoms)].x += force_k.x
								forceMap[index3+g+len(aminoA.Atoms)].y += force_k.y
								forceMap[index3+g+len(aminoA.Atoms)].z += force_k.z
							} else {
								forceMap[index3+g+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
								forceMap[index3+g+len(aminoA.Atoms)].x = force_k.x
								forceMap[index3+g+len(aminoA.Atoms)].y = force_k.y
								forceMap[index3+g+len(aminoA.Atoms)].z = force_k.z
							}

						}

					}
				}
			}

		}
		index3 += len(aminoA.Atoms)

	}

	// additional for dihedral
	index4 := 0
	for m, aminoA := range p.Residue {
		if m == len(p.Residue)-1 {
			break
		}

		for n := range aminoA.Atoms {
			if aminoA.Atoms[n].element == "CB" {
				atom1 := aminoA.Atoms[n]
				atom2 := aminoA.Atoms[n+1]
				for g := range p.Residue[m+1].Atoms {
					if p.Residue[m+1].Atoms[g].element == "N" {
						atom3 := p.Residue[m+1].Atoms[g]
						for h := g + 1; h < len(p.Residue[m+1].Atoms); h++ {
							if p.Residue[m+1].Atoms[h].element == "H" {
								atom4 := p.Residue[m+1].Atoms[h]
								phi := CalculateDihedralAngle(atom1, atom2, atom3, atom4)
								parameterList := []float64{180.0, 6.90360, 2}
								if len(parameterList) != 1 {

									force_i, force_j, force_k, force_l := CalculateProperDihedralsForce(parameterList[1], phi, parameterList[2], parameterList[0], atom1, atom2, atom3, atom4)

									dihedralEnergy += CalculateProperDihedralAngleEnergy(parameterList[1], phi, parameterList[2], parameterList[0])
									_, exist := forceMap[index4+n]
									if exist {
										forceMap[index4+n].x += force_i.x
										forceMap[index4+n].y += force_i.y
										forceMap[index4+n].z += force_i.z
									} else {
										forceMap[index4+n] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index4+n].x = force_i.x
										forceMap[index4+n].y = force_i.y
										forceMap[index4+n].z = force_i.z
									}

									_, exist1 := forceMap[index4+n+1]
									if exist1 {
										forceMap[index4+n+1].x += force_j.x
										forceMap[index4+n+1].y += force_j.y
										forceMap[index4+n+1].z += force_j.z
									} else {
										forceMap[index4+n+1] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index4+n+1].x = force_j.x
										forceMap[index4+n+1].y = force_j.y
										forceMap[index4+n+1].z = force_j.z
									}

									_, exist2 := forceMap[index4+g+len(aminoA.Atoms)]
									if exist2 {
										forceMap[index4+g+len(aminoA.Atoms)].x += force_k.x
										forceMap[index4+g+len(aminoA.Atoms)].y += force_k.y
										forceMap[index4+g+len(aminoA.Atoms)].z += force_k.z
									} else {
										forceMap[index4+g+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index4+g+len(aminoA.Atoms)].x = force_k.x
										forceMap[index4+g+len(aminoA.Atoms)].y = force_k.y
										forceMap[index4+g+len(aminoA.Atoms)].z = force_k.z
									}

									_, exist3 := forceMap[index4+h+len(aminoA.Atoms)]
									if exist3 {
										forceMap[index4+h+len(aminoA.Atoms)].x += force_l.x
										forceMap[index4+h+len(aminoA.Atoms)].y += force_l.y
										forceMap[index4+h+len(aminoA.Atoms)].z += force_l.z
									} else {
										forceMap[index4+h+len(aminoA.Atoms)] = &TriTuple{x: 0.0, y: 0.0, z: 0.0}
										forceMap[index4+h+len(aminoA.Atoms)].x = force_l.x
										forceMap[index4+h+len(aminoA.Atoms)].y = force_l.y
										forceMap[index4+h+len(aminoA.Atoms)].z = force_l.z
									}

								}
							}
						}
					}
				}
			}

		}
		index4 += len(aminoA.Atoms)

	}
	fmt.Println("bondEnergy:", bondEnergy)
	fmt.Println("angleEnergy:", angleEnergy)
	fmt.Println("dihedralEnergy:", dihedralEnergy)
	totalEnergy := bondEnergy + angleEnergy + dihedralEnergy
	return totalEnergy, forceMap
}

func SearchParameter(value int, parameterData parameterDatabase, atoms ...*Atom) []float64 {
	for i := range parameterData.atomPair {
		sym := 0
		for j := range parameterData.atomPair[i].atomName {
			if parameterData.atomPair[i].atomName[j] == "X" {
				continue
			}
			if atoms[j].element[0] != parameterData.atomPair[i].atomName[j][0] {
				break
			}
			sym += 1
		}

		if sym == value {
			return parameterData.atomPair[i].parameter
		}

	}

	for i := range parameterData.atomPair {
		sym := 0
		for j := range parameterData.atomPair[i].atomName {
			if parameterData.atomPair[i].atomName[value-j-1] == "X" {
				continue
			}
			if atoms[j].element[0] != parameterData.atomPair[i].atomName[value-j-1][0] {
				break
			}
			sym += 1
		}

		if sym == value {
			return parameterData.atomPair[i].parameter
		}

	}

	return []float64{0.0}
}
func CalculateNetForce(a int) TriTuple {

	return TriTuple{x: 1.0, y: 1.0, z: 1.0}
}

func SteepestDescent(protein *Protein, h float64, forceMap map[int]*TriTuple) *Protein {
	for i := range protein.Residue {
		for j := range protein.Residue[i].Atoms {
			_, exist := forceMap[protein.Residue[i].Atoms[j].index+1]
			if exist {
				force := forceMap[protein.Residue[i].Atoms[j].index+1]
				magn := magnitude(*force)
				protein.Residue[i].Atoms[j].position.x = protein.Residue[i].Atoms[j].position.x + (force.x*h)/magn
				protein.Residue[i].Atoms[j].position.y = protein.Residue[i].Atoms[j].position.y + (force.y*h)/magn
				protein.Residue[i].Atoms[j].position.z = protein.Residue[i].Atoms[j].position.z + (force.z*h)/magn
			}

		}
	}

	return protein

}

func CalculateBondForce(k, r, r_0 float64, atom1, atom2 *Atom) TriTuple {
	bondLen := Distance(atom1.position, atom2.position)
	unitVector := TriTuple{
		x: (atom1.position.x - atom2.position.x) / bondLen,
		y: (atom1.position.y - atom2.position.y) / bondLen,
		z: (atom1.position.z - atom2.position.z) / bondLen,
	}

	fScale := k * (r - r_0*10)

	// unit: g/s^2/mol/A
	force := TriTuple{
		x: fScale * unitVector.x / math.Pow10(6),
		y: fScale * unitVector.y / math.Pow10(6),
		z: fScale * unitVector.z / math.Pow10(6),
	}

	return force

}

func CalculateAngleForce(k, theta, theta_0 float64, atom1, atom2, atom3 *Atom) (TriTuple, TriTuple, TriTuple) {
	der_that_cos := (-1) * (1 / math.Sin(theta/360*math.Pi))
	der_U_thate := k * (theta - theta_0) / 360 * math.Pi

	if math.IsNaN(der_that_cos) {
		return TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}
	}
	der_theta_x_12 := DerivateAnglePositionX(atom1, atom2, atom3, theta)
	der_theta_x_32 := DerivateAnglePositionX(atom3, atom2, atom1, theta)

	der_theta_y_12 := DerivateAnglePositionY(atom1, atom2, atom3, theta)
	der_theta_y_32 := DerivateAnglePositionY(atom3, atom2, atom1, theta)

	der_theta_z_12 := DerivateAnglePositionZ(atom1, atom2, atom3, theta)
	der_theta_z_32 := DerivateAnglePositionZ(atom3, atom2, atom1, theta)

	force_i := TriTuple{
		x: der_U_thate * der_that_cos * der_theta_x_12 / math.Pow10(6),
		y: der_U_thate * der_that_cos * der_theta_y_12 / math.Pow10(6),
		z: der_U_thate * der_that_cos * der_theta_z_12 / math.Pow10(6),
	}

	force_k := TriTuple{
		x: der_U_thate * der_that_cos * der_theta_x_32 / math.Pow10(6),
		y: der_U_thate * der_that_cos * der_theta_y_32 / math.Pow10(6),
		z: der_U_thate * der_that_cos * der_theta_z_32 / math.Pow10(6),
	}

	force_j := TriTuple{
		x: -force_i.x - force_k.x,
		y: -force_i.y - force_k.y,
		z: -force_i.z - force_k.z,
	}

	if math.IsNaN(force_i.x) || math.IsNaN(force_j.x) || math.IsNaN(force_k.x) {
		return TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}, TriTuple{x: 0.0, y: 0.0, z: 0.0}
	}

	return force_i, force_j, force_k

}

func DerivateAnglePositionX(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.x-atom2.position.x)/Distance(atom3.position, atom2.position) - (atom1.position.x-atom2.position.x)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

func DerivateAnglePositionY(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.y-atom2.position.y)/Distance(atom3.position, atom2.position) - (atom1.position.y-atom2.position.y)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

func DerivateAnglePositionZ(atom1, atom2, atom3 *Atom, theta float64) float64 {
	return 1 / Distance(atom1.position, atom2.position) * ((atom3.position.z-atom2.position.z)/Distance(atom3.position, atom2.position) - (atom1.position.z-atom2.position.z)/Distance(atom1.position, atom2.position)*math.Cos(theta))
}

func CalculateProperDihedralsForce(kd, phi, pn, phase float64, atom1, atom2, atom3, atom4 *Atom) (TriTuple, TriTuple, TriTuple, TriTuple) {
	der_U_phi := -0.5 * kd * pn * math.Sin((pn*phi - phase/360*math.Pi))
	der_phi_cos := -1 / math.Sin(phi)

	vector12 := CalculateVector(atom1, atom2)
	vector32 := CalculateVector(atom3, atom2)
	vector43 := CalculateVector(atom4, atom3)

	v_t := Cross(vector12, vector32)
	v_u := Cross(vector43, vector32)

	der_cos_tx, der_cos_ty, der_cos_tz := CalculateDerivate(v_t, v_u, phi)
	der_cos_ux, der_cos_uy, der_cos_uz := CalculateDerivate(v_u, v_t, phi)

	force_i := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*(-vector32.z) + der_cos_tz*vector32.y) / math.Pow10(6),
		y: der_U_phi * der_phi_cos * (der_cos_tz*(-vector32.x) + der_cos_tx*vector32.z) / math.Pow10(6),
		z: der_U_phi * der_phi_cos * (der_cos_tx*(-vector32.y) + der_cos_ty*vector32.x) / math.Pow10(6),
	}

	force_j := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*(-vector12.z+vector32.z) + der_cos_tz*(-vector32.y+vector12.y)) / math.Pow10(6),
		y: der_U_phi * der_phi_cos * (der_cos_tz*(-vector12.x+vector32.x) + der_cos_tx*(-vector32.z+vector12.z)) / math.Pow10(6),
		z: der_U_phi * der_phi_cos * (der_cos_tx*(-vector12.y+vector32.y) + der_cos_ty*(-vector32.x+vector12.x)) / math.Pow10(6),
	}

	force_k := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_ty*vector12.z - der_cos_tz*vector12.y + der_cos_uy*(vector32.z+vector43.z) - der_cos_uz*(vector32.y+vector43.y)) / math.Pow10(6),
		y: der_U_phi * der_phi_cos * (der_cos_tz*vector12.x - der_cos_tx*vector12.z + der_cos_uz*(vector32.x+vector43.x) - der_cos_ux*(vector32.z+vector43.z)) / math.Pow10(6),
		z: der_U_phi * der_phi_cos * (der_cos_tx*vector12.y - der_cos_ty*vector12.z + der_cos_ux*(vector32.y+vector43.y) - der_cos_uy*(vector32.x+vector43.x)) / math.Pow10(6),
	}

	force_l := TriTuple{
		x: der_U_phi * der_phi_cos * (der_cos_uy*(-vector32.z) + der_cos_uz*vector32.y) / math.Pow10(6),
		y: der_U_phi * der_phi_cos * (der_cos_uz*(-vector32.x) + der_cos_ux*vector32.z) / math.Pow10(6),
		z: der_U_phi * der_phi_cos * (der_cos_ux*(-vector32.y) + der_cos_uy*vector32.x) / math.Pow10(6),
	}

	return force_i, force_j, force_k, force_l

}

func CalculateDerivate(v1, v2 TriTuple, phi float64) (float64, float64, float64) {
	return (1 / magnitude(v1)) * (v2.x/magnitude(v2) - v1.x/magnitude(v1)*math.Cos(phi)),
		(1 / magnitude(v1)) * (v2.y/magnitude(v2) - v1.y/magnitude(v1)*math.Cos(phi)),
		(1 / magnitude(v1)) * (v2.z/magnitude(v2) - v1.z/magnitude(v1)*math.Cos(phi))
}

func Cross(v1, v2 TriTuple) TriTuple {
	return TriTuple{
		x: v1.x * v2.x,
		y: v1.y * v2.y,
		z: v1.z * v2.z,
	}
}

func CopyProtein(currentProtein *Protein) *Protein {
	var newProtein Protein
	newProtein.Name = currentProtein.Name

	newProtein.Residue = make([]*Residue, len(currentProtein.Residue))
	for i := range currentProtein.Residue {
		newProtein.Residue[i] = CopyResidue(currentProtein.Residue[i])
	}

	return &newProtein
}

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

func CopyTriTuple(tri TriTuple) TriTuple {
	var newTri TriTuple
	newTri.x = tri.x
	newTri.y = tri.y
	newTri.z = tri.z

	return newTri
}

func CalculateVector(atom1, atom2 *Atom) TriTuple {
	var vector TriTuple
	vector.x = (atom2.position.x - atom1.position.x)
	vector.y = (atom2.position.y - atom1.position.y)
	vector.z = (atom2.position.z - atom1.position.z)

	return vector
}
