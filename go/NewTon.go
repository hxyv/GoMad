// This file contains all the functions related to NewTon equation
package main

import (
	"fmt"
	"math"
)

// SimulateMD
// Input: an initial Protein object, related parameter dataset, and a float time.
// Output: a slice of Protein resulting from MD simulation, where the time interval between generations is specified by time.
func SimulateMD(initialProtein Protein, time float64, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) []Protein {
	timePoints := make([]Protein, 0)
	cerition := 10000000.0
	timePoints = append(timePoints, initialProtein)
	totalTime := 0.0
	iteration := 100 // 100

	// range over each iteration
	for i := 0; i < iteration; i++ {
		newProtein, _ := UpdateProtein(timePoints[len(timePoints)-1], time, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
		timePoints = append(timePoints, newProtein)

		// cerition of limited time
		totalTime += time
		if totalTime > cerition {
			break
		}
		fmt.Println("Distance between input and new protein at 962 N:", Distance(initialProtein.Residue[65].Atoms[0].position, newProtein.Residue[65].Atoms[0].position))
	}

	return timePoints
}

// UpdateProtein
// Input: a Protein object, related parameter dataset and a float time.
// Output: a Protein object resulting from a single step according to bonded and unbonded interaction, using a time interval specified by time.
func UpdateProtein(currentProtein Protein, time float64, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (Protein, float64) {
	newProtein := CopyProtein(&currentProtein)

	energy, forceMap := CombineEnergyAndForce(newProtein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)

	forceIndex := 1

	// range and update every atom in the protein
	for i, b := range newProtein.Residue {
		for j, a := range b.Atoms {
			_, exist := forceMap[forceIndex]

			if exist {

				oldAcceleration, oldVelocity := a.accelerated, a.velocity // OK :)
				newProtein.Residue[i].Atoms[j].accelerated = UpdateAcceleration(forceMap[forceIndex], a)
				newProtein.Residue[i].Atoms[j].velocity = UpdateVelocity(a, oldAcceleration, time)
				newProtein.Residue[i].Atoms[j].position = UpdatePosition(a, oldAcceleration, oldVelocity, time)

			}
			forceIndex += 1
		}

	}
	return *newProtein, energy
}

// UpdateVelocity
// Input: a pointer of Atom object a, previous acceleration and velocity as TriTuple objects, and a float time.
// Output: Updated velocity vector of a as an TriTuple according to physics nerd velocity update equations.
func UpdateVelocity(a *Atom, oldAcceleration TriTuple, time float64) TriTuple {
	var currentVelocity TriTuple // starts at (0, 0)

	// apply equations. Note that b's acceleration has already been updated :)
	currentVelocity.x = a.velocity.x + 0.5*(a.accelerated.x+oldAcceleration.x)*time
	currentVelocity.y = a.velocity.y + 0.5*(a.accelerated.y+oldAcceleration.y)*time
	currentVelocity.z = a.velocity.z + 0.5*(a.accelerated.z+oldAcceleration.z)*time

	return currentVelocity
}

// UpdatePosition
// Input: a pointer of Atom object a, previous acceleration and velocity as OTriTuple objects, and a float time.
// Output: Updated position of a as an TriTuple according to physics nerd update equations.
func UpdatePosition(a *Atom, oldAcceleration, oldVelocity TriTuple, time float64) TriTuple {
	var pos TriTuple

	pos.x = a.position.x + oldVelocity.x*time + 0.5*oldAcceleration.x*time*time
	pos.y = a.position.y + oldVelocity.y*time + 0.5*oldAcceleration.y*time*time
	pos.z = a.position.z + oldVelocity.z*time + 0.5*oldAcceleration.z*time*time

	return pos
}

// UpdateAcceleration
// Input: a pointer of Atom object and a force in TriTuple.
// Output: The acceleration of a in the next generation after computing the net force acting on a from all the interactions.
func UpdateAcceleration(force *TriTuple, a *Atom) TriTuple {
	var accel TriTuple

	// F = m * a or a = F/m according to Newton law
	accel.x = force.x / a.mass
	accel.y = force.y / a.mass
	accel.z = force.z / a.mass

	return accel
}

// calculate RMSD in each timepoint
// Input: a series of Protein
// Output: a list of RMSD
func CalculateRMSD(timePoints []Protein) []float64 {
	var RMSDValue []float64

	for i := 1; i < len(timePoints); i++ {
		length := 0.0
		value := 0.0
		for j := range timePoints[i].Residue {
			for k := range timePoints[i].Residue[j].Atoms {
				dis := Distance(timePoints[0].Residue[j].Atoms[k].position, timePoints[i].Residue[j].Atoms[k].position)

				value += dis * dis
				length += 1.0
			}
		}
		RMSDValue = append(RMSDValue, math.Sqrt(value/length))

	}

	return RMSDValue
}
