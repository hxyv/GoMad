package main

import (
	"fmt"
	"math"
)

// SimulateGravity
// Input: an initial Universe object, a number of generations, and a float time.
// Output: a slice of numGens + 1 Universes resulting from simulating gravity over numGens generations, where the time interval between generations is specified by time.
func SimulateMD(initialProtein Protein, time float64, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) []Protein {
	timePoints := make([]Protein, 0)
	cerition := 100000000000.0
	timePoints = append(timePoints, initialProtein)
	totalTime := 0.0
	iteration := 100 // 100
	CheckPosition(timePoints[0])
	fmt.Println("after first check")
	for i := 0; i < iteration; i++ {
		newProtein, _ := UpdateProtein(timePoints[len(timePoints)-1], time, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)
		timePoints = append(timePoints, newProtein)
		CheckPosition(timePoints[len(timePoints)-1])
		totalTime += time
		if totalTime > cerition {
			break
		}
		fmt.Println("Distance between input and new protein at 962 N:", Distance(initialProtein.Residue[65].Atoms[0].position, newProtein.Residue[65].Atoms[0].position))
	}

	return timePoints
}

// UpdateUniverse
// Input: a Universe object and a float time.
// Output: a Universe object resulting from a single step according to the gravity simulation, using a time interval specified by time.
func UpdateProtein(currentProtein Protein, time float64, residueParameterBondValue, residueParameterOtherValue map[string]residueParameter, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter parameterDatabase) (Protein, float64) {
	newProtein := CopyProtein(&currentProtein)

	energy, forceMap := CombineEnergyAndForce(newProtein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondParameter, pairtypesParameter)

	forceIndex := 1

	// range and update every body in universe
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

	//AddSimpleBondConstraints(newProtein, bondParameter, residueParameterValue)

	return *newProtein, energy
}

// UpdateVelocity
// Input: Body object b, previous acceleration and velocity as OrderedPair objects, and a float time.
// Output: Updated velocity vector of b as an OrderedPair according to physics nerd velocity update equations.
func UpdateVelocity(a *Atom, oldAcceleration TriTuple, time float64) TriTuple {
	var currentVelocity TriTuple // starts at (0, 0)

	// apply equations. Note that b's acceleration has already been updated :)
	currentVelocity.x = a.velocity.x + 0.5*(a.accelerated.x+oldAcceleration.x)*time
	currentVelocity.y = a.velocity.y + 0.5*(a.accelerated.y+oldAcceleration.y)*time
	currentVelocity.z = a.velocity.z + 0.5*(a.accelerated.z+oldAcceleration.z)*time

	return currentVelocity
}

// UpdatePosition
// Input: Body object b, previous acceleration and velocity as OrderedPair objects, and a float time.
// Output: Updated position of b as an OrderedPair according to physics nerd update equations.
func UpdatePosition(a *Atom, oldAcceleration, oldVelocity TriTuple, time float64) TriTuple {
	var pos TriTuple

	pos.x = a.position.x + oldVelocity.x*time + 0.5*oldAcceleration.x*time*time
	pos.y = a.position.y + oldVelocity.y*time + 0.5*oldAcceleration.y*time*time
	pos.z = a.position.z + oldVelocity.z*time + 0.5*oldAcceleration.z*time*time

	return pos
}

// UpdateAcceleration
// Input: currentUniverse Universe object and a Body object b.
// Output: The acceleration of b in the next generation after computing the net force of gravity acting on b over all bodies in currentUniverse.
func UpdateAcceleration(force *TriTuple, a *Atom) TriTuple {
	var accel TriTuple

	// F = m * a or a = F/m according to Newton law
	accel.x = force.x / a.mass
	accel.y = force.y / a.mass
	accel.z = force.z / a.mass

	return accel
}

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
		//fmt.Println(RMSDValue[len(RMSDValue)-1])

	}

	return RMSDValue
}

func CheckPosition(p Protein) {
	fmt.Println("Check the nan position index")
	for _, residue := range p.Residue {
		for _, atom := range residue.Atoms {
			if math.IsNaN(atom.position.x) || math.IsNaN(atom.position.y) || math.IsNaN(atom.position.z) {
				fmt.Println(atom.index)
			}

		}
	}
}
