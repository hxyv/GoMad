package main

import (
	"math"
)

func CalculateBondStretchEnergy(k, r, r_0 float64) float64 {
	return 0.5 * k * (r - r_0) * (r - r_0)
}

func CalculateAnglePotentialEnergy(k, theta, theta_0 float64) float64 {
	return 0.5 * k * (theta - theta_0) * (theta - theta_0)
}

func CalculateProperDihedralAngleEnergy(kd, phi, pn, phase float64) float64 {
	return kd * (1 + math.Cos(pn*phi-phase))
}

func CalculateAngle(atom1, atom2, atom3 *Atom) float64 {
	var vector1 TriTuple
	var vector2 TriTuple
	vector1.x = atom2.position.x - atom1.position.x
	vector1.y = atom2.position.y - atom1.position.y
	vector1.z = atom2.position.z - atom1.position.z

	vector2.x = atom3.position.x - atom1.position.x
	vector2.y = atom3.position.y - atom1.position.y
	vector2.z = atom3.position.z - atom1.position.z

	upperValue := vector1.dot(vector2)
	lowerValue := Distance(atom1.position, atom2.position) * Distance(atom2.position, atom3.position)

	value := upperValue / lowerValue

	return math.Acos(value)
}

func Distance(p1, p2 TriTuple) float64 {
	deltaX := p1.x - p2.x
	deltaY := p1.y - p2.y
	deltaZ := p1.z - p2.z
	return math.Sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ)

}

func CalculateDihedralAngle(atom1, atom2, atom3, atom4 *Atom) float64 {
	var vector1 TriTuple
	var vector2 TriTuple
	var vector3 TriTuple

	vector1.x = atom2.position.x - atom1.position.x
	vector1.y = atom2.position.y - atom1.position.y

	vector2.x = atom3.position.x - atom2.position.x
	vector2.y = atom3.position.y - atom2.position.y

	vector3.x = atom4.position.x - atom4.position.x
	vector3.y = atom4.position.y - atom4.position.y

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

func CalculateTotoalBondStretchEnergy(k, r, r_0 float64) float64 {
	return 0.5 * k * (r - r_0) * (r - r_0)
}

func PerformEnergyMinimization(currentProtein *Protein) *Protein {

	iteration := 100
	// set maximum displacement
	h := 0.01

	for i := 0; i < iteration; i++ {
		// calculate total Energy of original protein
		initialEnergy := CalculateTotalEnergy(currentProtein)

		tempProtein := CopyProtein(currentProtein)

		// perform SteepestDescent, update positions in protein
		SteepestDescent(tempProtein, h)

		// calculate total Energy of updated protein
		updatedEnergy := CalculateTotalEnergy(tempProtein)

		// if Total energy decrease, accept the changes of positions and increase maximum displacement h
		// Otherwise, reject the changes in positions and decrease maximum displacement h
		if updatedEnergy < initialEnergy {
			currentProtein = tempProtein
			h *= 1.2
		} else {
			h *= 0.2 * h
		}

	}

	return currentProtein

}

func CalculateTotalEnergy()

func SteepestDescent()

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

	return &newAtom
}

func CopyTriTuple(tri TriTuple) TriTuple {
	var newTri TriTuple
	newTri.x = tri.x
	newTri.y = tri.y
	newTri.z = tri.z

	return newTri
}
