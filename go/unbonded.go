package main

import (
	"math"
)

// q1: the charge of ion 1
// q2: the charge of ion 2
// epsilon: the vacuum dielectric permittivity
// r: distance between q1 and q2
func CalculateElectricPotentialEnergy(q1, q2, epsilon, r float64) float64 {
	return (q1 * q2) / (4 * math.Pi * epsilon * r)
}

// A: coefficient 1
// B: coefficient 2
// r: distance between atom 1 and atom 2
func CalculateLJPotentialEnergy(A, B, r float64) float64 {
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
