// This file contains all the functions tests
package main

import (
	"bufio"
	"fmt"
	"io/fs"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

// //////////
// datatype area
// //////////

// SeparationForceTest is a struct holds the information for a test that take three float64 and two Boid
type UpdateAccelerationTest = struct {
	a           Atom
	force       TriTuple
	resultaccel TriTuple
}

type UpdatePositionTest = struct {
	a               Atom
	oldAcceleration TriTuple
	oldVelocity     TriTuple
	time            float64
	result          TriTuple
}

type UpdateVelocityTest = struct {
	a               Atom
	oldAcceleration TriTuple
	time            float64
	result          TriTuple
}

type CalculateRMSDTest = struct {
	timePoints []Protein
	result     []float64
}

type UpdateProteinTest = struct {
	currentProtein Protein
	time           float64
}

type CalculateElectricPotentialEnergyTest struct {
	a1, a2 Atom
	r      float64
	result float64
}

type CalculateLJPotentialEnergyTest struct {
	A, B, r float64
	result  float64
}

type CalculateElectricForceTest struct {
	a1, a2 Atom
	r      float64
	result TriTuple
}

type CalculateLJForceTest struct {
	A, B, r float64
	a1, a2  Atom
	result  TriTuple
}

// //////////
// Test area
// //////////

// Test CalculateElectricPotentialEnergy
func TestCalculateElectricPotentialEnergy(t *testing.T) {
	tests := ReadCalculateElectricPotentialEnergyTests("Tests/CalculateElectricPotentialEnergy/")

	for _, test := range tests {
		result := CalculateElectricPotentialEnergy(&test.a1, &test.a2, test.r)
		if result != test.result {
			t.Errorf("CalculateElectricPotentialEnergy(%v, %v, %v) = %v; expected %v", test.a1, test.a2, test.r, result, test.result)
		}
	}
}

// Test CalculateLJPotentialEnergy
func TestCalculateLJPotentialEnergy(t *testing.T) {
	tests := ReadCalculateLJPotentialEnergyTest("Tests/CalculateLJPotentialEnergy/")

	for _, test := range tests {
		result := CalculateLJPotentialEnergy(test.B, test.A, test.r)
		if result != test.result {
			t.Errorf("CalculateLJPotentialEnergy(%v, %v, %v) = %v; expected %v", test.A, test.B, test.r, result, test.result)
		}
	}
}

// Test CalculateElectricForce
func TestCalculateElectricForce(t *testing.T) {
	tests := ReadCalculateElectricForceTests("Tests/CalculateElectricForce/")

	for _, test := range tests {
		result := CalculateElectricForce(&test.a1, &test.a2, test.r)
		if result != test.result {
			t.Errorf("CalculateElectricForce(%v, %v, %v) = %v; expected %v", test.a1, test.a2, test.r, result, test.result)
		}
	}
}

// Test CalculateLJForce
func TestCalculateLJForce(t *testing.T) {
	tests := ReadCalculateLJForceTest("Tests/CalculateLJForce/")

	for _, test := range tests {
		result := CalculateLJForce(&test.a1, &test.a2, test.B, test.A, test.r)
		if result != test.result {
			t.Errorf("CalculateLJForce(%v, %v, %v, %v, %v) = %v; expected %v", test.A, test.B, test.r, test.A, test.B, result, test.result)
		}
	}
}

func TestCalculateTotalUnbondedEnergyForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateTotalUnbondedEnergyForce" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateTotalUnbondedEnergyForce" + "/output")

	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/CalculateTotalUnbondedEnergyForce/" + "input/" + inputFile.Name())

		// function
		result1, _ := CalculateTotalUnbondedEnergyForce(&protein, nonbondedParameter)

		// real output
		realResult1, _, _ := ReadFloatMapIntTriTuple("Tests/CalculateTotalUnbondedEnergyForce/" + "output/" + outputFiles[i].Name())

		// compare
		if realResult1 != result1 {
			t.Errorf("Energy in CalculateTotalUnbondedEnergyForce() = %v, want %v in %v", result1, realResult1, outputFiles[i])

		}

	}
}

// Test CalculateRMSD
func TestCalculateRMSDTest(t *testing.T) {
	// Read in all tests from the Tests/Distance directory and run them
	tests := ReadCalculateRMSDTestTests("Tests/CalculateRMSD/")
	for _, test := range tests {
		// Run the test
		ourAnswer := CalculateRMSD(test.timePoints)
		// Check the result
		for i := range ourAnswer {
			if ourAnswer[i] != test.result[i] {
				t.Errorf("CalculateRMSD(%v) = %v, want %v", i, ourAnswer[i], test.result[i])
			}
		}
	}
}

func TestUpdateProteinTest(t *testing.T) {
	inputFiles := ReadDirectory("Tests/UpdateProtein" + "/input")
	outputFiles := ReadDirectory("Tests/UpdateProtein" + "/output")

	time := 1.0
	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/UpdateProtein/" + "input/" + inputFile.Name())

		// function
		result, _ := UpdateProtein(protein, time, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)

		// real output
		realResult, _ := ReadOneProtein("Tests/UpdateProtein" + "/output/" + outputFiles[i].Name())

		// compare
		if realResult.Name != result.Name {
			t.Errorf("UpdateProtein() = %v, want %v", result, realResult)
		}

		// compare residues
		for i := range realResult.Residue {
			if len(realResult.Residue[i].Atoms) == len(result.Residue[i].Atoms) {
				if realResult.Residue[i].ChainID != result.Residue[i].ChainID || realResult.Residue[i].Name != result.Residue[i].Name || realResult.Residue[i].ID != result.Residue[i].ID {
					t.Errorf("UpdateProtein() = %v, want %v", result, realResult)
				}
				for j := range realResult.Residue[i].Atoms {
					if realResult.Residue[i].Atoms[j].position.x != result.Residue[i].Atoms[j].position.x {
						t.Errorf("UpdateProtein() = %v, want %v, in %v.", realResult.Residue[i].Atoms[j], result.Residue[i].Atoms[j], outputFiles[i].Name())
					}
				}
			} else {
				t.Errorf("Number mismatch")
			}
		}

	}
}

func ReadCalculateRMSDTestTests(directory string) []CalculateRMSDTest {

	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CalculateRMSDTest, numFiles)
	for i, inputFile := range inputFiles {
		tests[i].timePoints, _ = ReadProteins(directory + "input/" + inputFile.Name())
	}
	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		out, _ := readFileline(directory + "output/" + outputFile.Name())
		tests[i].result = convertStringToFloatSlice(out[0])
	}

	return tests
}

func TestUpdateAcceleration(t *testing.T) {
	// Read in all tests from the Tests/UpdateAcceleration directory and run them
	tests := ReadUpdateAccelerationTests("Tests/UpdateAcceleration/")
	for _, test := range tests {
		// Run the test
		ourAnswer := UpdateAcceleration(&test.force, &test.a)
		// Check the result
		if ourAnswer != test.resultaccel {
			t.Errorf("UpdateAcceleration(%v, %v) = %v, want %v", test.force, test.a, ourAnswer, test.resultaccel)
		}
	}
}

// Test UpdatePosition
func TestUpdatePosition(t *testing.T) {
	// Read in all tests from the Tests/UpdatePosition directory and run them
	tests := ReadUpdatePositionTests("Tests/UpdatePosition/")
	for _, test := range tests {
		// Run the test
		ourAnswer := UpdatePosition(&test.a, test.oldAcceleration, test.oldVelocity, test.time)
		// Check the result
		if ourAnswer != test.result {
			t.Errorf("UpdateAcceleration() = %v, want %v", ourAnswer, test.result)
		}
	}
}

// Test UpdateVelocity
func TestUpdateVelocity(t *testing.T) {
	// Read in all tests from the Tests/UpdateVelocity directory and run them
	tests := ReadUpdatePositionTests("Tests/UpdateVelocity/")
	for _, test := range tests {
		// Run the test
		ourAnswer := UpdateVelocity(&test.a, test.oldAcceleration, test.time)
		// Check the result
		if ourAnswer != test.result {
			t.Errorf("UpdateAcceleration() = %v, want %v", ourAnswer, test.result)
		}
	}
}

// test Distance function
func TestDistance(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/Distance" + "/input")
	outputFiles := ReadDirectory("Tests/Distance" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/Distance/" + "input/" + inputFile.Name())
		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[1])[0]
		atom1.position.y = convertStringToFloatSlice(pair[1])[1]
		atom1.position.z = convertStringToFloatSlice(pair[1])[2]
		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[0])[0]
		atom2.position.y = convertStringToFloatSlice(pair[0])[1]
		atom2.position.z = convertStringToFloatSlice(pair[0])[2]
		// function
		result := ConvertFloat(Distance(atom1.position, atom2.position))

		// read output
		out, _ := readFileline("Tests/Distance" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult, _ = strconv.ParseFloat(out[0], 64)
		realResult = ConvertFloat(realResult)
		// compare
		if realResult != result {
			t.Errorf("Distance() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateVector
func TestCalculateVector(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CalculateVector" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateVector" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateVector/" + "input/" + inputFile.Name())
		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[1])[0]
		atom1.position.y = convertStringToFloatSlice(pair[1])[1]
		atom1.position.z = convertStringToFloatSlice(pair[1])[2]
		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[0])[0]
		atom2.position.y = convertStringToFloatSlice(pair[0])[1]
		atom2.position.z = convertStringToFloatSlice(pair[0])[2]

		// function
		result := CalculateVector(&atom1, &atom2)

		// read output
		out, _ := readFileline("Tests/CalculateVector" + "/output/" + outputFiles[i].Name())
		var realResult TriTuple
		realResult.x = convertStringToFloatSlice(out[0])[0]
		realResult.y = convertStringToFloatSlice(out[0])[1]
		realResult.z = convertStringToFloatSlice(out[0])[2]

		// compare
		if realResult != result {
			t.Errorf("CalculateVector() = %v, want %v", result, realResult)
		}

	}
}

// Test CopyTriTuple
func TestCopyTriTuple(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CopyTriTuple" + "/input")
	outputFiles := ReadDirectory("Tests/CopyTriTuple" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CopyTriTuple/" + "input/" + inputFile.Name())
		var temp TriTuple
		temp.x = convertStringToFloatSlice(pair[0])[0]
		temp.y = convertStringToFloatSlice(pair[0])[1]
		temp.z = convertStringToFloatSlice(pair[0])[2]

		// function
		result := CopyTriTuple(temp)

		// read output
		out, _ := readFileline("Tests/CopyTriTuple" + "/output/" + outputFiles[i].Name())
		var realResult TriTuple
		realResult.x = convertStringToFloatSlice(out[0])[0]
		realResult.y = convertStringToFloatSlice(out[0])[1]
		realResult.z = convertStringToFloatSlice(out[0])[2]

		// compare
		if realResult != result {
			t.Errorf("TriTuple() = %v, want %v", result, realResult)
		}

	}
}

// Test CopyAtom
func TestCopyAtom(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CopyAtom" + "/input")
	outputFiles := ReadDirectory("Tests/CopyAtom" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		atomlist, _ := ReadAtoms("Tests/CopyAtom/" + "input/" + inputFile.Name())
		atom := atomlist[0]

		// function
		result := CopyAtom(&atom)

		// read output
		atomlist1, _ := ReadAtoms("Tests/CopyAtom/" + "output/" + outputFiles[i].Name())
		realResult := atomlist1[0]

		// compare
		if realResult != *result {
			t.Errorf("CopyAtom() = %v, want %v", result, realResult)
		}

	}
}

// Test CopyResidue
func TestCopyResidue(t *testing.T) {
	inputFiles := ReadDirectory("Tests/CopyResidue" + "/input")
	outputFiles := ReadDirectory("Tests/CopyResidue" + "/output")

	for i, inputFile := range inputFiles {
		// read input

		residueList, _ := ReadResidues("Tests/CopyResidue/" + "input/" + inputFile.Name())
		residue := residueList[0]
		// function
		result := CopyResidue(&residue)

		// read output
		residueList1, _ := ReadResidues("Tests/CopyResidue/" + "output/" + outputFiles[i].Name())
		realResult := residueList1[0]

		// compare
		if realResult.ChainID != result.ChainID || realResult.Name != result.Name || realResult.ID != result.ID {
			t.Errorf("CopyResidue() = %v, want %v", result, realResult)
		}

		// compare atoms
		if len(realResult.Atoms) == len(result.Atoms) {
			for i := range realResult.Atoms {
				if (*realResult.Atoms[i]) != (*result.Atoms[i]) {
					t.Errorf("CopyResidue() = %v, want %v, in %v.", result, realResult, outputFiles[i].Name())
				}
			}
		} else {
			t.Errorf("Number mismatch")
		}

	}
}

// Test CopyProtein
func TestCopyProtein(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CopyProtein" + "/input")
	outputFiles := ReadDirectory("Tests/CopyProtein" + "/output")

	for i, inputFile := range inputFiles {
		// read input

		protein, _ := ReadOneProtein("Tests/CopyProtein/" + "input/" + inputFile.Name())

		// function
		result := CopyProtein(&protein)

		// read output

		realResult, _ := ReadOneProtein("Tests/CopyProtein/" + "output/" + outputFiles[i].Name())

		// compare
		if realResult.Name != result.Name {
			t.Errorf("CopyProtein() = %v, want %v", result, realResult)
		}

		// compare residues
		for i := range protein.Residue {
			if len(protein.Residue[i].Atoms) == len(result.Residue[i].Atoms) {
				if protein.Residue[i].ChainID != result.Residue[i].ChainID || protein.Residue[i].Name != result.Residue[i].Name || protein.Residue[i].ID != result.Residue[i].ID {
					t.Errorf("CopyProtein() = %v, want %v", result, realResult)
				}
				for j := range protein.Residue[i].Atoms {
					if (*protein.Residue[i].Atoms[j]) != (*result.Residue[i].Atoms[j]) {
						t.Errorf("CopyProtein() = %v, want %v, in %v.", result, realResult, outputFiles[i].Name())
					}
				}
			} else {
				t.Errorf("Number mismatch")
			}
		}

	}
}

// Test Cross
func TestCross(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/Cross" + "/input")
	outputFiles := ReadDirectory("Tests/Cross" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/Cross/" + "input/" + inputFile.Name())
		var v1 TriTuple
		v1.x = convertStringToFloatSlice(pair[0])[0]
		v1.y = convertStringToFloatSlice(pair[0])[1]
		v1.z = convertStringToFloatSlice(pair[0])[2]

		var v2 TriTuple
		v2.x = convertStringToFloatSlice(pair[1])[0]
		v2.y = convertStringToFloatSlice(pair[1])[1]
		v2.z = convertStringToFloatSlice(pair[1])[2]
		// function
		result := Cross(v1, v2)

		// read output
		out, _ := readFileline("Tests/Cross" + "/output/" + outputFiles[i].Name())
		var realResult TriTuple
		realResult.x = convertStringToFloatSlice(out[0])[0]
		realResult.y = convertStringToFloatSlice(out[0])[1]
		realResult.z = convertStringToFloatSlice(out[0])[2]

		// compare
		if realResult != result {
			t.Errorf("Cross() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateDerivate
func TestCalculateDerivate(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CalculateDerivate" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateDerivate" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateDerivate/" + "input/" + inputFile.Name())
		var v1 TriTuple
		v1.x = convertStringToFloatSlice(pair[0])[0]
		v1.y = convertStringToFloatSlice(pair[0])[1]
		v1.z = convertStringToFloatSlice(pair[0])[2]

		var v2 TriTuple
		v2.x = convertStringToFloatSlice(pair[1])[0]
		v2.y = convertStringToFloatSlice(pair[1])[1]
		v2.z = convertStringToFloatSlice(pair[1])[2]

		phi := convertStringToFloatSlice(pair[0])[0]

		// function
		der1, der2, der3 := CalculateDerivate(v1, v2, phi)
		der1 = math.Round(der1*math.Pow10(4)) / math.Pow10(4)
		der2 = math.Round(der2*math.Pow10(4)) / math.Pow10(4)
		der3 = math.Round(der3*math.Pow10(4)) / math.Pow10(4)

		// read output
		out, _ := readFileline("Tests/CalculateDerivate" + "/output/" + outputFiles[i].Name())
		result1 := convertStringToFloatSlice(out[0])[0]
		result2 := convertStringToFloatSlice(out[0])[1]
		result3 := convertStringToFloatSlice(out[0])[2]

		// compare
		if der1 != result1 || der2 != result2 || der3 != result3 {
			t.Errorf("CalculateDerivate() = %v, %v, %v, want (%v, %v, %v)", der1, der2, der3, result1, result2, result3)
		}

	}
}

// Test Magnitude
func TestMagnitude(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/magnitude" + "/input")
	outputFiles := ReadDirectory("Tests/magnitude" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/magnitude/" + "input/" + inputFile.Name())
		var v1 TriTuple
		v1.x = convertStringToFloatSlice(pair[0])[0]
		v1.y = convertStringToFloatSlice(pair[0])[1]
		v1.z = convertStringToFloatSlice(pair[0])[2]

		// function
		result := magnitude(v1)

		// read output
		out, _ := readFileline("Tests/magnitude" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("magnitude() = %v, want %v", result, realResult)
		}

	}
}

// Test Dot
func TestDot(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/Dot" + "/input")
	outputFiles := ReadDirectory("Tests/Dot" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/Dot/" + "input/" + inputFile.Name())
		var v1 TriTuple
		v1.x = convertStringToFloatSlice(pair[0])[0]
		v1.y = convertStringToFloatSlice(pair[0])[1]
		v1.z = convertStringToFloatSlice(pair[0])[2]

		var v2 TriTuple
		v2.x = convertStringToFloatSlice(pair[1])[0]
		v2.y = convertStringToFloatSlice(pair[1])[1]
		v2.z = convertStringToFloatSlice(pair[1])[2]

		// function
		result := v1.dot(v2)

		// read output
		out, _ := readFileline("Tests/Dot" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("Dot() = %v, want %v", result, realResult)
		}

	}
}

// Test BuildNormalVector
func TestBuildNormalVector(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/BuildNormalVector" + "/input")
	outputFiles := ReadDirectory("Tests/BuildNormalVector" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/BuildNormalVector/" + "input/" + inputFile.Name())
		var v1 TriTuple
		v1.x = convertStringToFloatSlice(pair[0])[0]
		v1.y = convertStringToFloatSlice(pair[0])[1]
		v1.z = convertStringToFloatSlice(pair[0])[2]

		var v2 TriTuple
		v2.x = convertStringToFloatSlice(pair[1])[0]
		v2.y = convertStringToFloatSlice(pair[1])[1]
		v2.z = convertStringToFloatSlice(pair[1])[2]

		// function
		result := BuildNormalVector(v1, v2)
		// read output
		out, _ := readFileline("Tests/BuildNormalVector" + "/output/" + outputFiles[i].Name())
		var realResult TriTuple
		realResult.x = convertStringToFloatSlice(out[0])[0]
		realResult.y = convertStringToFloatSlice(out[0])[1]
		realResult.z = convertStringToFloatSlice(out[0])[2]

		// compare
		if realResult != result {
			t.Errorf("BuildNormalVector() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateAngle
func TestCalculateAngle(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CalculateAngle" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateAngle" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateAngle/" + "input/" + inputFile.Name())
		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[0])[0]
		atom1.position.y = convertStringToFloatSlice(pair[0])[1]
		atom1.position.z = convertStringToFloatSlice(pair[0])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[1])[0]
		atom2.position.y = convertStringToFloatSlice(pair[1])[1]
		atom2.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[2])[0]
		atom3.position.y = convertStringToFloatSlice(pair[2])[1]
		atom3.position.z = convertStringToFloatSlice(pair[2])[2]

		// function
		result := CalculateAngle(&atom1, &atom2, &atom3)
		// read output
		out, _ := readFileline("Tests/CalculateAngle" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != math.Round(result) {
			t.Errorf("CalculateAngle() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateDihedralAngle
func TestCalculateDihedralAngle(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CalculateDihedralAngle" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateDihedralAngle" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateDihedralAngle/" + "input/" + inputFile.Name())
		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[0])[0]
		atom1.position.y = convertStringToFloatSlice(pair[0])[1]
		atom1.position.z = convertStringToFloatSlice(pair[0])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[1])[0]
		atom2.position.y = convertStringToFloatSlice(pair[1])[1]
		atom2.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[2])[0]
		atom3.position.y = convertStringToFloatSlice(pair[2])[1]
		atom3.position.z = convertStringToFloatSlice(pair[2])[2]

		var atom4 Atom
		atom4.position.x = convertStringToFloatSlice(pair[3])[0]
		atom4.position.y = convertStringToFloatSlice(pair[3])[1]
		atom4.position.z = convertStringToFloatSlice(pair[3])[2]

		// function
		result := CalculateDihedralAngle(&atom1, &atom2, &atom3, &atom4)
		// read output
		out, _ := readFileline("Tests/CalculateDihedralAngle" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != math.Round(result) {
			t.Errorf("CalculateDihedralAngle() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateBondStretchEnergy
func TestCalculateBondStretchEnergy(t *testing.T) {
	// input and output files
	inputFiles := ReadDirectory("Tests/CalculateBondStretchEnergy" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateBondStretchEnergy" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateBondStretchEnergy/" + "input/" + inputFile.Name())

		k := convertStringToFloatSlice(pair[0])[0]
		r := convertStringToFloatSlice(pair[0])[1]
		r_0 := convertStringToFloatSlice(pair[0])[2]

		// function
		result := CalculateBondStretchEnergy(k, r, r_0)

		// read output
		out, _ := readFileline("Tests/CalculateBondStretchEnergy" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != math.Round(result) {
			t.Errorf("CalculateBondStretchEnergy() = %v, want %v", result, realResult)
		}

	}
}

// Test AnglePotentialEnergy
func TestCalculateAnglePotentialEnergy(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateAnglePotentialEnergy" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateAnglePotentialEnergy" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateAnglePotentialEnergy/" + "input/" + inputFile.Name())

		k := convertStringToFloatSlice(pair[0])[0]
		theta := convertStringToFloatSlice(pair[0])[1]
		theta_0 := convertStringToFloatSlice(pair[0])[2]

		// function
		result := CalculateAnglePotentialEnergy(k, theta, theta_0)

		// read output
		out, _ := readFileline("Tests/CalculateAnglePotentialEnergy" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult*math.Pi*math.Pi != result {
			t.Errorf("CalculateAnglePotentialEnergy() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateProperDihedralAngleEnergy
func TestCalculateProperDihedralAngleEnergy(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateProperDihedralAngleEnergy" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateProperDihedralAngleEnergy" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateProperDihedralAngleEnergy/" + "input/" + inputFile.Name())

		kd := convertStringToFloatSlice(pair[0])[0]
		phi := convertStringToFloatSlice(pair[0])[1] * math.Pi
		pn := convertStringToFloatSlice(pair[0])[2]
		phase := convertStringToFloatSlice(pair[0])[3]

		// function
		result := CalculateProperDihedralAngleEnergy(kd, phi, pn, phase)

		// read output
		out, _ := readFileline("Tests/CalculateProperDihedralAngleEnergy" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("CalculateProperDihedralAngleEnergy() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateBondForce
func TestCalculateBondForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateBondForce" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateBondForce" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateBondForce/" + "input/" + inputFile.Name())
		k := convertStringToFloatSlice(pair[0])[0]
		r_0 := convertStringToFloatSlice(pair[0])[1]

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[1])[0]
		atom1.position.y = convertStringToFloatSlice(pair[1])[1]
		atom1.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[2])[0]
		atom2.position.y = convertStringToFloatSlice(pair[2])[1]
		atom2.position.z = convertStringToFloatSlice(pair[2])[2]

		r := Distance(atom1.position, atom2.position)
		// function
		result := CalculateBondForce(k, r, r_0, &atom1, &atom2)
		// read output
		out, _ := readFileline("Tests/CalculateBondForce" + "/output/" + outputFiles[i].Name())
		var realResult TriTuple
		realResult.x = math.Round(convertStringToFloatSlice(out[0])[0] * math.Pow10(4))
		realResult.y = math.Round(convertStringToFloatSlice(out[0])[1] * math.Pow10(4))
		realResult.z = math.Round(convertStringToFloatSlice(out[0])[2] * math.Pow10(4))

		// compare
		if realResult.x != math.Round(result.x*math.Pow10(4)) || realResult.y != math.Round(result.y*math.Pow10(4)) || realResult.z != math.Round(result.z*math.Pow10(4)) {
			t.Errorf("CalculateBondForce() = %v(%v ,%v, %v), want %v", result, math.Round(result.x*math.Pow10(4)), math.Round(result.y*math.Pow10(4)), math.Round(result.z*math.Pow10(4)), realResult)
		}

	}
}

// Test CalculateAngleForce
func TestCalculateAngleForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateAngleForce" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateAngleForce" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateAngleForce/" + "input/" + inputFile.Name())
		k := convertStringToFloatSlice(pair[0])[0]
		theta_0 := convertStringToFloatSlice(pair[0])[1]

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[1])[0]
		atom1.position.y = convertStringToFloatSlice(pair[1])[1]
		atom1.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[2])[0]
		atom2.position.y = convertStringToFloatSlice(pair[2])[1]
		atom2.position.z = convertStringToFloatSlice(pair[2])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[3])[0]
		atom3.position.y = convertStringToFloatSlice(pair[3])[1]
		atom3.position.z = convertStringToFloatSlice(pair[3])[2]

		theta := CalculateAngle(&atom1, &atom2, &atom3)
		// function
		result1, result2, result3 := CalculateAngleForce(k, theta, theta_0, &atom1, &atom2, &atom3)
		result1.x = math.Round(result1.x*math.Pow10(5)) / math.Pow10(5)
		result1.y = math.Round(result1.y*math.Pow10(5)) / math.Pow10(5)
		result1.z = math.Round(result1.z*math.Pow10(5)) / math.Pow10(5)
		result2.x = math.Round(result2.x*math.Pow10(5)) / math.Pow10(5)
		result2.y = math.Round(result2.y*math.Pow10(5)) / math.Pow10(5)
		result2.z = math.Round(result2.z*math.Pow10(5)) / math.Pow10(5)
		result3.x = math.Round(result3.x*math.Pow10(5)) / math.Pow10(5)
		result3.y = math.Round(result3.y*math.Pow10(5)) / math.Pow10(5)
		result3.z = math.Round(result3.z*math.Pow10(5)) / math.Pow10(5)

		// read output
		out, _ := readFileline("Tests/CalculateAngleForce" + "/output/" + outputFiles[i].Name())
		var realResult1, realResult2, realResult3 TriTuple
		realResult1.x = convertStringToFloatSlice(out[0])[0]
		realResult1.y = convertStringToFloatSlice(out[0])[1]
		realResult1.z = convertStringToFloatSlice(out[0])[2]

		realResult2.x = convertStringToFloatSlice(out[1])[0]
		realResult2.y = convertStringToFloatSlice(out[1])[1]
		realResult2.z = convertStringToFloatSlice(out[1])[2]

		realResult3.x = convertStringToFloatSlice(out[2])[0]
		realResult3.y = convertStringToFloatSlice(out[2])[1]
		realResult3.z = convertStringToFloatSlice(out[2])[2]

		// compare
		if result1 != realResult1 || result2 != realResult2 || result3 != realResult3 {
			t.Errorf("CalculateAngleForce() = (%v ,%v, %v), want (%v, %v, %v)", result1, result2, result3, realResult1, realResult2, realResult3)
		}

	}
}

// Test DerivateAnglePositionX
func TestDerivateAnglePositionX(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/DerivateAnglePositionX" + "/input")
	outputFiles := ReadDirectory("Tests/DerivateAnglePositionX" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/DerivateAnglePositionX/" + "input/" + inputFile.Name())

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[0])[0]
		atom1.position.y = convertStringToFloatSlice(pair[0])[1]
		atom1.position.z = convertStringToFloatSlice(pair[0])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[1])[0]
		atom2.position.y = convertStringToFloatSlice(pair[1])[1]
		atom2.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[2])[0]
		atom3.position.y = convertStringToFloatSlice(pair[2])[1]
		atom3.position.z = convertStringToFloatSlice(pair[2])[2]

		theta := CalculateAngle(&atom1, &atom2, &atom3)
		// function
		result := math.Round(DerivateAnglePositionX(&atom1, &atom2, &atom3, theta)*math.Pow10(4)) / math.Pow10(4)
		// read output
		out, _ := readFileline("Tests/DerivateAnglePositionX" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("DerivateAnglePositionX() = %v, want %v", result, realResult)
		}

	}
}

// Test DerivateAnglePositionY
func TestDerivateAnglePositionY(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/DerivateAnglePositionY" + "/input")
	outputFiles := ReadDirectory("Tests/DerivateAnglePositionY" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/DerivateAnglePositionY/" + "input/" + inputFile.Name())

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[0])[0]
		atom1.position.y = convertStringToFloatSlice(pair[0])[1]
		atom1.position.z = convertStringToFloatSlice(pair[0])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[1])[0]
		atom2.position.y = convertStringToFloatSlice(pair[1])[1]
		atom2.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[2])[0]
		atom3.position.y = convertStringToFloatSlice(pair[2])[1]
		atom3.position.z = convertStringToFloatSlice(pair[2])[2]

		theta := CalculateAngle(&atom1, &atom2, &atom3)
		// function
		result := math.Round(DerivateAnglePositionY(&atom1, &atom2, &atom3, theta)*math.Pow10(4)) / math.Pow10(4)
		// read output
		out, _ := readFileline("Tests/DerivateAnglePositionY" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("DerivateAnglePositionY() = %v, want %v", result, realResult)
		}

	}
}

// Test DerivateAnglePositionZ
func TestDerivateAnglePositionZ(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/DerivateAnglePositionZ" + "/input")
	outputFiles := ReadDirectory("Tests/DerivateAnglePositionZ" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/DerivateAnglePositionZ/" + "input/" + inputFile.Name())

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[0])[0]
		atom1.position.y = convertStringToFloatSlice(pair[0])[1]
		atom1.position.z = convertStringToFloatSlice(pair[0])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[1])[0]
		atom2.position.y = convertStringToFloatSlice(pair[1])[1]
		atom2.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[2])[0]
		atom3.position.y = convertStringToFloatSlice(pair[2])[1]
		atom3.position.z = convertStringToFloatSlice(pair[2])[2]

		theta := CalculateAngle(&atom1, &atom2, &atom3)
		// function
		result := math.Round(DerivateAnglePositionZ(&atom1, &atom2, &atom3, theta)*math.Pow10(4)) / math.Pow10(4)
		// read output
		out, _ := readFileline("Tests/DerivateAnglePositionZ" + "/output/" + outputFiles[i].Name())
		var realResult float64
		realResult = convertStringToFloatSlice(out[0])[0]

		// compare
		if realResult != result {
			t.Errorf("DerivateAnglePositionZ() = %v, want %v", result, realResult)
		}

	}
}

// Test CalculateProperDihedralsForce
func TestCalculateProperDihedralsForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateProperDihedralsForce" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateProperDihedralsForce" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		pair, _ := readFileline("Tests/CalculateProperDihedralsForce/" + "input/" + inputFile.Name())
		//kd, phi, pn, phase
		kd := convertStringToFloatSlice(pair[0])[0]
		pn := convertStringToFloatSlice(pair[0])[1]
		phase := convertStringToFloatSlice(pair[0])[2]

		var atom1 Atom
		atom1.position.x = convertStringToFloatSlice(pair[1])[0]
		atom1.position.y = convertStringToFloatSlice(pair[1])[1]
		atom1.position.z = convertStringToFloatSlice(pair[1])[2]

		var atom2 Atom
		atom2.position.x = convertStringToFloatSlice(pair[2])[0]
		atom2.position.y = convertStringToFloatSlice(pair[2])[1]
		atom2.position.z = convertStringToFloatSlice(pair[2])[2]

		var atom3 Atom
		atom3.position.x = convertStringToFloatSlice(pair[3])[0]
		atom3.position.y = convertStringToFloatSlice(pair[3])[1]
		atom3.position.z = convertStringToFloatSlice(pair[3])[2]

		var atom4 Atom
		atom4.position.x = convertStringToFloatSlice(pair[4])[0]
		atom4.position.y = convertStringToFloatSlice(pair[4])[1]
		atom4.position.z = convertStringToFloatSlice(pair[4])[2]

		phi := CalculateDihedralAngle(&atom1, &atom2, &atom3, &atom4)
		// function
		result1, result2, result3, result4 := CalculateProperDihedralsForce(kd, phi, pn, phase, &atom1, &atom2, &atom3, &atom4)

		result1.x = math.Round(result1.x*math.Pow10(5)) / math.Pow10(5)
		result1.y = math.Round(result1.y*math.Pow10(5)) / math.Pow10(5)
		result1.z = math.Round(result1.z*math.Pow10(5)) / math.Pow10(5)
		result2.x = math.Round(result2.x*math.Pow10(5)) / math.Pow10(5)
		result2.y = math.Round(result2.y*math.Pow10(5)) / math.Pow10(5)
		result2.z = math.Round(result2.z*math.Pow10(5)) / math.Pow10(5)
		result3.x = math.Round(result3.x*math.Pow10(5)) / math.Pow10(5)
		result3.y = math.Round(result3.y*math.Pow10(5)) / math.Pow10(5)
		result3.z = math.Round(result3.z*math.Pow10(5)) / math.Pow10(5)
		result4.x = math.Round(result4.x*math.Pow10(5)) / math.Pow10(5)
		result4.y = math.Round(result4.y*math.Pow10(5)) / math.Pow10(5)
		result4.z = math.Round(result4.z*math.Pow10(5)) / math.Pow10(5)

		// read output
		out, _ := readFileline("Tests/CalculateProperDihedralsForce" + "/output/" + outputFiles[i].Name())
		var realResult1, realResult2, realResult3, realResult4 TriTuple
		realResult1.x = convertStringToFloatSlice(out[0])[0]
		realResult1.y = convertStringToFloatSlice(out[0])[1]
		realResult1.z = convertStringToFloatSlice(out[0])[2]

		realResult2.x = convertStringToFloatSlice(out[1])[0]
		realResult2.y = convertStringToFloatSlice(out[1])[1]
		realResult2.z = convertStringToFloatSlice(out[1])[2]

		realResult3.x = convertStringToFloatSlice(out[2])[0]
		realResult3.y = convertStringToFloatSlice(out[2])[1]
		realResult3.z = convertStringToFloatSlice(out[2])[2]

		realResult4.x = convertStringToFloatSlice(out[3])[0]
		realResult4.y = convertStringToFloatSlice(out[3])[1]
		realResult4.z = convertStringToFloatSlice(out[3])[2]

		// compare
		if result1 != realResult1 || result2 != realResult2 || result3 != realResult3 || result4 != realResult4 {
			t.Errorf("CalculateProperDihedralsForce() = (%v ,%v, %v, %v), want (%v, %v, %v, %v)", result1, result2, result3, result4, realResult1, realResult2, realResult3, realResult4)
		}

	}
}

// Test SearchParameter
func TestSearchParameter(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/SearchParameter" + "/input")
	outputFiles := ReadDirectory("Tests/SearchParameter" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		var parameterset parameterDatabase
		atomlist, _ := ReadAtoms("Tests/SearchParameter/" + "input/" + inputFile.Name())
		var result []float64
		if len(atomlist) == 2 {
			parameterset, _ = ReadParameterFile("../data/ffbonded_bondtypes.itp")
			result = SearchParameter(len(atomlist), parameterset, &atomlist[0], &atomlist[1])
		} else if len(atomlist) == 3 {
			parameterset, _ = ReadParameterFile("../data/ffbonded_angletypes.itp")
			result = SearchParameter(len(atomlist), parameterset, &atomlist[0], &atomlist[1], &atomlist[2])
		} else if len(atomlist) == 4 {
			parameterset, _ = ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
			result = SearchParameter(len(atomlist), parameterset, &atomlist[0], &atomlist[1], &atomlist[2], &atomlist[3])
		}

		out, _ := readFileline("Tests/SearchParameter" + "/output/" + outputFiles[i].Name())
		realResult := convertStringToFloatSlice(out[0])

		// compare
		if len(result) == len(realResult) {
			for j := range realResult {
				if realResult[j] != result[j] {
					t.Errorf("SearchParameter() = %v, want %v", result, realResult)
				}
			}
		} else {
			t.Errorf("Mismatch length, (%v, %v) in %v with %v and %v", len(realResult), len(result), outputFiles[i].Name(), realResult, result)
		}
	}

}

// Test CalculateTotalEnergyForce
func TestCalculateTotalEnergyForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CalculateTotalEnergyForce" + "/input")
	outputFiles := ReadDirectory("Tests/CalculateTotalEnergyForce" + "/output")

	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/CalculateTotalEnergyForce/" + "input/" + inputFile.Name())

		// function
		result1, result2 := CalculateTotalEnergyForce(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)
		result1 = ConvertFloat(result1)

		for t := range result2 {
			res := ConvertTriTuple(*result2[t])
			result2[t] = &res
		}

		// real output
		realResult1, realResult2, _ := ReadFloatMapIntTriTuple("Tests/CalculateTotalEnergyForce" + "/output/" + outputFiles[i].Name())
		realResult1 = ConvertFloat(realResult1)
		for t := range realResult2 {
			res := ConvertTriTuple(*realResult2[t])
			realResult2[t] = &res
		}

		// compare
		if realResult1 != result1 {
			t.Errorf("Energy in CalculateTotalEnergyForce() = %v, want %v in %v", result1, realResult1, outputFiles[i])

		}

		if len(result2) == len(realResult2) {
			for key := range result2 {
				if (*realResult2[key]) != (*result2[key]) {
					t.Errorf("Energy in CalculateTotalEnergyForce() = %v, want %v in %v", result2, realResult2, outputFiles[i])
				}
			}
		} else {
			t.Errorf("Mismatch")
		}

	}
}

// Test CombineEnergyAndForce
func TestCombineEnergyAndForce(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/CombineEnergyAndForce" + "/input")
	outputFiles := ReadDirectory("Tests/CombineEnergyAndForce" + "/output")

	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/CombineEnergyAndForce/" + "input/" + inputFile.Name())

		// function
		result1, result2 := CombineEnergyAndForce(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)

		result1 = ConvertFloat(result1)

		for t := range result2 {
			res := ConvertTriTuple(*result2[t])
			result2[t] = &res
		}

		// real output
		realResult1, realResult2, _ := ReadFloatMapIntTriTuple("Tests/CombineEnergyAndForce" + "/output/" + outputFiles[i].Name())

		realResult1 = ConvertFloat(realResult1)
		for t := range realResult2 {
			res := ConvertTriTuple(*realResult2[t])
			realResult2[t] = &res
		}

		// compare
		if realResult1 != result1 {
			t.Errorf("Energy in CombineEnergyAndForce() = %v, want %v in %v", result1, realResult1, outputFiles[i])

		}

		if len(result2) == len(realResult2) {
			for key := range result2 {
				if (*realResult2[key]) != (*result2[key]) {
					t.Errorf("Energy in CombineEnergyAndForce() = %v, want %v in %v", result2, realResult2, outputFiles[i])
				}
			}
		} else {
			t.Errorf("Mismatch")
		}

	}
}

// Test SteepestDescent
func TestSteepestDescent(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/SteepestDescent" + "/input")
	outputFiles := ReadDirectory("Tests/SteepestDescent" + "/output")
	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/SteepestDescent/" + "input/" + inputFile.Name())
		h := 0.01
		_, forceMap := CombineEnergyAndForce(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter)

		result := CopyProtein(&protein)
		// function
		SteepestDescent(result, h, forceMap)
		temp := ConvertProtein(*result)
		result = &temp

		//real output
		realResult, _ := ReadOneProtein("Tests/SteepestDescent/" + "output/" + outputFiles[i].Name())
		realResult = ConvertProtein(realResult)

		// compare

		if realResult.Name != result.Name {
			t.Errorf("SteepestDescent() = %v, want %v", result, realResult)
		}

		// compare residues
		for i := range realResult.Residue {
			if len(realResult.Residue[i].Atoms) == len(result.Residue[i].Atoms) {
				if realResult.Residue[i].ChainID != result.Residue[i].ChainID || realResult.Residue[i].Name != result.Residue[i].Name || realResult.Residue[i].ID != result.Residue[i].ID {
					t.Errorf("SteepestDescent() = %v, want %v", result, realResult)
				}
				for j := range realResult.Residue[i].Atoms {
					if (*realResult.Residue[i].Atoms[j]) != (*result.Residue[i].Atoms[j]) {
						t.Errorf("SteepestDescent() = %v, want %v, in %v .", result, realResult, outputFiles[i].Name())
					}
				}
			} else {
				t.Errorf("Number mismatch")
			}
		}

	}

}

// Test PerformEnergyMinimization
func TestPerformEnergyMinimization(t *testing.T) {
	// input and outputfiles
	inputFiles := ReadDirectory("Tests/PerformEnergyMinimization" + "/input")
	outputFiles := ReadDirectory("Tests/PerformEnergyMinimization" + "/output")

	residueParameterBondValue, error := ReadAminoAcidsPara("../data/aminoacids_revised.rtp")
	Check(error)
	residueParameterOtherValue, error := ReadAminoAcidsPara("../data/aminoacids.rtp")
	Check(error)
	bondParameter, error := ReadParameterFile("../data/ffbonded_bondtypes.itp")
	Check(error)
	angleParameter, error := ReadParameterFile("../data/ffbonded_angletypes.itp")
	Check(error)
	dihedralParameter, error := ReadParameterFile("../data/ffbonded_dihedraltypes.itp")
	Check(error)
	nonbondedParameter, error := ReadParameterFile("../data/ffnonbonded_nonbond_params.itp")
	Check(error)
	pairtypesParameter, error := ReadParameterFile("../data/ffnonbonded_pairtypes.itp")
	Check(error)

	for i, inputFile := range inputFiles {
		// read input
		protein, _ := ReadOneProtein("Tests/PerformEnergyMinimization/" + "input/" + inputFile.Name())
		iteration := 1
		// function
		result := PerformEnergyMinimization(&protein, residueParameterBondValue, residueParameterOtherValue, bondParameter, angleParameter, dihedralParameter, nonbondedParameter, pairtypesParameter, iteration)
		temp := ConvertProtein(*result)
		result = &temp

		// real output
		realResult, _ := ReadOneProtein("Tests/PerformEnergyMinimization" + "/output/" + outputFiles[i].Name())
		realResult = ConvertProtein(realResult)

		// compare
		if realResult.Name != result.Name {
			t.Errorf("PerformEnergyMinimization() = %v, want %v", result, realResult)
		}

		// compare residues
		for i := range realResult.Residue {
			if len(realResult.Residue[i].Atoms) == len(result.Residue[i].Atoms) {
				if realResult.Residue[i].ChainID != result.Residue[i].ChainID || realResult.Residue[i].Name != result.Residue[i].Name || realResult.Residue[i].ID != result.Residue[i].ID {
					t.Errorf("PerformEnergyMinimization() = %v, want %v", result, realResult)
				}
				for j := range realResult.Residue[i].Atoms {
					if (*realResult.Residue[i].Atoms[j]) != (*result.Residue[i].Atoms[j]) {
						t.Errorf("PerformEnergyMinimization() = %v, want %v, in %v with atom %v in %v.", result, realResult, outputFiles[i].Name(), *realResult.Residue[i].Atoms[j], *result.Residue[i].Atoms[j])
					}
				}
			} else {
				t.Errorf("Number mismatch")
			}
		}

	}
}

// //////////
// Readtest area
// //////////

// ReadCalculateElectricPotentialEnergyTests reads the input and output file for CalculateElectricPotentialEnergyTest
func ReadCalculateElectricPotentialEnergyTests(directory string) []CalculateElectricPotentialEnergyTest {

	// Read input file
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CalculateElectricPotentialEnergyTest, numFiles)
	for i, inputFile := range inputFiles {
		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].a1.charge = convertStringToFloatSlice(pair[0])[0]
		tests[i].a2.charge = convertStringToFloatSlice(pair[0])[1]
		tests[i].r = convertStringToFloatSlice(pair[1])[0]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFiles := range outputFiles {
		out, _ := readFileline(directory + "output/" + outputFiles.Name())
		tests[i].result = convertStringToFloatSlice(out[0])[0]
	}

	return tests
}

// ReadCalculateLJPotentialEnergyTest reads the input and output file for CalculateLJPotentialEnergyTest
func ReadCalculateLJPotentialEnergyTest(directory string) []CalculateLJPotentialEnergyTest {

	// Read input file
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CalculateLJPotentialEnergyTest, numFiles)
	for i, inputFile := range inputFiles {
		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].B = convertStringToFloatSlice(pair[0])[0]
		tests[i].A = convertStringToFloatSlice(pair[0])[1]
		tests[i].r = convertStringToFloatSlice(pair[1])[0]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFiles := range outputFiles {
		out, _ := readFileline(directory + "output/" + outputFiles.Name())
		tests[i].result = convertStringToFloatSlice(out[0])[0]
	}

	return tests
}

// ReadCalculateElectricForceTests reads the input and output file for CalculateElectricForceTest
func ReadCalculateElectricForceTests(directory string) []CalculateElectricForceTest {

	// Read input file
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CalculateElectricForceTest, numFiles)
	for i, inputFile := range inputFiles {
		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].a1.charge = convertStringToFloatSlice(pair[0])[0]
		tests[i].a2.charge = convertStringToFloatSlice(pair[0])[1]
		tests[i].r = convertStringToFloatSlice(pair[1])[0]
		tests[i].a1.position.x = convertStringToFloatSlice(pair[2])[0]
		tests[i].a1.position.y = convertStringToFloatSlice(pair[2])[1]
		tests[i].a1.position.z = convertStringToFloatSlice(pair[2])[2]
		tests[i].a2.position.x = convertStringToFloatSlice(pair[3])[0]
		tests[i].a2.position.y = convertStringToFloatSlice(pair[3])[1]
		tests[i].a2.position.z = convertStringToFloatSlice(pair[3])[2]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFiles := range outputFiles {
		out, _ := readFileline(directory + "output/" + outputFiles.Name())
		tests[i].result.x = convertStringToFloatSlice(out[0])[0]
		tests[i].result.y = convertStringToFloatSlice(out[1])[0]
		tests[i].result.z = convertStringToFloatSlice(out[2])[0]
	}

	return tests
}

// ReadCalculateLJEnergyTest reads the input and output file for CalculateLJEnergyTest
func ReadCalculateLJForceTest(directory string) []CalculateLJForceTest {

	// Read input file
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CalculateLJForceTest, numFiles)
	for i, inputFile := range inputFiles {
		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].B = convertStringToFloatSlice(pair[0])[0]
		tests[i].A = convertStringToFloatSlice(pair[0])[1]
		tests[i].r = convertStringToFloatSlice(pair[1])[0]
		tests[i].a1.position.x = convertStringToFloatSlice(pair[2])[0]
		tests[i].a1.position.y = convertStringToFloatSlice(pair[2])[1]
		tests[i].a1.position.z = convertStringToFloatSlice(pair[2])[2]
		tests[i].a2.position.x = convertStringToFloatSlice(pair[3])[0]
		tests[i].a2.position.y = convertStringToFloatSlice(pair[3])[1]
		tests[i].a2.position.z = convertStringToFloatSlice(pair[3])[2]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFiles := range outputFiles {
		out, _ := readFileline(directory + "output/" + outputFiles.Name())
		tests[i].result.x = convertStringToFloatSlice(out[0])[0]
		tests[i].result.y = convertStringToFloatSlice(out[1])[0]
		tests[i].result.z = convertStringToFloatSlice(out[2])[0]
	}

	return tests
}

// func ReadUpdateAcceleration read the input and output file for UpdateAccelerationTest
func ReadUpdateAccelerationTests(directory string) []UpdateAccelerationTest {

	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]UpdateAccelerationTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's map

		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].force.x = convertStringToFloatSlice(pair[0])[0]
		tests[i].force.y = convertStringToFloatSlice(pair[0])[1]
		tests[i].force.z = convertStringToFloatSlice(pair[0])[2]
		tests[i].a.mass = convertStringToFloatSlice(pair[1])[0]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		out, _ := readFileline(directory + "output/" + outputFile.Name())
		tests[i].resultaccel.x = convertStringToFloatSlice(out[0])[0]
		tests[i].resultaccel.y = convertStringToFloatSlice(out[0])[1]
		tests[i].resultaccel.z = convertStringToFloatSlice(out[0])[2]

	}

	return tests
}

// func ReadUpdatePositionTests read the input and output file for UpdatePositionTest
func ReadUpdatePositionTests(directory string) []UpdatePositionTest {

	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]UpdatePositionTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's map

		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].a.position.x = convertStringToFloatSlice(pair[0])[0]
		tests[i].a.position.y = convertStringToFloatSlice(pair[0])[1]
		tests[i].a.position.z = convertStringToFloatSlice(pair[0])[2]
		tests[i].oldAcceleration.x = convertStringToFloatSlice(pair[1])[0]
		tests[i].oldAcceleration.y = convertStringToFloatSlice(pair[1])[1]
		tests[i].oldAcceleration.z = convertStringToFloatSlice(pair[1])[2]
		tests[i].oldVelocity.x = convertStringToFloatSlice(pair[2])[0]
		tests[i].oldVelocity.y = convertStringToFloatSlice(pair[2])[1]
		tests[i].oldVelocity.z = convertStringToFloatSlice(pair[2])[2]
		tests[i].time = convertStringToFloatSlice(pair[3])[0]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		out, _ := readFileline(directory + "output/" + outputFile.Name())
		tests[i].result.x = convertStringToFloatSlice(out[0])[0]
		tests[i].result.y = convertStringToFloatSlice(out[0])[1]
		tests[i].result.z = convertStringToFloatSlice(out[0])[2]

	}

	return tests
}

// func ReadUpdateVelocityTests read the input and output file for UpdateVelocityTest
func ReadUpdateVelocityTests(directory string) []UpdateVelocityTest {

	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]UpdateVelocityTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's map

		pair, _ := readFileline(directory + "input/" + inputFile.Name())
		tests[i].a.velocity.x = convertStringToFloatSlice(pair[0])[0]
		tests[i].a.velocity.y = convertStringToFloatSlice(pair[0])[1]
		tests[i].a.velocity.z = convertStringToFloatSlice(pair[0])[2]
		tests[i].a.accelerated.x = convertStringToFloatSlice(pair[1])[0]
		tests[i].a.accelerated.y = convertStringToFloatSlice(pair[1])[1]
		tests[i].a.accelerated.z = convertStringToFloatSlice(pair[1])[2]
		tests[i].oldAcceleration.x = convertStringToFloatSlice(pair[2])[0]
		tests[i].oldAcceleration.y = convertStringToFloatSlice(pair[2])[1]
		tests[i].oldAcceleration.z = convertStringToFloatSlice(pair[2])[2]
		tests[i].time = convertStringToFloatSlice(pair[3])[0]
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")

	//ensure same number of input and output files
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		out, _ := readFileline(directory + "output/" + outputFile.Name())
		tests[i].result.x = convertStringToFloatSlice(out[0])[0]
		tests[i].result.y = convertStringToFloatSlice(out[0])[1]
		tests[i].result.z = convertStringToFloatSlice(out[0])[2]

	}

	return tests
}

// read atoms in the file
func ReadAtoms(filename string) ([]Atom, error) {
	// read lines
	lines, err := readFileline(filename)
	if err != nil {
		return nil, err
	}

	// read atom in each line
	var atoms []Atom
	for _, line := range lines {
		// Convert the line to a slice of float64
		line := strings.Split(line, " ")
		if len(line) != 16 {
			panic("Error: number of vaiables in atom input do not match!")
		}
		atom := ReadOneAtom(line)

		atoms = append(atoms, atom)
	}
	return atoms, nil
}

// read one atom in a line
func ReadOneAtom(line []string) Atom {
	var atom Atom
	atom.index, _ = strconv.Atoi(line[0])
	atom.position.x, _ = strconv.ParseFloat(line[1], 64)
	atom.position.y, _ = strconv.ParseFloat(line[2], 64)
	atom.position.z, _ = strconv.ParseFloat(line[3], 64)
	atom.velocity.x, _ = strconv.ParseFloat(line[4], 64)
	atom.velocity.y, _ = strconv.ParseFloat(line[5], 64)
	atom.velocity.z, _ = strconv.ParseFloat(line[6], 64)
	atom.force.x, _ = strconv.ParseFloat(line[7], 64)
	atom.force.y, _ = strconv.ParseFloat(line[8], 64)
	atom.force.z, _ = strconv.ParseFloat(line[9], 64)
	atom.accelerated.x, _ = strconv.ParseFloat(line[10], 64)
	atom.accelerated.y, _ = strconv.ParseFloat(line[11], 64)
	atom.accelerated.z, _ = strconv.ParseFloat(line[12], 64)
	atom.mass, _ = strconv.ParseFloat(line[13], 64)
	atom.element = line[14]
	atom.charge, _ = strconv.ParseFloat(line[15], 64)

	return atom
}

// read a residue in some lines
func ReadOneResidue(lines []string) Residue {
	var residue Residue
	var atoms []*Atom

	// read fields except atoms
	for i, _ := range lines {
		line := strings.Split(lines[i], " ")
		if i == 0 {
			residue.Name = line[0]
			residue.ChainID = line[1]
			residue.ID, _ = strconv.Atoi(line[2])

			continue
		}

		// read atoms
		atom := ReadOneAtom(line)
		atoms = append(atoms, &atom)

	}
	residue.Atoms = atoms
	return residue
}

// read residues in one file
func ReadResidues(filename string) ([]Residue, error) {
	// read lines
	lines, err := readFileline(filename)
	if err != nil {
		return nil, err
	}

	var residues []Residue

	for i := 0; i < len(lines); i++ {
		line1 := strings.Split(lines[i], " ")
		AtomLen, _ := strconv.Atoi(line1[3])
		residue := ReadOneResidue(lines[i : i+AtomLen+1])
		residues = append(residues, residue)

	}

	return residues, err
}

// read one protein
func ReadOneProtein(filename string) (Protein, error) {
	var protein Protein
	// read lines in the file
	lines, err := readFileline(filename)
	if err != nil {
		return protein, err
	}

	// read protein name
	firstLine := strings.Split(lines[0], " ")
	protein.Name = firstLine[0]
	residueLen, _ := strconv.Atoi(firstLine[1])
	startIndex := 1
	// read residues
	for j := 0; j < residueLen; j++ {
		line := strings.Split(lines[startIndex], " ")
		AtomLen, _ := strconv.Atoi(line[3])
		// read one residue in a line
		residue := ReadOneResidue(lines[startIndex : startIndex+AtomLen+1])
		protein.Residue = append(protein.Residue, &residue)
		startIndex += AtomLen + 1
	}

	return protein, err
}

// read proteins
func ReadProteins(filename string) ([]Protein, error) {
	var proteins []Protein
	// read lines
	lines, err := readFileline(filename)
	if err != nil {
		return proteins, err
	}

	i := 0
	for i < len(lines) {
		// Check if the current line starts a new protein
		if strings.HasPrefix(lines[i], "Protein") {
			firstLine := strings.Split(lines[i], " ")
			protein := Protein{Name: firstLine[0]}
			residueLen, _ := strconv.Atoi(firstLine[1])
			i++ // Move to the first residue line
			for j := 0; j < residueLen; j++ {
				line := strings.Split(lines[i], " ")
				AtomLen, _ := strconv.Atoi(line[3])
				residue := ReadOneResidue(lines[i : i+AtomLen+1])
				protein.Residue = append(protein.Residue, &residue)
				i += AtomLen + 1
			}
			proteins = append(proteins, protein)
		} else {
			i++ // Skip any unrelated lines
		}
	}
	return proteins, nil
}

// read directory
func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// func readFileline take file as input
// and return  each line as []string
func readFileline(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())

	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return lines, nil
}

// func convertStringToFloatSlice takes string as input
// return a slice of float64
func convertStringToFloatSlice(input string) []float64 {
	// Split the input string by spaces
	fields := strings.Fields(input)
	var floatValues []float64

	// Convert each field to a float64
	for _, field := range fields {
		floatValue, _ := strconv.ParseFloat(field, 64)
		floatValues = append(floatValues, floatValue)
	}

	return floatValues
}

// ReadTriTuples reads a file and parses its contents into a slice of TriTuple.
// Each line in the file should contain three float values, separated by spaces.
func ReadTriTuples(filename string) ([]TriTuple, error) {
	lines, err := readFileline(filename)
	if err != nil {
		return nil, err
	}

	var triTuples []TriTuple
	for _, line := range lines {
		// Convert the line to a slice of float64
		floatValues := convertStringToFloatSlice(line)
		if len(floatValues) != 3 {
			return nil, fmt.Errorf("invalid line format: %s", line)
		}

		// Create a TriTuple and add it to the slice
		triTuple := TriTuple{
			x: floatValues[0],
			y: floatValues[1],
			z: floatValues[2],
		}
		triTuples = append(triTuples, triTuple)
	}
	return triTuples, nil
}

// read a float and a map containing pointers of TriTuple with int as key
func ReadFloatMapIntTriTuple(filename string) (float64, map[int]*TriTuple, error) {
	// read lines
	lines, err := readFileline(filename)
	if err != nil {
		return math.NaN(), nil, err
	}
	var energy float64
	Map := make(map[int]*TriTuple)
	for i := range lines {
		line := strings.Split(lines[i], " ")
		if len(line) == 1 { // read float64
			energy, _ = strconv.ParseFloat(line[0], 64)
		}
		if len(line) == 4 { // read TriTuple
			key, _ := strconv.Atoi(line[0])
			var value TriTuple
			value.x, _ = strconv.ParseFloat(line[1], 64)
			value.y, _ = strconv.ParseFloat(line[2], 64)
			value.z, _ = strconv.ParseFloat(line[3], 64)
			Map[key] = &TriTuple{x: value.x, y: value.y, z: value.z}
		}
	}
	return energy, Map, nil
}

// convert

func ConvertFloat(n float64) float64 {
	return math.Round(n*math.Pow10(5)) / math.Pow10(5)
}

func ConvertTriTuple(tri TriTuple) TriTuple {
	return TriTuple{x: ConvertFloat(tri.x), y: ConvertFloat(tri.y), z: ConvertFloat(tri.z)}
}

func ConvertProtein(p Protein) Protein {
	for i := range p.Residue {
		temp := ConvertResidue(*p.Residue[i])
		p.Residue[i] = &temp
	}

	return p
}

func ConvertResidue(r Residue) Residue {
	for i := range r.Atoms {
		temp := ConvertAtom(*r.Atoms[i])
		r.Atoms[i] = &temp
	}

	return r
}

func ConvertAtom(a Atom) Atom {
	a.force = ConvertTriTuple(a.force)
	a.position = ConvertTriTuple(a.position)
	a.velocity = ConvertTriTuple(a.velocity)
	a.accelerated = ConvertTriTuple(a.accelerated)
	return a
}
