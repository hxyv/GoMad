package main

import (
	"bufio"
	"fmt"
	"io/fs"
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

// //////////
// Test area
// //////////
func TestUpdateAcceleration(t *testing.T) {
	// Read in all tests from the Tests/Distance directory and run them
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

func TestUpdatePosition(t *testing.T) {
	// Read in all tests from the Tests/Distance directory and run them
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

func TestUpdateVelocity(t *testing.T) {
	// Read in all tests from the Tests/Distance directory and run them
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

func TestCalculateVector(t *testing.T) {
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

		if realResult != result {
			t.Errorf("CalculateVector() = %v, want %v", result, realResult)
		}

	}
}

func TestCopyTriTuple(t *testing.T) {
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

		if realResult != result {
			t.Errorf("TriTuple() = %v, want %v", result, realResult)
		}

	}
}

func TestCopyAtom(t *testing.T) {
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

		if realResult != *result {
			t.Errorf("CopyAtom() = %v, want %v", result, realResult)
		}

	}
}

/*
func TestCopyResidue(t *testing.T) {
	inputFiles := ReadDirectory("Tests/CopyResidue" + "/input")
	outputFiles := ReadDirectory("Tests/CopyResidue" + "/output")

	for i, inputFile := range inputFiles {
		// read input
		atomlist, _ := ReadAtoms("Tests/CopyResidue/" + "input/" + inputFile.Name())
		atom := atomlist[0]

		// function
		result := CopyAtom(&atom)

		// read output
		atomlist1, _ := ReadAtoms("Tests/CopyResidue/" + "output/" + outputFiles[i].Name())
		realResult := atomlist1[0]

		if realResult != *result {
			t.Errorf("CopyAtom() = %v, want %v", result, realResult)
		}

	}
}
*/

func TestCross(t *testing.T) {
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

		if realResult != result {
			t.Errorf("Cross() = %v, want %v", result, realResult)
		}

	}
}

func TestCalculateDerivate(t *testing.T) {
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

		phi := convertStringToFloatSlice(pair[2])[0]
		// function
		der1, der2, der3 := CalculateDerivate(v1, v2, phi)

		// read output
		out, _ := readFileline("Tests/CalculateDerivate" + "/output/" + outputFiles[i].Name())
		result1 := convertStringToFloatSlice(out[0])[0]
		result2 := convertStringToFloatSlice(out[0])[1]
		result3 := convertStringToFloatSlice(out[0])[2]

		if der1 != result1 || der2 != result2 || der3 != result3 {
			t.Errorf("CalculateDerivate() = %v. %v, %v, want (%v, %v, %v)", der1, der2, der3, result1, result2, result3)
		}

	}
}

func TestMagnitude(t *testing.T) {
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

		if realResult != result {
			t.Errorf("magnitude() = %v, want %v", result, realResult)
		}

	}
}

func TestDot(t *testing.T) {
	inputFiles := ReadDirectory("Tests/Dot" + "/input")
	outputFiles := ReadDirectory("Tests/Dot" + "/output")

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

		if realResult != result {
			t.Errorf("magnitude() = %v, want %v", result, realResult)
		}

	}
}

// //////////
// Readtest area
// //////////
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

func ReadAtoms(filename string) ([]Atom, error) {
	lines, err := readFileline(filename)
	if err != nil {
		return nil, err
	}

	var atoms []Atom
	for _, line := range lines {
		// Convert the line to a slice of float64
		line := strings.Split(line, " ")
		if len(line) != 16 {
			panic("Error: number of vaiables in atom input do not match!")
		}
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

		atoms = append(atoms, atom)
	}
	return atoms, nil
}

/*
func ReadResidues(filename string) ([]Residue, error) {
	lines, err := readFileline(filename)
	if err != nil {
		return nil, err
	}

	for i := 0; i < len(lines)/2; i++{
		var residue Residue
		line1 := strings.Split(lines[2*i], " ")
		residue.ChainID = line1[0]
        residue.ID, _ = strconv.Atoi(line1[1])

	}
}
*/

// /////
// /////
// Read func
// ////
// ////
// ReadDirectory reads in a directory and returns a slice of fs.DirEntry objects containing file info for the directory
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
