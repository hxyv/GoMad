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
