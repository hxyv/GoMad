package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// ///////////////
// ////These function are used for read protein from PDB
// ///////////////

// readProteinFromFile take a fileName as example
// return the Protein structure using the informtion of file
func readProteinFromFile(filepath string) (Protein, error) {
	file, err := os.Open(filepath)
	if err != nil {
		return Protein{}, err
	}
	defer file.Close()

	var protein Protein
	var currentResidue *Residue

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "ATOM") {
			parts := strings.Fields(line)
			if len(parts) < 11 {
				continue
			}

			atomIndex, _ := strconv.Atoi(parts[1])
			element := parts[2]
			residueName := parts[3]
			chainID := parts[4]
			residueID, _ := strconv.Atoi(parts[5])
			x, _ := strconv.ParseFloat(parts[6], 64)
			y, _ := strconv.ParseFloat(parts[7], 64)
			z, _ := strconv.ParseFloat(parts[8], 64)

			if currentResidue == nil || currentResidue.ID != residueID {
				currentResidue = &Residue{
					Name:    residueName,
					ID:      residueID,
					ChainID: chainID,
					Atoms:   []*Atom{},
				}
				protein.Residue = append(protein.Residue, currentResidue)
			}

			atom := &Atom{
				index:    atomIndex,
				position: TriTuple{x: x, y: y, z: z},
				element:  element,
			}
			currentResidue.Atoms = append(currentResidue.Atoms, atom)
		}
	}

	if err := scanner.Err(); err != nil {
		return Protein{}, err
	}

	// upload weight of each atoms
	protein.UpdateMasses(massTable)

	return protein, nil
}

func (p *Protein) UpdateMasses(massTable map[string]float64) {
	for _, residue := range p.Residue {
		for _, atom := range residue.Atoms {
			// Extract the first character of the element to match in the mass table
			baseElement := string(atom.element[0])

			if mass, found := massTable[baseElement]; found {
				atom.mass = mass
			} else {
				fmt.Printf("Warning: Mass not found for element %s (using base element %s)\n", atom.element, baseElement)
				atom.mass = 0.0 //
			}
		}
	}
}

// the mass table for the common atoms in protein
var massTable = map[string]float64{
	"H": 1.0079,
	"C": 12.0107,
	"N": 14.0067,
	"O": 15.9994,
	"S": 32.065,
	// Add more elements as needed
}

// ///////////////
// ////These function are used for read parameter for MDsimulation
// ///////////////

// **** highest level function ****
// Function ReadParameterFile take a filePath as example
// and return parameterDatabase that contains the information between each atom pairs
func ReadParameterFile(filePath string) (parameterDatabase, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return parameterDatabase{}, err
	}
	defer file.Close()

	var pairs parameterDatabase
	Firstline, _ := GetFirstLine(filePath)
	funcPosition, len, _ := FindPosition(Firstline)
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		pair, err := ParseParameterPairLine(line, funcPosition, len)
		if err != nil {
			continue
		}
		pairPointer := &pair
		pairs.atomPair = append(pairs.atomPair, pairPointer)
	}

	if err := scanner.Err(); err != nil {
		return parameterDatabase{}, err
	}

	return pairs, nil
}

// Function to parse a single line and return a parameterPair struct
func ParseParameterPairLine(line string, funcPosition, length int) (parameterPair, error) {
	// Remove any leading/trailing whitespace
	line = strings.TrimSpace(line)

	// Skip empty lines or lines that start with a comment
	if line == "" || strings.HasPrefix(line, ";") {
		return parameterPair{}, fmt.Errorf("empty or comment line")
	}

	// check the number of line
	fields := strings.Fields(line)
	// Initialize parameterPair struct and parse atom names
	atomName := fields[:funcPosition-1]
	parameter := make([]float64, length-funcPosition-1)

	var pair parameterPair
	pair.atomName = atomName
	pair.parameter = parameter

	// Parse function
	function, err := strconv.Atoi(fields[funcPosition-1])
	if err != nil {
		return parameterPair{}, err
	}
	pair.Function = function

	// Parse parameters
	for i := 0; i < len(pair.parameter); i++ {
		param, err := strconv.ParseFloat(fields[funcPosition+i], 64)
		if err != nil {
			return parameterPair{}, err
		}
		pair.parameter[i] = param
	}

	return pair, nil
}

// Function FindPosition to judge the length of the parameter and the position of func
func FindPosition(line string) (int, int, error) {
	// Check if the line starts with ";"
	if !strings.HasPrefix(line, ";") {
		return -1, 0, fmt.Errorf("line does not start with a comment: %s", line)
	}

	// Split the line into fields
	fields := strings.Fields(line)
	len := len(fields)

	// Look for "func" in the fields
	for i, field := range fields {
		if field == "func" {
			// Return the index of "func"
			return i, len, nil
		}
	}

	return -1, 0, fmt.Errorf("'func' not found in the line")
}

// function GetFirstLine used to retrive the first line of a file
func GetFirstLine(filePath string) (string, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return "", err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read and return the first line
	if scanner.Scan() {
		return scanner.Text(), nil
	}

	return "", fmt.Errorf("file does not have any lines")
}

// ///////////////
// ////These function are used for read parameter for aminoacids.rtp
// ///////////////

// ****highest level function****
func ReadAminoAcidsPara(fileName string) (map[string]residueParameter, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	//creata a map for residueParameter
	residues := make(map[string]residueParameter)
	var currentResidue *residueParameter
	section := ""

	// scan the line of the file
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, ";") {
			continue
		}

		if strings.HasPrefix(line, "[") && strings.HasSuffix(line, "]") {
			section = strings.Trim(line, "[] ")
			if section != "atoms" && section != "bonds" && section != "angles" && section != "dihedrals" && section != "impropers" &&
				section != "all_dihedrals" && section != "HH14" && section != "RemoveDih" && section != "bondedtypes" {
				currentResidue = &residueParameter{name: section}
				residues[section] = *currentResidue
			}
			continue
		}

		if currentResidue == nil {
			continue
		}
		// append the information we need for the residueParameter
		//atoms,bonds,angles,dihedrals
		switch section {
		case "atoms":
			parts := strings.Fields(line)
			if len(parts) >= 4 {
				x, _ := strconv.ParseFloat(parts[2], 64)
				y, _ := strconv.ParseFloat(parts[3], 64)
				currentResidue.atoms = append(currentResidue.atoms, &atoms{
					atoms: []string{parts[0], parts[1]},
					x:     x,
					y:     y,
				})
			}
		case "bonds":
			parts := strings.Fields(line)
			if len(parts) >= 3 {
				currentResidue.bonds = append(currentResidue.bonds, &bonds{
					atoms: []string{parts[0], parts[1]},
					para:  parts[2],
				})
			}
		case "angles":
			parts := strings.Fields(line)
			if len(parts) >= 4 {
				currentResidue.angles = append(currentResidue.angles, &angles{
					atoms:       []string{parts[0], parts[1], parts[2]},
					gromos_type: parts[3],
				})
			}
		case "dihedrals":
			parts := strings.Fields(line)
			if len(parts) >= 5 {
				currentResidue.dihedrals = append(currentResidue.dihedrals, &dihedrals{
					atoms:       []string{parts[0], parts[1], parts[2], parts[3]},
					gromos_type: parts[4],
				})
			}
		}
		residues[currentResidue.name] = *currentResidue
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return residues, nil
}

// ///////////////
// ////These function are used for read parameter for charge
// ///////////////
// ****highest level function****
func parseChargeFile(filename string) (map[string]map[string]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	chargeData := make(map[string]map[string]float64)
	scanner := bufio.NewScanner(file)
	var currentResidue string

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// Skip empty lines
		if line == "" {
			continue
		}

		// Check for residue header lines like "[ ALA ]"
		if strings.HasPrefix(line, "[") && strings.HasSuffix(line, "]") {
			currentResidue = strings.TrimSpace(line[1 : len(line)-1])
			chargeData[currentResidue] = make(map[string]float64)
			continue
		}

		// Parse atom data lines
		fields := strings.Fields(line)

		// Ensure that we have at least three columns
		if len(fields) < 3 {
			return nil, fmt.Errorf("invalid line format: %s", line)
		}

		atomName := fields[0]
		atomChargeStr := fields[2]

		// Parse atom charge
		atomCharge, err := strconv.ParseFloat(atomChargeStr, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid atom charge '%s' in line: %s", atomChargeStr, line)
		}

		// Store the charge data
		if currentResidue == "" {
			return nil, fmt.Errorf("atom data without residue header: %s", line)
		}
		chargeData[currentResidue][atomName] = atomCharge
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return chargeData, nil
}

func TemporaryPlot(RMSD []float64, time float64) {
	p := plot.New()

	p.Title.Text = "Trajectory analysis"
	p.X.Label.Text = "Time(fs)"
	p.Y.Label.Text = "RMSD"

	points := make(plotter.XYs, len(RMSD))
	for i := range points {
		points[i].X = float64(i) * time
		points[i].Y = RMSD[i]
	}

	s, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}

	p.Add(s)

	if err := p.Save(4*vg.Inch, 4*vg.Inch, "visualization/scatter.png"); err != nil {
		panic(err)
	}

}

// function WriteProteinToPDB takes a protein structure as input
// return a pdb file
func WriteProteinToPDB(protein *Protein, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// Write the title section
	if _, err := writer.WriteString(fmt.Sprintf("HEADER    %s\n", protein.Name)); err != nil {
		return err
	}

	atomIndex := 1
	for _, residue := range protein.Residue {
		for _, atom := range residue.Atoms {
			// Format atom data according to the PDB file format
			_, err := writer.WriteString(fmt.Sprintf(
				"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %-2s\n",
				atomIndex,                                         // Atom serial number
				atom.element,                                      // Atom name
				residue.Name,                                      // Residue name
				residue.ChainID,                                   // Chain identifier
				residue.ID,                                        // Residue sequence number
				atom.position.x, atom.position.y, atom.position.z, // Atom coordinates
				atom.element, // Element symbol
			))
			if err != nil {
				return err
			}
			atomIndex++
		}
	}

	// Write the termination line
	if _, err := writer.WriteString("END\n"); err != nil {
		return err
	}

	return writer.Flush()
}

// writeRMSD writes a slice of float64 values to a CSV file.
func writeRMSD(slice []float64) error {
	// Open the file for writing
	outFile, err := os.Create("result/RMSD.csv")
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	defer outFile.Close()

	writer := csv.NewWriter(outFile)
	defer writer.Flush()

	for i, value := range slice {
		err := writer.Write([]string{strconv.Itoa(i), strconv.FormatFloat(value, 'f', 2, 64)})
		if err != nil {
			log.Fatalf("Failed to write to output file: %v", err)
		}
	}

	return nil
}
