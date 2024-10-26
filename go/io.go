package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Function to parse a PDB line based on spaces
func parsePDBLine(line string) (Atom, string, error) {
	fields := strings.Fields(line)
	var atom Atom

	// Parse coordinates
	x, err := strconv.ParseFloat(fields[6], 64)
	if err != nil {
		return Atom{}, "", fmt.Errorf("error parsing x position: %v", err)
	}
	y, err := strconv.ParseFloat(fields[7], 64)
	if err != nil {
		return Atom{}, "", fmt.Errorf("error parsing y position: %v", err)
	}
	z, err := strconv.ParseFloat(fields[8], 64)
	if err != nil {
		return Atom{}, "", fmt.Errorf("error parsing z position: %v", err)
	}

	// Parse element symbol
	element := fields[2]

	// pass value to atom object
	atom.Position.X = x
	atom.Position.Y = y
	atom.Position.Z = z
	atom.Element = element

	// Extract residue name
	residueName := fields[3]

	return atom, residueName, nil
}

func readProteinFromFile(fileName string) (Protein, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return Protein{}, err
	}
	defer file.Close()

	var residues []*Residue
	var currentResidue *Residue
	var protein Protein

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()

		// Process lines that start with "ATOM"
		if strings.HasPrefix(line, "ATOM") {
			atom, residueName, err := parsePDBLine(line)
			if err != nil {
				return Protein{}, err
			}

			// If it's a new residue or first row, start a new Residue object
			if currentResidue == nil || currentResidue.Name != residueName {
				// Add the previous residue to the list if it exists
				if currentResidue != nil {
					residues = append(residues, currentResidue)
				}

				// otherwise,Create a new Residue object
				currentResidue = &Residue{
					Name:  residueName,
					Atoms: []*Atom{&atom},
				}
			} else {
				// If it's the same residue, add the atom to the current residue
				currentResidue.Atoms = append(currentResidue.Atoms, &atom)
			}
		}
	}

	// Add the last residue to the list
	if currentResidue != nil {
		residues = append(residues, currentResidue)
	}

	if err := scanner.Err(); err != nil {
		return Protein{}, err
	}

	// Set the residues in the protein
	protein.Residue = residues

	return protein, nil
}
