package main

import (
	"fmt"
	"math"
)

func AddSimpleBondConstraints(p *Protein, bondParameter parameterDatabase, residueParameterValue map[string]residueParameter) {
	tolerance := 2

	coverage := false
	for !coverage {
		coverage = true

		for q, residue := range p.Residue {

			// Calculate bondstretch energy
			// range over each bond

			for _, bondPairs := range residueParameterValue[residue.Name].bonds {
				for i := 0; i < len(residue.Atoms)-1; i++ {
					atom1 := residue.Atoms[i]
					if atom1.element[0] == (*bondPairs).atoms[0][0] {
						for j := i + 1; j < len(residue.Atoms); j++ {
							if residue.Atoms[j].element[0] == (*bondPairs).atoms[1][0] {

								atom2 := residue.Atoms[j]
								parameterList := SearchParameter(2, bondParameter, atom1, atom2)
								if len(parameterList) != 1 {
									r0 := parameterList[0]

									r := Distance(atom1.position, atom2.position)
									if r == 0 {
										continue
									}
									maximalIteration := 5
									remark := 0
									for r > float64(tolerance) && remark < maximalIteration {
										correction := (r - r0*10) / 100000
										remark += 1
										var tempAtom Atom
										tempAtom.position.x = atom1.position.x
										tempAtom.position.y = atom1.position.y
										tempAtom.position.z = atom1.position.z
										if math.IsNaN((tempAtom.position.x - atom2.position.x) * correction) {
											break
										}

										if math.IsNaN((tempAtom.position.y - atom2.position.y) * correction) {
											break
										}
										if math.IsNaN((tempAtom.position.z - atom2.position.z) * correction) {
											break
										}

										if math.IsNaN((atom2.position.x - atom1.position.x) * correction) {
											break
										}

										if math.IsNaN((atom2.position.y - atom1.position.y) * correction) {
											break
										}
										if math.IsNaN((atom2.position.z - atom1.position.z) * correction) {
											break
										}

										if atom1.index == 233 {
											fmt.Println(atom1.position)
										}
										if atom2.index == 233 {
											fmt.Println(atom1.position)
										}
										p.Residue[q].Atoms[i].position.x += (atom2.position.x - atom1.position.x) * correction
										p.Residue[q].Atoms[i].position.y += (atom2.position.y - atom1.position.y) * correction
										p.Residue[q].Atoms[i].position.z += (atom2.position.z - atom1.position.z) * correction

										p.Residue[q].Atoms[j].position.x += (tempAtom.position.x - atom2.position.x) * correction
										p.Residue[q].Atoms[j].position.y += (tempAtom.position.y - atom2.position.y) * correction
										p.Residue[q].Atoms[j].position.z += (tempAtom.position.z - atom2.position.z) * correction
										coverage = false
										if atom1.index == 233 {
											fmt.Println(atom1.position)
										}
										if atom2.index == 233 {
											fmt.Println(atom1.position)
										}
									}

								}

							}
						}
					}
				}
			}
		}

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

							parameterList := SearchParameter(2, bondParameter, atom1, atom2)
							if len(parameterList) != 1 {
								r0 := parameterList[0]
								r := Distance(atom1.position, atom2.position)
								if r == 0 {
									continue
								}
								maximalIteration := 5
								remark := 0
								for r > float64(tolerance) || remark < maximalIteration {
									remark += 1
									correction := (r - r0*10) / 100000
									if math.IsNaN((atom1.position.x - atom2.position.x) * correction) {
										break
									}

									if math.IsNaN((atom1.position.y - atom2.position.y) * correction) {
										break
									}
									if math.IsNaN((atom1.position.z - atom2.position.z) * correction) {
										break
									}

									if math.IsNaN((atom2.position.x - atom1.position.x) * correction) {
										break
									}

									if math.IsNaN((atom2.position.y - atom1.position.y) * correction) {
										break
									}
									if math.IsNaN((atom2.position.z - atom1.position.z) * correction) {
										break
									}
									p.Residue[m].Atoms[n].position.x += (atom1.position.x - atom2.position.x) * correction
									p.Residue[m].Atoms[n].position.y += (atom1.position.y - atom2.position.y) * correction
									p.Residue[m].Atoms[n].position.z += (atom1.position.z - atom2.position.z) * correction

									p.Residue[m+1].Atoms[t].position.x += (atom2.position.x - atom1.position.x) * correction
									p.Residue[m+1].Atoms[t].position.y += (atom2.position.y - atom1.position.y) * correction
									p.Residue[m+1].Atoms[t].position.z += (atom2.position.z - atom1.position.z) * correction
									coverage = false
								}

							}
						}
					}
				}
				break
			}
		}

	}

}
