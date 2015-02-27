package main

import "fmt"
import "os"
import "bufio"
import "strings"
import "strconv"

func WriteFrequencies(fl *os.File, sample *Sample, step int64) {
	for locus_idx := 0; locus_idx < sample.n_loci; locus_idx++ {
		
		locus_frequencies := sample.locus_frequencies[locus_idx]
		
		for pop_idx := 0; pop_idx < locus_frequencies.n_populations; pop_idx++ {
			
			fl.WriteString(fmt.Sprintf("%d,%d,%d,%d", locus_idx, pop_idx, step, locus_frequencies.n_haplotypes))
			
			for haplo_idx := 0; haplo_idx < locus_frequencies.n_haplotypes; haplo_idx++ {
				fl.WriteString(fmt.Sprintf(",%f", locus_frequencies.frequencies[haplo_idx]))
			}

			fl.WriteString("\n")
		}
	}
}

func WriteLocusPhiValues(fl *os.File, sample *Sample, step int64) {
	for locus_idx := 0; locus_idx < sample.n_loci; locus_idx++ {
		fl.WriteString(fmt.Sprintf("%d,%d,%f\n", step, locus_idx, sample.locus_phi_values[locus_idx]))
	}
}

func ReadOccurrences(flname string) (*ObservedData, error) {
	file, err := os.Open(flname)

	if err != nil {
		return nil, err
	}
	defer file.Close()

	populations := make([][]uint64, 0)

	loci := make([][][]uint64, 0)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		ln := scanner.Text()

		// new locus
		if strings.Contains(ln, "Marker") {
			if len(populations) != 0 {
				loci = append(loci, populations)
				populations = make([][]uint64, 0)
			}

			continue
		} else {
			cols := strings.Split(ln, " ")

			haplotype_counts := make([]uint64, len(cols) - 2)

			for i, col := range cols[2:] {
				value, _ := strconv.ParseUint(col, 10, 64)

				haplotype_counts[i] = value
			}

			populations = append(populations, haplotype_counts)
		}
	}

	fmt.Printf("Loci: %d, Populations: %d, Haplotypes: %d\n", len(loci), len(loci[0]), len(loci[0][0]))

	observed := NewObservedData(loci)

	return observed, err
}