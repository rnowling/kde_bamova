package main

import "fmt"
import "os"
import "bufio"
import "strings"
import "strconv"

func ReadOccurrences(flname string) (*Occurrences, error) {
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

	occur := Occurrences{n_loci: len(loci), n_populations: len(loci[0]), n_haplotypes: len(loci[0][0]), observed_counts: &loci}

	return &occur, err
}