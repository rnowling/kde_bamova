package main

type ObservedData struct {
	n_loci int
	n_populations int
	n_haplotypes int

	observed_counts [][][]uint64

	observed_frequencies *PopulationFrequencies

	individual_counts [][]uint64
}

type PopulationFrequencies struct {
	n_loci int
	n_populations int
	n_haplotypes int

	population_frequencies [][][]float64
}

func NewObservedData(loci_counts *[][][]uint64) *ObservedData {
	n_loci := len(*loci_counts)
	n_populations := len((*loci_counts)[0])
	n_haplotypes := len((*loci_counts)[0][0])

	observed := ObservedData{n_loci: n_loci, n_populations: n_populations, n_haplotypes: n_haplotypes, observed_counts: *loci_counts}

	observed.countIndividuals()
	observed.normalizeCounts()

	return &observed
}

func NewPopulationFrequencies(freq *[][][]float64) *PopulationFrequencies {
	n_loci := len(*freq)
	n_populations := len((*freq)[0])
	n_haplotypes := len((*freq)[0][0])

	pop_freq := PopulationFrequencies{n_loci: n_loci, n_populations: n_populations, n_haplotypes: n_haplotypes, population_frequencies: *freq}

	return &pop_freq
}

func (observed *ObservedData) countIndividuals() {
	individuals := make([][]uint64, observed.n_loci)

	for i := 0; i < observed.n_loci; i++ {
		pop_counts := make([]uint64, observed.n_populations)

		for j := 0; j < observed.n_populations; j++ {
			sum := uint64(0)

			for k := 0; k < observed.n_haplotypes; k++ {
				sum += observed.observed_counts[i][j][k]
			}

			pop_counts[j] = sum
		}

		individuals[i] = pop_counts
	}

	observed.individual_counts = individuals
}

func (observed *ObservedData) normalizeCounts() {
	observed_frequencies := make([][][]float64, observed.n_loci)

	for i := 0; i < observed.n_loci; i++ {
		loci := observed.observed_counts[i]
		pop_freq := make([][]float64, observed.n_populations)

		for j := 0; j < observed.n_populations; j++ {
			population := loci[j]
			loci_freq := make([]float64, observed.n_haplotypes)

			for k, haplotype_cnt := range population {
				loci_freq[k] = float64(haplotype_cnt) / float64(observed.individual_counts[i][j])
			}

			pop_freq[j] = loci_freq
		}

		observed_frequencies[i] = pop_freq
	}

	observed.observed_frequencies = NewPopulationFrequencies(&observed_frequencies)
}
