package main

type ObservedData struct {
	n_loci int

	locus_counts []*LocusCounts
	locus_frequencies []*LocusFrequencies
}

type LocusCounts struct {
	n_populations int
	n_haplotypes int

	counts [][]uint64
	individuals []uint64
}

type LocusFrequencies struct {
	n_populations int
	n_haplotypes int

	frequencies [][]float64
}

func (locusFreq *LocusFrequencies) Copy() *LocusFrequencies {
	copy_locus := make([][]float64, locusFreq.n_populations)
	for pop_idx := 0; pop_idx < locusFreq.n_populations; pop_idx++ {
		copy_haplo := make([]float64, locusFreq.n_haplotypes)
		copy(copy_haplo, locusFreq.frequencies[pop_idx])
		copy_locus[pop_idx] = copy_haplo
	}

	new_locus_freq := LocusFrequencies{n_populations: locusFreq.n_populations, 
		n_haplotypes: locusFreq.n_haplotypes, frequencies: copy_locus}

	return &new_locus_freq
}

func NewLocusCounts(locus_counts [][]uint64) *LocusCounts {
	n_populations := len(locus_counts)
	n_haplotypes := len(locus_counts[0])

	individuals := make([]uint64, n_populations)
	for i := 0; i < n_populations; i++ {
		for j := 0; j < n_haplotypes; j++ {
			individuals[i] += locus_counts[i][j]
		}
	}

	counts := LocusCounts{n_populations: n_populations, n_haplotypes: n_haplotypes, counts: locus_counts,
		individuals: individuals}

	return &counts
}


func NewObservedData(loci_counts [][][]uint64) *ObservedData {
	n_loci := len(loci_counts)

	observed_counts := make([]*LocusCounts, n_loci)
	observed_freq := make([]*LocusFrequencies, n_loci)

	for i := 0; i < n_loci; i++ {
		locus_counts := NewLocusCounts(loci_counts[i])
		locus_freq := NewLocusFrequenciesFromCounts(locus_counts)
		observed_counts[i] = locus_counts
		observed_freq[i] = locus_freq
	}

	observed := ObservedData{n_loci: n_loci, locus_counts: observed_counts,
		locus_frequencies: observed_freq}

	return &observed
}

func NewLocusFrequencies(freq [][]float64) *LocusFrequencies {
	n_populations := len((freq))
	n_haplotypes := len((freq)[0])

	locus_freq := LocusFrequencies{n_populations: n_populations, n_haplotypes: n_haplotypes, frequencies: freq}

	return &locus_freq
}

func NewLocusFrequenciesFromCounts(locus_counts *LocusCounts) *LocusFrequencies {
	frequencies := make([][]float64, locus_counts.n_populations)

	for i := 0; i < locus_counts.n_populations; i++ {
		haplo_freq := make([]float64, locus_counts.n_haplotypes)

		for j := 0; j < locus_counts.n_haplotypes; j++ {
			haplo_freq[j] = float64(locus_counts.counts[i][j]) / float64(locus_counts.individuals[i])
		}

		frequencies[i] = haplo_freq
	}

	return NewLocusFrequencies(frequencies)
}
