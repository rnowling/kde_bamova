package main

type Occurrences struct {
	n_loci int
	n_populations int
	n_haplotypes int

	observed_counts *[][][]uint64
}