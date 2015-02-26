package main

import "fmt"

func main() {
	fmt.Printf("Hello, world.\n")

	observed, _ := ReadOccurrences("../apple_bamovacomp_multinomial_3000.txt")

	fmt.Printf("Loci: %d, Populations: %d, Haplotypes: %d\n", observed.n_loci, observed.locus_counts[0].n_populations, observed.locus_counts[0].n_haplotypes)

	sampler := NewSampler(observed, float64(0.01))
	for step := 0; step < 10; step++ {
		sample := sampler.Sample()
		fmt.Printf("Step %d, Log Prob %f\n", step, sample.log_probability)
	}
	
}