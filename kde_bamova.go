package main

import "fmt"

func main() {
	fmt.Printf("Hello, world.\n")

	observed, _ := ReadOccurrences("../apple_bamovacomp_multinomial.txt")

	fmt.Printf("Loci: %d, Populations: %d, Haplotypes: %d\n", observed.n_loci, observed.n_populations, observed.n_haplotypes)
}