package main

import "fmt"

func main() {
	fmt.Printf("Hello, world.\n")

	occur, _ := ReadOccurrences("../apple_bamovacomp_multinomial.txt")

	fmt.Printf("Loci: %d, Populations: %d, Haplotypes: %d\n", occur.n_loci, occur.n_populations, occur.n_haplotypes)
}