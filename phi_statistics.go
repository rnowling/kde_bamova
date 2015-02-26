package main

func ssdwpfunc(locusFreq *LocusFrequencies, locusCounts *LocusCounts) float64 {
	ssdwp := float64(0.0)

	for pop_idx := 0; pop_idx < locusFreq.n_populations; pop_idx++ {
		gene_copies := locusCounts.individuals[pop_idx] * uint64(2) // DIPLOID

		pop_ssd := float64(0.0)
		for i := 0; i < locusFreq.n_haplotypes - 1; i++ {
			for j := i + 1; j < locusFreq.n_haplotypes; j++ {
				pop_ssd += locusFreq.frequencies[pop_idx][i] * locusFreq.frequencies[pop_idx][j]
			}
		}

		ssdwp += float64(gene_copies) * pop_ssd
	}

	return ssdwp
}

func ssdtotalfunc(locusFreq *LocusFrequencies, locusCounts *LocusCounts, locus_gene_copies uint64) float64 {
	total_freq := make([]float64, locusFreq.n_haplotypes)
	for haplo_idx := 0; haplo_idx < locusFreq.n_haplotypes; haplo_idx++ {
		top := float64(0.0)
		for pop_idx := 0; pop_idx < locusFreq.n_populations; pop_idx++ {
			pop_gene_copies := locusCounts.individuals[pop_idx] * uint64(2) // DIPLOID 
			top += float64(pop_gene_copies) * locusFreq.frequencies[pop_idx][haplo_idx]
		}
		total_freq[haplo_idx] = top / float64(locus_gene_copies)
	}

	ssd := float64(0.0)
	for haplo1 := 0; haplo1 < locusFreq.n_haplotypes - 1; haplo1++ {
		for haplo2 := haplo1 + 1; haplo2 < locusFreq.n_haplotypes; haplo2++ {
			ssd += float64(locus_gene_copies) * total_freq[haplo1] * total_freq[haplo2]
		}
	}

	return ssd
}

func CalculateLocusPhi(locusFreq *LocusFrequencies, locusCounts *LocusCounts) float64 {
	locus_gene_copies := uint64(0)
	for pop_idx := 0; pop_idx < locusCounts.n_populations; pop_idx++ {
		locus_gene_copies += uint64(2) * locusCounts.individuals[pop_idx] // DIPLOID
	}

	ssdwp := ssdwpfunc(locusFreq, locusCounts)

	ssdtotal := ssdtotalfunc(locusFreq, locusCounts, locus_gene_copies)

	dfb := locusFreq.n_populations - 1

	dfw := locus_gene_copies - uint64(locusFreq.n_populations)

	if dfw == 0 {
		return float64(0.0)
	}

	msdwp := ssdwp / float64(dfw)

	msdap := (ssdtotal - ssdwp) / float64(dfb)

	varAP := (msdap - msdwp) / (float64(locus_gene_copies) / float64(locusFreq.n_populations))

	if (varAP + msdwp) == 0.0 {
		return float64(0.0)
	}

	phi_st := varAP/(msdwp + varAP)

	return phi_st
}

func CalculatePhis(locusFreq []*LocusFrequencies, locusCounts []*LocusCounts) *[]float64 {
	n_loci := len(locusFreq)
	phi_st := make([]float64, n_loci)

	for i := 0; i < n_loci; i++ {
		phi_st[i] = CalculateLocusPhi(locusFreq[i], locusCounts[i])
	}

	return &phi_st
}