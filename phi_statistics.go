package main

func ssdwpfunc(popFreq *PopulationFrequencies, observed *ObservedData, locus_idx int) float64 {
	ssdwp := float64(0.0)

	for pop_idx := 0; pop_idx < observed.n_populations; pop_idx++ {
		gene_copies := observed.individual_counts[locus_idx][pop_idx] * uint64(2) // DIPLOID

		pop_ssd := float64(0.0)
		for i := 0; i < observed.n_haplotypes - 1; i++ {
			for j := i + 1; j < observed.n_haplotypes; j++ {
				pop_ssd += popFreq.population_frequencies[locus_idx][pop_idx][i] * popFreq.population_frequencies[locus_idx][pop_idx][j]
			}
		}

		ssdwp += float64(gene_copies) * pop_ssd
	}

	return ssdwp
}

func ssdtotalfunc(popFreq *PopulationFrequencies, observed *ObservedData, locus_gene_copies uint64, locus_idx int) float64 {
	total_freq := make([]float64, observed.n_haplotypes)
	for haplo_idx := 0; haplo_idx < observed.n_haplotypes; haplo_idx++ {
		top := float64(0.0)
		for pop_idx := 0; pop_idx < observed.n_populations; pop_idx++ {
			pop_gene_copies := observed.individual_counts[locus_idx][pop_idx] * uint64(2) // DIPLOID 
			top += float64(pop_gene_copies) * popFreq.population_frequencies[locus_idx][pop_idx][haplo_idx]
		}
		total_freq[haplo_idx] = top / float64(locus_gene_copies)
	}

	ssd := float64(0.0)
	for haplo1 := 0; haplo1 < observed.n_haplotypes - 1; haplo1++ {
		for haplo2 := haplo1 + 1; haplo2 < observed.n_haplotypes; haplo2++ {
			ssd += float64(locus_gene_copies) * total_freq[haplo1] * total_freq[haplo2]
		}
	}

	return ssd
}

func CalculateLocusPhi(popFreq *PopulationFrequencies, observed *ObservedData, locus_idx int) float64 {
	locus_gene_copies := uint64(0)
	for pop_idx := 0; pop_idx < observed.n_populations; pop_idx++ {
		locus_gene_copies += uint64(2) * observed.individual_counts[locus_idx][pop_idx] // DIPLOID
	}

	ssdwp := ssdwpfunc(popFreq, observed, locus_idx)

	ssdtotal := ssdtotalfunc(popFreq, observed, locus_gene_copies, locus_idx)

	dfb := popFreq.n_populations - 1

	dfw := locus_gene_copies - uint64(popFreq.n_populations)

	if dfw == 0 {
		return float64(0.0)
	}

	msdwp := ssdwp / float64(dfw)

	msdap := (ssdtotal - ssdwp) / float64(dfb)

	varAP := (msdap - msdwp) / (float64(locus_gene_copies) / float64(popFreq.n_populations))

	if (varAP + msdwp) == 0.0 {
		return float64(0.0)
	}

	phi_st := varAP/(msdwp + varAP)

	return phi_st
}

func CalculatePhis(popFreq *PopulationFrequencies, observed *ObservedData) *[]float64 {
	phi_st := make([]float64, popFreq.n_loci)

	for i := 0; i < popFreq.n_loci; i++ {
		phi_st[i] = CalculateLocusPhi(popFreq, observed, i)
	}

	return &phi_st
}