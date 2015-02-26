package main

import "math"
import "testing"

func TestNewLocusFrequenciesFromCounts(t *testing.T) {
	tol := 1e-5
	pop1_indiv := []uint64{0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	pop2_indiv := []uint64{0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	counts := [][]uint64{pop1_indiv, pop2_indiv}

	locus_counts := NewLocusCounts(counts)
	locus_frequencies := NewLocusFrequenciesFromCounts(locus_counts)

	expected_freq1 := []float64{0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	expected_freq2 := []float64{0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	expected_freq := [][]float64{expected_freq1, expected_freq2}

	for pop_idx := 0; pop_idx < 2; pop_idx++ {
		for haplo_idx := 0; haplo_idx < 16; haplo_idx++ {
			if math.Abs(expected_freq[pop_idx][haplo_idx] - locus_frequencies.frequencies[pop_idx][haplo_idx]) > tol {
				t.Errorf("Incorrect locus frequency found. Found: %f, Expected %f", locus_frequencies.frequencies[pop_idx][haplo_idx], 
					expected_freq[pop_idx][haplo_idx])
			}
		}
	}

	
}