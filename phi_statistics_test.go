package main

import "math"
import "testing"

func TestCalculateLocusPhi(t *testing.T) {
	expected_phi_st := 0.431818181818
	tol := 1e-5
	pop1_indiv := []uint64{0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	pop2_indiv := []uint64{0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	counts := [][]uint64{pop1_indiv, pop2_indiv}

	locus_counts := NewLocusCounts(counts)
	locus_frequencies := NewLocusFrequenciesFromCounts(locus_counts)

	phi_st := CalculateLocusPhi(locus_frequencies, locus_counts)

	if math.Abs(expected_phi_st - phi_st) > tol {
			t.Errorf("Incorrect phi_st found. Found: %f, Expected %f", phi_st, expected_phi_st)
	}
}