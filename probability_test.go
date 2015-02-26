package main

import "math"
import "testing"

func TestLogMultinomial(t *testing.T) {
	counts := []uint64{1, 1, 1, 1, 1}
	prob := []float64{0.0, 0.0, 0.0, 0.0, 0.0}

	result := LogMultinomial(counts, prob)
	expected := math.Log(float64(0.0000000000001))

	if result != expected {
		t.Errorf("LogMultinomial returned incorrect result when given positive counts and all zero probabilities. Found: %f, Expected %f", result, expected)
	}

	prob = []float64{0.25, 0.25, 0.25, 0.25, 0.0}
	result = LogMultinomial(counts, prob)

	if result != expected {
		t.Errorf("LogMultinomial returned incorrect result when given positive counts and a zero probability. Found: %f, Expected %f", result, expected)
	}

	counts = []uint64{1, 1, 1, 1, 0}
	prob = []float64{0.2, 0.2, 0.2, 0.2, 0.2}
	result = LogMultinomial(counts, prob)
	expected = float64(-3.2596978193884554)

	if (result - expected) > 1e-10 {
		t.Errorf("LogMultinomial returned incorrect result. Found: %f, Expected %f", result, expected)
	}

	counts = []uint64{1, 2, 3, 4, 5}
	prob = []float64{0.06666666666666667, 0.13333333333333333, 0.2, 0.26666666666666666, 0.3333333333333333}
	result = LogMultinomial(counts, prob)
	expected = float64(-4.8974356218359674)

	if (result - expected) > 1e-10 {
		t.Errorf("LogMultinomial returned incorrect result. Found: %f, Expected %f", result, expected)
	}
}