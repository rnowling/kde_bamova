package main

import "fmt"
import "math"
import "math/rand"

type Model int

const (
	KDE_MODEL Model = iota
	UNIFORM_MODEL Model = iota
	SINGLE_MODEL Model = iota
	KDE_ONLY_MODEL Model = iota
)

type Sampler struct {
	n_loci int
	// thread pool
	sampler_channels []chan *LocusSample
}

type LocusSimulationParameters struct {
	kde *KernelDensityEstimate // nil if model is not KDE_MODEL
	counts *LocusCounts
	model Model
}

type LocusSample struct {
	n_populations int
	n_haplotypes int

	freq *LocusFrequencies
	phi_st float64
	log_probability float64
}

type Sample struct {
	n_loci int

	locus_frequencies []*LocusFrequencies
	locus_phi_values []float64
	log_probability float64
}



func NewLocusSample(freq *LocusFrequencies, phi_st float64, log_probability float64) *LocusSample {
	sample := LocusSample{n_populations: freq.n_populations, n_haplotypes:
		freq.n_haplotypes, freq: freq, phi_st: phi_st, log_probability: log_probability}

	return &sample
}

func NewLocusSimulationParameters(model Model, kde *KernelDensityEstimate, counts *LocusCounts) *LocusSimulationParameters {
	params := LocusSimulationParameters{model: model, kde: kde, counts: counts}

	return &params
}

func NewSampler(model Model, observed *ObservedData, bandwidth float64) *Sampler {
	n_loci := observed.n_loci

	fmt.Printf("Calculating phis\n")
	locus_phi_values := CalculatePhis(observed.locus_frequencies, observed.locus_counts)

	var kde *KernelDensityEstimate = nil
	if model == KDE_MODEL || model == KDE_ONLY_MODEL {
		fmt.Printf("Training KDE\n")
		kde = NewKDE(*locus_phi_values, bandwidth)
	}

	fmt.Printf("Creating go routines\n")
	sampler_channels := make([]chan *LocusSample, n_loci)
	for i := 0; i < n_loci; i++ {
		params := NewLocusSimulationParameters(model, kde, observed.locus_counts[i])
		sampler_channels[i] = make(chan *LocusSample)

		go sampleLocus(params, observed.locus_frequencies[i], sampler_channels[i])
	}

	sampler := Sampler{n_loci: n_loci, sampler_channels: sampler_channels}

	return &sampler
}

func (sampler *Sampler) Sample() *Sample {
	fmt.Printf("Executing Sample process\n")
	locus_frequencies := make([]*LocusFrequencies, sampler.n_loci)
	locus_phi_values := make([]float64, sampler.n_loci)
	log_probability := float64(0.0)

	for i := 0; i < sampler.n_loci; i++ {
		locus_sample := <- sampler.sampler_channels[i]
		locus_phi_values[i] = locus_sample.phi_st
		log_probability += locus_sample.log_probability
		locus_frequencies[i] = locus_sample.freq
	}

	sample := Sample{n_loci: sampler.n_loci, locus_frequencies: locus_frequencies,
		locus_phi_values: locus_phi_values, log_probability: log_probability}

	return &sample
}

func sampleLocusFreq(prev_freq *LocusFrequencies, pop_idx int, params *LocusSimulationParameters) *LocusFrequencies {
	alphas := make([]float64, prev_freq.n_haplotypes)

	DWEIGHT := float64(1.0)
	DADD := float64(1.0)
	SMALLNUM := float64(0.0000000000001)

	for i := 0; i < prev_freq.n_haplotypes; i++ {
		alphas[i] = DWEIGHT * prev_freq.frequencies[pop_idx][i] + DADD + SMALLNUM
	}

	new_freq := SampleDirichlet(alphas)

	new_locus_freq := prev_freq.Copy()
	new_locus_freq.frequencies[pop_idx] = new_freq

	return new_locus_freq
}

func logProbability(freq *LocusFrequencies, params *LocusSimulationParameters) (float64, float64) {
	phi_st := CalculateLocusPhi(freq, params.counts)
	log_prob := float64(0.0)

	if params.model == KDE_MODEL || params.model == KDE_ONLY_MODEL {
		log_prob += params.kde.LogProb(phi_st)
	}

	if params.model != KDE_ONLY_MODEL {
		for pop_idx := 0; pop_idx < params.counts.n_populations; pop_idx++ {
			log_prob += LogMultinomial(params.counts.counts[pop_idx],
				freq.frequencies[pop_idx])
		}
	}

	return log_prob, phi_st
}

func sampleLocus(params *LocusSimulationParameters, initialFreq *LocusFrequencies, output chan *LocusSample) {
	locus_freq := sampleLocusFreq(initialFreq, 0, params)
	for pop_idx := 1; pop_idx < initialFreq.n_populations; pop_idx++ {
		locus_freq = sampleLocusFreq(locus_freq, pop_idx, params)
	}

	log_prob, phi_st := logProbability(locus_freq, params)

	for true {
		sample := NewLocusSample(locus_freq, phi_st, log_prob)
		output <- sample

		for pop_idx := 0; pop_idx < initialFreq.n_populations; pop_idx++ {
			proposed_freq := sampleLocusFreq(locus_freq, pop_idx, params)

			proposed_log_prob, proposed_phi_st := logProbability(proposed_freq, params)

			log_r := math.Log(rand.Float64())
			log_ratio := proposed_log_prob - log_prob

			if proposed_log_prob >= log_prob || log_r < log_ratio {
				phi_st = proposed_phi_st
				log_prob = proposed_log_prob
				locus_freq = proposed_freq
			}
		}
	}
}