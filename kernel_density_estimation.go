package main

import "math"

type KernelDensityEstimate struct {
	model_values []float64
	bandwidth float64
}

func NewKDE(model_values *[]float64, bandwidth float64) *KernelDensityEstimate {
	model_copy := make([]float64, len(*model_values))
	copy(model_copy, *model_values)

	kde := KernelDensityEstimate{model_values: model_copy, bandwidth: bandwidth}

	return &kde
}

func gaussKernel(x float64) float64 {
	num := math.Exp(float64(-0.5) * x * x)
	denom := math.Sqrt(float64(2.0) * math.Pi)

	return num / denom
}

func (kde *KernelDensityEstimate) LogProb(value float64) float64 {
	log_prob := float64(0.0)

	nh := float64(len(kde.model_values)) * kde.bandwidth

	for i := 0; i < len(kde.model_values); i++ {
		x := (value - kde.model_values[i]) / kde.bandwidth
		log_prob += gaussKernel(x) / (nh)
	}

	return log_prob
}