package main

import "math"

func logFactorial(n uint64) float64 {
	if n == 0 || n == 1 {
		return float64(0.0)
	}

	log_fact := float64(0.0)
	for i := uint64(1); i <= n; i++ {
		log_fact += math.Log(float64(i))
	}

	return log_fact
}

func LogMultinomial(counts []uint64, probabilities []float64) float64 {
	total_count := uint64(0)
	mult_log_prob := float64(0.0)

	fact_denum := float64(0.0)

	for i := 0; i < len(counts); i++ {
		// log(x^0) = log(1) = 0
		// log(fact(0)) = 0
		if counts[i] == 0 {
			continue
		}


		// if counts > 0 and prob == 0.0, entire prob is 0
		if probabilities[i] == 0.0 {
			return math.Log(float64(0.0000000000001))
		}

		
		total_count += counts[i]
		fact_denum += logFactorial(counts[i])
		mult_log_prob += float64(counts[i]) * math.Log(probabilities[i])
	}	

	
	fact_num := logFactorial(total_count)

	return fact_num + mult_log_prob - fact_denum
}