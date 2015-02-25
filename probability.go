package main

import "math"
import "math/rand"

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

func SampleStandardNormal() float64 {
	u1 := rand.Float64()
	u2 := rand.Float64()

	z0 := math.Sqrt(-2.0 * math.Log(u1)) * math.Cos(2.0 * math.Pi * u2)
	//z1 := math.Sqrt(-2.0 * math.Log(u1)) * math.Sin(2.0 * math.Pi * u2)

	return z0
}

func SampleGamma(shape float64, scale float64) float64 {
	if shape < 1 {
		return sampleGammaSmallShape(shape, scale)
	} else {
		return sampleGammaLargeShape(shape, scale)
	}

}

func sampleGammaSmallShape(shape float64, scale float64) float64 {
	for true {
		u := rand.Float64()
		bGS := 1.0 + shape / math.E
		p := bGS * u

		if (p <= 1.0) {
			x := math.Pow(p, 1.0 / shape)
			u2 := rand.Float64()

			if u2 > math.Exp(-x) {
				continue
			} else {
				return scale * x
			}
		} else {
			x := -1.0 * math.Log((bGS - p) / shape)
			u2 := rand.Float64()

			if u2 > math.Pow(x, shape - 1.0) {
				continue
			} else {
				return scale * x
			}
		}
	}

	return -1.0
}

func sampleGammaLargeShape(shape float64, scale float64) float64 {
	d := shape - 0.333333333333333333
	c := 1.0 / (3.0 * math.Sqrt(d))

	for true {
		x := rand.Float64() // TODO Gaussian
		base := (1 + c * x)
		v :=  base * base * base

		if v <= 0 {
			continue
		} 

		x2 := x * x
		u := rand.Float64()

		if u < 1.0 - 0.0331 * x2 * x2 {
			return scale * d * v
		}

		if math.Log(u) < 0.5 * x2 + d * (1 - v + math.Log(v)) {
			return scale * d * v
		}
	}

	return -1.0
}