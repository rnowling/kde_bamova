package main

import "fmt"
import "flag"
import "os"

var seed int64
var steps int64
var inputPath string
var phiOutputPath string
var freqOutputPath string
var bandwidth float64
var outputPeriod int64

func parseArgs() {
	flag.Int64Var(&seed, "seed", -1, "Seed for RNG")
	flag.Int64Var(&steps, "steps", -1, "Number of steps to run")
	flag.StringVar(&inputPath, "input", "", "Input file")
	flag.StringVar(&phiOutputPath, "phiOut", "", "Phi output file")
	flag.StringVar(&freqOutputPath, "freqOut", "", "Frequency output file")
	flag.Float64Var(&bandwidth, "bandwidth", 0.01, "Bandwidth for KDE")
	flag.Int64Var(&outputPeriod, "outputPeriod", 10, "Output period")

	flag.Parse()

	if inputPath == "" {
		fmt.Printf("Path to input file not specified.\n")
		os.Exit(1)
	} 

	if steps == -1 {
		fmt.Printf("Number of steps not specified.\n")
		os.Exit(1)
	}
}

func main() {
	parseArgs()

	observed, _ := ReadOccurrences(inputPath)

	fmt.Printf("Loci: %d, Populations: %d, Haplotypes: %d\n", observed.n_loci, observed.locus_counts[0].n_populations, observed.locus_counts[0].n_haplotypes)

	sampler := NewSampler(observed, bandwidth)

	var phiOutputFile *os.File
	var freqOutputFile *os.File
	var err error

	if phiOutputPath != "" {
		phiOutputFile, err = os.Create(phiOutputPath)

		if err != nil {
			fmt.Printf("Failed to open phi output file for writing.\n")
			os.Exit(1)
		}
	}

	if freqOutputPath != "" {
		freqOutputFile, err = os.Create(freqOutputPath)

		if err != nil {
			fmt.Printf("Failed to open frequency output file for writing.\n")
			os.Exit(1)
		}
	}

	for step := int64(0); step < steps; step++ {
		sample := sampler.Sample()
		fmt.Printf("Step %d, Log Prob %f\n", step, sample.log_probability)

		if step % outputPeriod == 0 {
			if phiOutputPath != "" {
				WriteLocusPhiValues(phiOutputFile, sample, step)
			}

			if freqOutputPath != "" {
				WriteFrequencies(freqOutputFile, sample, step)
			}
		}
	}

	if phiOutputPath != "" {
		phiOutputFile.Close()
	}

	if freqOutputPath != "" {
		freqOutputFile.Close()
	}
	
}