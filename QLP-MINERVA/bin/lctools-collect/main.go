// Package main contains code for the utility lctools-collect.
// lctools-collect takes the outputs of fiphot, which contains the photometry
// results of every target in one frame, and outputs one file for each target
// containing photometry results from all frames, i.e. light curves.
// It performs this task in-memory efficiently using the model of MapReduce.
// Nevertheless, its performance is bound by the I/O speed of the medium that
// stores its inputs and outputs.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math/bits"
	"math/rand"
	"os"
	"path"
	"runtime"
	"sync"
)

var (
	outputDir = ""
)

// Type Tic represents a TIC ID.
// It is currently just a int, which is sufficient on modern 64-bit arch.
type Tic int

// Type FiphotLine is a data structure passed between map and reduce workers.
// It stores the TIC ID of the target and the exact line of fiphot output.
type FiphotLine struct {
	Tic  Tic
	Data string
}

// mapWorker takes a file from queue, reads its lines, and send each line to the appropriate reduceWorker.
func mapWorker(id int, queue <-chan string, toReduce [](chan *FiphotLine), wg *sync.WaitGroup, hashFunc func(Tic) int) {
	for filename := range queue {
		fmt.Printf("[%v] Mapping %v.\n", id, filename)
		f, err := os.Open(filename)
		if err != nil {
			fmt.Printf("Open %v error: %v\n", filename, err)
			continue
		}
		scanner := bufio.NewScanner(f)
		var tic Tic
		// We want to send FiphotLines from this file to the reduceWorkers only when
		// no error has occured, thus the sending operation atomic on read and scan
		// success.
		var fiphotLines []*FiphotLine
		for scanner.Scan() {
			data := scanner.Text()
			if data == "" || data[0] == '#' {
				continue
			}
			_, err := fmt.Sscanf(data, "%d", &tic)
			if err != nil {
				fmt.Printf("Malformed line in %v: %v\n%v\n", filename, err, data)
				continue
			}
			fiphot := FiphotLine{tic, data}
			fiphotLines = append(fiphotLines, &fiphot)
		}
		if err := scanner.Err(); err != nil {
			fmt.Printf("Error scanning %v: %v\n", filename, err)
		}
		for _, fiphot := range fiphotLines {
			toReduce[hashFunc(fiphot.Tic)] <- fiphot
		}
		f.Close()
	}
	wg.Done()
}

// reduceWorker takes a FiphotLine from queue, groups them by TIC, and write out buffered lines for each TIC.
func reduceWorker(id int, queue <-chan *FiphotLine, wg *sync.WaitGroup) {
	// Buffers lines to be writen, grouped by TIC.
	fiphotMap := make(map[Tic][]string)
	for fiphot := range queue {
		if lines, ok := fiphotMap[fiphot.Tic]; ok {
			fiphotMap[fiphot.Tic] = append(lines, fiphot.Data)
		} else {
			fiphotMap[fiphot.Tic] = []string{fiphot.Data}
		}
	}
	fmt.Printf("[reduce %v]: Receiving channel closed; writing %v light curves to disk.\n", id, len(fiphotMap))
	for tic, lines := range fiphotMap {
		err := writeLines(tic, lines)
		if err != nil {
			fmt.Println(err)
		}
	}
	wg.Done()
}

// Takes a slice of strings and write them to a file whose name is determined by the TIC.
func writeLines(tic Tic, lines []string) error {
	f, err := os.Create(path.Join(outputDir, fmt.Sprintf("%v.rlc", tic)))
	if err != nil {
		return fmt.Errorf("Error creating final phot file for %v: %v", tic, err)
	}
	defer f.Close()
	writer := bufio.NewWriter(f)
	for _, line := range lines {
		_, err = fmt.Fprintln(writer, line)
		if err != nil {
			return fmt.Errorf("Error writing to buffer for %v: %v", tic, err)
		}
	}
	err = writer.Flush()
	if err != nil {
		return fmt.Errorf("Error flushing buffer for %v: %v", tic, err)
	}
	return nil
}

// Creates the hash function for TICs based on the number of workers available.
// A universal hashing function using multiply-add-shift is used.
// The number of workers MUST BE A POWER OF TWO. Bad things happen if it isn't.
func hashTic(numWorkers int) func(Tic) int {
	// The log base-2 of the number of workers
	logN := uint(bits.TrailingZeros(uint(numWorkers)))
	// Parameter a is an odd positive integer modulo 1 << (word size).
	a := (uint(rand.Int()) << 1) + 1
	// Parameter b is a non-negative integer less than numWorkers.
	b := uint(rand.Intn(numWorkers))
	// We want the highest logN bits. Thus, right shift (word size) - logN.
	shiftBits := bits.UintSize - logN
	return func(id Tic) int {
		return int((a*uint(id) + b) >> shiftBits)
	}
}

// Sanity checks number of reduce workers
func setNumReduceWorker(numReduceWorker *int) error {
	if *numReduceWorker <= 0 {
		return fmt.Errorf("Error: number of reduce workers must be positive")
	}
	n := uint(*numReduceWorker)
	// If n is not a power of two, return the largest power of two less than n
	logN := bits.UintSize - 1 - uint(bits.LeadingZeros(n))
	if (n & (n - 1)) != 0 {
		*numReduceWorker = 1 << logN
	}
	return nil
}

// Parse the input file into slice of filenames.
func parseInFile(inFile string) ([]string, error) {
	f, err := os.Open(inFile)
	if err != nil {
		return nil, fmt.Errorf("Error opening input file %v: %v", inFile, err)
	}
	reader := bufio.NewScanner(f)
	var filenames []string
	for reader.Scan() {
		filenames = append(filenames, reader.Text())
	}
	if err := reader.Err(); err != nil {
		return nil, fmt.Errorf("Error scanning input file %v: %v", inFile, err)
	}
	return filenames, nil
}

func main() {
	// Set up and parse flags.
	inFile := flag.String("infile", "phot.in", "File containing filenames of input files.")
	defaultNumWorkers := runtime.GOMAXPROCS(0) / 2
	var numMapWorkers, numReduceWorkers int
	flag.StringVar(&outputDir, "outdir", "", "Output directory. Will be created if missing.")
	flag.IntVar(&numMapWorkers, "mapworker", defaultNumWorkers, "Default number of map workers.")
	flag.IntVar(&numReduceWorkers, "reduceworker", defaultNumWorkers, "Default number of reduce workers. Must be a power of two.")
	flag.Parse()

	// Check output directory.
	if err := os.MkdirAll(outputDir, os.ModePerm); err != nil {
		fmt.Printf("Error creating output directory: %v\n", err)
		return
	}

	// Set number of reduce workers
	if err := setNumReduceWorker(&numReduceWorkers); err != nil {
		fmt.Println(err)
	}
	fmt.Printf("Setting reduce worker count to %v.\n", numReduceWorkers)

	// Parse input file.
	inFiles, err := parseInFile(*inFile)
	if err != nil {
		fmt.Println(err)
		return
	}

	// Set up channels.
	inputQueue := make(chan string)
	reduceChannels := make([]chan *FiphotLine, numReduceWorkers)
	for i := 0; i < numReduceWorkers; i++ {
		reduceChan := make(chan *FiphotLine)
		reduceChannels[i] = reduceChan
	}

	// Set up work workers.
	var wgMap sync.WaitGroup
	wgMap.Add(numMapWorkers)
	hashFunc := hashTic(numReduceWorkers)
	for i := 0; i < numMapWorkers; i++ {
		go mapWorker(i, inputQueue, reduceChannels, &wgMap, hashFunc)
	}

	// Set up reduce workers.
	var wgReduce sync.WaitGroup
	wgReduce.Add(numReduceWorkers)
	for i, reduceChan := range reduceChannels {
		go reduceWorker(i, reduceChan, &wgReduce)
	}

	// Dispatch jobs to map workers.
	for _, v := range inFiles {
		inputQueue <- v
	}

	// Clean up and wait for all workers to finish.
	close(inputQueue)
	wgMap.Wait()
	for _, v := range reduceChannels {
		close(v)
	}
	wgReduce.Wait()
}
