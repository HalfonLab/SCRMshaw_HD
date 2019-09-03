/* Copyright (c) 2000 by Martin Tompa and Saurabh Sinha.
 * All rights reserved.  Redistribution is not permitted without the 
 * express written permission of the authors.
 * The program YMF implements an algorithm for identifying
 * likely transcription factor binding sites in yeast, described in the
 * following paper:
 "A Statistical Method for Finding Transcription Factor Binding Sites"
 by Saurabh Sinha and Martin Tompa,
 Eighth International Conference on Intelligent Systems for
 Molecular Biology, San Diego, USA, August 2000, 344-354.
 */

#include <stdio.h>
#include "preproc-genome.h"

int main(int argc, char **argv)
{
	Genome genome;

	genome.ReadAllChromosomesAndComputeP(argc, argv);
	genome.ComputeStationaryDistribution();

	return 0;

}
