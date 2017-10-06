/**The MIT License (MIT)
 *
 * Copyright (c) 2016 Dylan
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

//Standard includes
#include <stdlib.h>
#include <iostream>
//MKL includes
#include "mkl.h"
//Header file includes
#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/PotentialData.hpp"
#include "../include/Solve.hpp"
#include "../include/Diagonalize.hpp"

#if !defined(MKL_ILP64)
	#define LI "%li"
	#else
	#define LE "%lli"
#endif

int main() {

	putenv("MKP_BLOCKTIME=0");
	putenv("KMP_AFFINITY=verbose,granularity=fine,compact,norespect");
	mkl_set_num_threads(mkl_get_max_threads());
	mkl_disable_fast_mm();

	bool load_bin = false;

	SimulationData sim_data(128, 128);
	PotentialData pot_data(sim_data);
	WaveFunction psi(sim_data, pot_data.harmonic_trap);

	//Find the ground state
	calculate_ground_state(sim_data, psi, pot_data);
	diagonalize_hamiltonian(sim_data, psi, pot_data);
	calculate_time_evolution(sim_data, psi, pot_data);
	return 0;

}

