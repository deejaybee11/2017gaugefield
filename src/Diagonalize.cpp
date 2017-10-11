/**MIT License
*
*Copyright (c) 2016 Dylan
*
*Permission is hereby granted, free of charge, to any person obtaining a copy
*of this software and associated documentation files (the "Software"), to deal
*in the Software without restriction, including without limitation the rights
*to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the Software is
*furnished to do so, subject to the following conditions:
*
*The above copyright notice and this permission notice shall be included in all
*copies or substantial portions of the Software.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*SOFTWARE.
*/
#include "../include/Diagonalize.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <float.h>
#include <iomanip>
#include <fstream>

#include "mkl.h"
#include "mkl_dfti.h"
#include "mkl_lapacke.h"


#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/PotentialData.hpp"
#include "../include/SaveData.hpp"

#define Nm 3

int diagonalize_hamiltonian(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data, int count) {
	
	//The three variables for storing Hamiltonian diagonals
	double H1, H2, H3;
	double min;
	double *delta;
	delta = (double*)mkl_malloc(sim_data.get_num_y() * sizeof(double), 64);
	double min_eigenvalue = 0;
	double *p_minus_a = 0;
	p_minus_a = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);

	if (count >= sim_data.detuning_ramp_time) {
		goto cleanup;
	}
	else {
		for (int i = 0; i < sim_data.get_num_y(); ++i) {
			delta[i] = sim_data.detuning_ramp_shape[count] * sim_data.detuning_gradient * sim_data.y[i];
		}

	}
		

	//diagonalization variables. Refer to documentation for what they are (I forgot)
	MKL_INT info, n = Nm, lda = 3, ldvl = 3, ldvr = 3;
	MKL_Complex16 eigenvalues[Nm], vl[Nm * Nm], vr[Nm * Nm];
	//detuning vector
	

	int index;
	#pragma omp parallel for private(index, min, H1, H2, H3, eigenvalues)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {

			index = i*sim_data.get_num_y() + j;
			H1 = 0.5 * pow(sim_data.px[i] + 2 * sim_data.recoil_k, 2.0) - delta[j];
			H2 = 0.5 * pow(sim_data.px[i], 2.0);
			H3 = 0.5 * pow(sim_data.px[i] - 2 * sim_data.recoil_k, 2.0) + delta[j];

			MKL_Complex16 hamiltonian[9] = {{H1, 0}, {0.5 * sim_data.omega_r, 0}, {0, 0}, {0.5 * sim_data.omega_r, 0}, {H2, 0}, {0.5 * sim_data.omega_r, 0}, {0, 0}, {0.5 * sim_data.omega_r, 0}, {H3, 0}};

			info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'N', n, hamiltonian, lda, eigenvalues, vl, ldvl, vr, ldvr);

			if (eigenvalues[0].real < eigenvalues[1].real) {
				min = std::min(eigenvalues[0].real, eigenvalues[2].real);
			}
			else {
				min = std::min(eigenvalues[1].real, eigenvalues[2].real);
			}
			
//			std::cout << "eval =  [" << eigenvalues[0].real << ", " << eigenvalues[1].real << ", " << eigenvalues[2].real << "]" << std::endl;
			p_minus_a[index] = min;
			
		}
	}
	
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		pot_data.kinetic_energy_y[i] = p_minus_a[i];
	}

	system("rm test.fit");
	system("rm deltatest.txt");
	save_fits_potential(sim_data, p_minus_a, "test.fit");
	save_binary(delta, "deltatest.txt", sim_data.get_num_y());

	goto cleanup;

	cleanup:
		mkl_free(delta);
		mkl_free(p_minus_a);
		if (info != 0) exit(EXIT_FAILURE);
		return 0;
}
