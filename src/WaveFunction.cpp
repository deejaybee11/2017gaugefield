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

#include "../include/WaveFunction.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <sys/stat.h>

#include "../include/SimulationData.hpp"
#include "../include/PotentialData.hpp"
#include "../include/SaveData.hpp"

WaveFunction::WaveFunction(SimulationData &sim_data, double *harmonic_trap) {

	struct stat sb;
	if (stat("InitialState.fit", &sb) == 0) {
		system("rm InitialState.fit");
		std::cout << "InitialState.fit deleted" << std::endl;
	}

	this->psi = (MKL_Complex16*)mkl_malloc(sim_data.get_total_pts() * sizeof(MKL_Complex16), 64);
	this->psi_conj = (MKL_Complex16*)mkl_malloc(sim_data.get_total_pts() * sizeof(MKL_Complex16), 64);
	this->abs_psi = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);

	int index = 0;

	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {

			index = i*sim_data.get_num_y() + j;
			if (sqrt(pow(sim_data.x[i], 2.0) + pow(sim_data.y[j], 2.0) < 5)) {
				this->psi[index].real = 1;
				this->psi[index].imag = 0;
			}
			else {
				this->psi[index].real = 0;
				this->psi[index].imag = 0;
			}
		}
	}
	calc_abs(sim_data);
	save_fits_potential(sim_data, this->abs_psi, "InitialState.fit");

}	

void WaveFunction::calc_abs(SimulationData &sim_data) {

	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		this->abs_psi[i] = this->psi[i].real * this->psi[i].real + this->psi[i].imag * this->psi[i].imag;
	}
}

void WaveFunction::calc_norm(SimulationData &sim_data) {
	
	double psi_sum = 0;
	double temp_real = 0;
	double temp_imag = 0;
	double norm_fac = 0;
	
	#pragma omp parallel for reduction(+:psi_sum)
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		psi_sum += this->abs_psi[i];
	}

	norm_fac = sqrt(1.0 / (psi_sum * sim_data.dx * sim_data.dy));
	std::cout << norm_fac << std::endl;
	#pragma omp parallel for private(temp_real, temp_imag)
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		temp_real = this->psi[i].real * norm_fac;
		temp_imag = this->psi[i].imag * norm_fac;
		this->psi[i].real = temp_real;
		this->psi[i].real = temp_imag;
	}
}

WaveFunction::~WaveFunction() {

	mkl_free(psi);
	mkl_free(abs_psi);
	mkl_free(psi_conj);
}

