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

//Header include
#include "../include/PotentialData.hpp"
//Standard includes
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
//MKL includes
#include "mkl.h"
//Other includes
#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/SaveData.hpp"

PotentialData::PotentialData(SimulationData &sim_data) {
	//Allocate memory for the potential arrays
	this->harmonic_trap = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);
	this->non_linear = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);
	this->kinetic_energy = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);
	this->pos_operator = (MKL_Complex16*)mkl_malloc(sim_data.get_total_pts() * sizeof(MKL_Complex16), 64);
	this->mom_operator = (MKL_Complex16*)mkl_malloc(sim_data.get_total_pts() * sizeof(MKL_Complex16), 64);

	//Fill arrays
	double harmonic_val = 0;
	double kin_energy_val = 0;
	int index;
	#pragma omp parallel for private(index, harmonic_val)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			index = i * sim_data.get_num_y() + j;
			harmonic_val = 0.5 * (pow(sim_data.gamma_x, 2.0) * pow(sim_data.x[i], 2.0) + pow(sim_data.gamma_y, 2.0) * pow(sim_data.y[j], 2.0));
			kin_energy_val = 0.5 * (pow(sim_data.px[i], 2.0) + pow(sim_data.py[j], 2.0));
			
			this->harmonic_trap[index] = harmonic_val;
			this->kinetic_energy[index] = kin_energy_val;
		}
	} 

}

void PotentialData::assign_position_operator(SimulationData &sim_data, WaveFunction &psi, bool trap_on, bool real_time) {

	double theta = 0;
	double harm = 0;

	if (trap_on) harm = 1;

	if (real_time) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_total_pts(); ++i) {
			theta = (this->non_linear[i] + harm * this->harmonic_trap[i]) * 0.5 * sim_data.dt;
			this->pos_operator[i].real = cos(theta);
			this->pos_operator[i].imag = -1.0 * sin(theta);
		}
	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_total_pts(); ++i) {
			theta = (this->non_linear[i] + harm * this->harmonic_trap[i]) * 0.5 * sim_data.dt;
			//std::cout << std::setprecision(16) << "theta=" << theta << " dt=" << sim_data.dt <<  std::endl;
			this->pos_operator[i].real = exp(-1.0 * theta);
			this->pos_operator[i].imag = 0;
		}
//	std::cout << pos_operator[1562].real << " " << pos_operator[1562].imag << std::endl;
	
	}
}

void PotentialData::assign_momentum_operator(SimulationData &sim_data, WaveFunction &psi, bool real_time) {

	double theta = 0;

	if (real_time) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_total_pts(); ++i) {
			theta = (kinetic_energy[i]) * sim_data.dt;
			this->mom_operator[i].real = cos(theta);
			this->mom_operator[i].imag = -1.0 * sin(theta);
		}
	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_total_pts(); ++i) {
			theta = (kinetic_energy[i]) * sim_data.dt;
			this->mom_operator[i].real = exp(-1.0 * theta);
			this->mom_operator[i].imag = 0;
		}
//	std::cout << mom_operator[1562].real << " " << mom_operator[1562].imag << std::endl;

	}
}

void PotentialData::calculate_non_linear_energy(SimulationData &sim_data, WaveFunction &psi) {

	double nonlinearval = 0;
	psi.calc_abs(sim_data);
	#pragma omp parallel for private(nonlinearval)
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		nonlinearval = sim_data.beta * psi.abs_psi[i];
		this->non_linear[i] = nonlinearval;
	}
}

PotentialData::~PotentialData() {
	mkl_free(harmonic_trap);
	mkl_free(non_linear);
	mkl_free(mom_operator);
	mkl_free(pos_operator);
}
