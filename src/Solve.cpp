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

#include "../include/Solve.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <cstring>

#include "mkl.h"

#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/PotentialData.hpp"
#include "../include/SaveData.hpp"
#include "../include/Diagonalize.hpp"

void calculate_ground_state(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data) {

	//Checks if the groundstate file already exists
	struct stat sb;
	if (stat("GroundState.fit", &sb) == 0) {
		system("rm GroundState.fit");
		std::cout << "GroundState.fit deleted" << std::endl;
	}

	//Creates all neccessary handles for the FFT routine
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG N[2]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y();
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0 / (N[0]*N[1])));
	status = DftiCommitDescriptor(handle);

	// Build dressed dispersion at δ=0 so the imaginary-time ground state has
	// the correct anisotropic effective mass from the Raman coupling.
	diagonalize_hamiltonian(sim_data, psi, pot_data, 0);
	pot_data.assign_momentum_operator(sim_data, psi, false);

	// Same stride convention as real-time handles above.
	DFTI_DESCRIPTOR_HANDLE handle_x_im = 0, handle_y_im = 0;
	MKL_LONG x_strides_im[] = {0, N[1]};
	DftiCreateDescriptor(&handle_x_im, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[0]);
	DftiSetValue(handle_x_im, DFTI_NUMBER_OF_TRANSFORMS, N[1]);
	DftiSetValue(handle_x_im, DFTI_BACKWARD_SCALE, 1.0 / N[0]);
	DftiSetValue(handle_x_im, DFTI_INPUT_DISTANCE, 1);
	DftiSetValue(handle_x_im, DFTI_INPUT_STRIDES, x_strides_im);
	DftiCommitDescriptor(handle_x_im);

	DftiCreateDescriptor(&handle_y_im, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[1]);
	DftiSetValue(handle_y_im, DFTI_NUMBER_OF_TRANSFORMS, N[0]);
	DftiSetValue(handle_y_im, DFTI_BACKWARD_SCALE, 1.0 / N[1]);
	DftiSetValue(handle_y_im, DFTI_INPUT_DISTANCE, N[1]);
	DftiCommitDescriptor(handle_y_im);

	//Now we can find the ground state!
	//
	double mu = 0, mu_prev = 0;
	const double mu_tol = 1e-6;
	for (int i = 0; i < sim_data.num_i_steps; ++i) {

		psi.calc_abs(sim_data);
		psi.calc_norm(sim_data);
		mu = -log(psi.last_norm_sq) / (2.0 * sim_data.dt);

		if ((i%1000 == 0)&& (i>0)) {
			std::cout << "Istep " << i << "  mu = " << mu << std::endl;
		}
		if (i > 0 && fabs(mu - mu_prev) < mu_tol * fabs(mu)) {
			std::cout << "Converged at step " << i << "  mu = " << mu << std::endl;
			break;
		}
		mu_prev = mu;
		//Nonlinear calculation
		pot_data.calculate_non_linear_energy(sim_data, psi);
		pot_data.assign_position_operator(sim_data, psi, true, false);
		//Multiply psi by the position operator. dt = dt/2
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
		// x-step: dressed dispersion E_min(px, y) in imaginary time
		DftiComputeForward(handle_x_im, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_x, psi.psi);
		DftiComputeBackward(handle_x_im, psi.psi, psi.psi);
		// y-step: bare p_y²/2 in imaginary time
		DftiComputeForward(handle_y_im, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_y, psi.psi);
		DftiComputeBackward(handle_y_im, psi.psi, psi.psi);
		//second half-step position operator
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
	}

	DftiFreeDescriptor(&handle_x_im);
	DftiFreeDescriptor(&handle_y_im);
	DftiFreeDescriptor(&handle);
	save_wavefunction_binary(psi, sim_data, "GroundState.bin");
	save_fits_wavefunction(sim_data, psi, "GroundState.fit");
}


void calculate_time_evolution(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data) {

	//Creates all neccessary handles for the FFT routine
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle_x = 0;
	DFTI_DESCRIPTOR_HANDLE handle_y = 0;
	MKL_LONG N[2]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y();
	// handle_x: Ny x-direction FFTs. Memory layout psi[ix*Ny+iy] means x-elements
	// are Ny apart; Ny transforms (one per y-slice) are distance 1 apart.
	MKL_LONG x_strides[] = {0, N[1]};
	DftiCreateDescriptor(&handle_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[0]);
	DftiCreateDescriptor(&handle_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[1]);

	DftiSetValue(handle_x, DFTI_NUMBER_OF_TRANSFORMS, N[1]);
	DftiSetValue(handle_x, DFTI_BACKWARD_SCALE, 1.0 / N[0]);
	DftiSetValue(handle_x, DFTI_INPUT_DISTANCE, 1);
	DftiSetValue(handle_x, DFTI_INPUT_STRIDES, x_strides);

	// handle_y: Nx y-direction FFTs. y-elements are contiguous (stride 1);
	// Nx transforms (one per x-slice) are Ny apart.
	DftiSetValue(handle_y, DFTI_NUMBER_OF_TRANSFORMS, N[0]);
	DftiSetValue(handle_y, DFTI_BACKWARD_SCALE, 1.0 / N[1]);
	DftiSetValue(handle_y, DFTI_INPUT_DISTANCE, N[1]);

	DftiCommitDescriptor(handle_x);
	DftiCommitDescriptor(handle_y);

	//Now we can find the ground state!
	//

	// Seed symmetry breaking — needed for vortex nucleation in a deterministic GPE
	srand(42);
	const double noise_amp = 1e-2;
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		psi.psi[i].real += noise_amp * ((double)rand() / RAND_MAX - 0.5);
		psi.psi[i].imag += noise_amp * ((double)rand() / RAND_MAX - 0.5);
	}

	for (int i = 0; i < sim_data.num_r_steps; ++i) {

		diagonalize_hamiltonian(sim_data, psi, pot_data, i);
		pot_data.assign_momentum_operator(sim_data, psi, true);  // update mom_operator_x from kinetic_energy_x each step

		if (i % sim_data.save_interval == 0) {
			int frame = i / sim_data.save_interval;
			std::cout << "Rstep " << i << " out of " << sim_data.num_r_steps << std::endl;
			char psi_path[256], kinx_path[256];
			sprintf(psi_path,  "%s/psi%d.fit",  sim_data.folder, frame);
			sprintf(kinx_path, "%s/kinx%d.fit", sim_data.folder, frame);
			save_fits_wavefunction(sim_data, psi, psi_path);
			save_fits_potential(sim_data, pot_data.kinetic_energy_x, kinx_path);
		}

		//First half-step position operator (trap + nonlinear)
		pot_data.calculate_non_linear_energy(sim_data, psi);
		pot_data.assign_position_operator(sim_data, psi, true, true);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
		//Full-step momentum operator: FFT x, then y
		DftiComputeForward(handle_x, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_x, psi.psi);
		DftiComputeBackward(handle_x, psi.psi, psi.psi);
		DftiComputeForward(handle_y, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_y, psi.psi);
		DftiComputeBackward(handle_y, psi.psi, psi.psi);
		//Second half-step position operator (recompute density after kinetic step)
		pot_data.calculate_non_linear_energy(sim_data, psi);
		pot_data.assign_position_operator(sim_data, psi, true, true);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
	}

	DftiFreeDescriptor(&handle_x);
	DftiFreeDescriptor(&handle_y);

}





