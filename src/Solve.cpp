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

	pot_data.assign_momentum_operator(sim_data, psi, false);
	//Now we can find the ground state!
	//
	for (int i = 0; i < sim_data.num_i_steps; ++i) {

		if (i%1000 == 0) {
			std::cout << "Istep " << i << " out of " << sim_data.num_i_steps << std::endl;
		}
		
		psi.calc_abs(sim_data);
		psi.calc_norm(sim_data);
		//Nonlinear calculation
		pot_data.calculate_non_linear_energy(sim_data, psi);
		pot_data.assign_position_operator(sim_data, psi, true, false);
		//Multiply psi by the position operator. dt = dt/2
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
		//Perform fft and multiply by the momentum operator
		status = DftiComputeForward(handle, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator, psi.psi);
		//transform back and multiply by pos operator again
		status = DftiComputeBackward(handle, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
	}

	DftiFreeDescriptor(&handle);
	save_fits_wavefunction(sim_data, psi, "GroundState.fit");
}


void calculate_time_evolution(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data) {

	//Creates all neccessary handles for the FFT routine
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle_x = 0;
	DFTI_DESCRIPTOR_HANDLE handle_y = 0;
	MKL_LONG N[2]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y();
	MKL_LONG strides[] = {0, N[0]};
	DftiCreateDescriptor(&handle_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[0]);
	DftiCreateDescriptor(&handle_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, N[1]);
	DftiSetValue(handle_x, DFTI_NUMBER_OF_TRANSFORMS, N[1]); 
	DftiSetValue(handle_x, DFTI_BACKWARD_SCALE, (1.0 / (N[0])));
	DftiSetValue(handle_x, DFTI_INPUT_DISTANCE, N[0]);
	DftiSetValue(handle_y, DFTI_NUMBER_OF_TRANSFORMS, N[0]); 
	DftiSetValue(handle_y, DFTI_BACKWARD_SCALE, (1.0 / (N[1])));
	DftiSetValue(handle_y, DFTI_INPUT_DISTANCE, 1);
	DftiSetValue(handle_y, DFTI_INPUT_STRIDES, strides);
	DftiCommitDescriptor(handle_x);
	DftiCommitDescriptor(handle_y);

	pot_data.assign_momentum_operator(sim_data, psi, true);
	//Now we can find the ground state!
	//

	for (int i = 0; i < sim_data.num_r_steps; ++i) {

		diagonalize_hamiltonian(sim_data, psi, pot_data, i);
		if (i%1000 == 0) {
			std::cout << "Rstep " << i << " out of " << sim_data.num_r_steps << std::endl;
			char buf1[200];
			char buf2[200];
			strcpy(buf1, sim_data.folder);
			sprintf(buf2, "/psi%d.fit\0", i/1000);
			strcat(buf1, buf2);
			int length = strlen(buf1);
			char *full_filename;
			full_filename = (char*)mkl_malloc(length * sizeof(char), 64);
			strcpy(full_filename, buf1);
			save_fits_wavefunction(sim_data, psi, full_filename);
			mkl_free(full_filename);

		}
		
		psi.calc_abs(sim_data);
		psi.calc_norm(sim_data);
		//Nonlinear calculation
		pot_data.calculate_non_linear_energy(sim_data, psi);
		pot_data.assign_position_operator(sim_data, psi, false, true);
		//Multiply psi by the position operator. dt = dt/2
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
		//Perform fft in x and multiply by the momentum operator
		DftiComputeForward(handle_x, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_x, psi.psi);
		DftiComputeBackward(handle_x, psi.psi, psi.psi);
		//then do y transform
		DftiComputeForward(handle_y, psi.psi, psi.psi);
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.mom_operator_y, psi.psi);
		DftiComputeBackward(handle_y, psi.psi, psi.psi);
		//transform back and multiply by pos operator again
		vzMul(sim_data.get_total_pts(), psi.psi, pot_data.pos_operator, psi.psi);
	}

	DftiFreeDescriptor(&handle_x);
	DftiFreeDescriptor(&handle_y);

}





