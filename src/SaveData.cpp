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

#include <stdlib.h>
#include <fstream>

#include "mkl.h"

#include "../include/SaveData.hpp"

void save_fits_wavefunction(SimulationData &sim_data, WaveFunction &psi, const char *fits_file_name) {

	double *save_data = 0;
	save_data = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);

	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};
	//Calc abs psi
	psi.calc_abs(sim_data);
	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		save_data[i] = psi.abs_psi[i];
	}
	fits_create_file(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	nelements = naxes[0]*naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);

	mkl_free(save_data);
}

void save_fits_potential(SimulationData &sim_data, double *potential, const char * fits_file_name) {

	fitsfile *fptr;
	int status = 0;
	double *save_data = 0;
	save_data = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};
	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		save_data[i] = potential[i];
	}
	fits_create_file(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	nelements = naxes[0]*naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);

	mkl_free(save_data);
}


void save_binary(double *data, const char *filename, int length) {
	
	remove(filename);
	std::ofstream output(filename);
	output << length << "\n";
	for (int i = 0; i < length; ++i) {
		output << data[i] << "\n";
	}
	output.close();
}


