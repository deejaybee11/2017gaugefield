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

static void fits_open_overwrite(fitsfile **fptr, const char *filename, int *status) {
	char buf[512];
	snprintf(buf, sizeof(buf), "!%s", filename);
	fits_create_file(fptr, buf, status);
}

static void write_sim_header(fitsfile *fptr, SimulationData &sim_data, int *status) {
	fits_write_key(fptr, TDOUBLE, "ELL_UM",   &sim_data.ell_um,            "oscillator length [um]",       status);
	fits_write_key(fptr, TDOUBLE, "TRAP_HZ",  &sim_data.trap_hz,           "trap frequency [Hz]",          status);
	fits_write_key(fptr, TDOUBLE, "K_REC",    &sim_data.recoil_k,          "recoil momentum [hbar/ell]",   status);
	fits_write_key(fptr, TDOUBLE, "OMEGA_R",  &sim_data.omega_r,           "Rabi coupling [hbar*omega]",   status);
	fits_write_key(fptr, TDOUBLE, "DET_GRAD", &sim_data.detuning_gradient, "detuning gradient [hbar*w/l]", status);
	fits_write_key(fptr, TDOUBLE, "HALFBOX",  &sim_data.length_x,          "full box length [ell]",        status);
	fits_write_key(fptr, TDOUBLE, "DT_US",    &sim_data.dt_us,             "timestep [us]",                status);
	int savintv = sim_data.save_interval;
	fits_write_key(fptr, TINT,    "SAVINTV",  &savintv,                    "save interval [steps]",        status);
}

void save_fits_wavefunction(SimulationData &sim_data, WaveFunction &psi, const char *fits_file_name) {

	double *save_data = 0;
	save_data = (double*)mkl_malloc(sim_data.get_total_pts() * sizeof(double), 64);

	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};
	psi.calc_abs(sim_data);
	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_total_pts(); ++i) {
		save_data[i] = psi.abs_psi[i];
	}
	fits_open_overwrite(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	write_sim_header(fptr, sim_data, &status);
	nelements = naxes[0]*naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	if (status) fits_report_error(stderr, status);

	mkl_free(save_data);
}

void save_fits_potential(SimulationData &sim_data, double *potential, const char *fits_file_name) {

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
	fits_open_overwrite(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	write_sim_header(fptr, sim_data, &status);
	nelements = naxes[0]*naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	if (status) fits_report_error(stderr, status);

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

void save_wavefunction_binary(WaveFunction &psi, SimulationData &sim_data, const char *filename) {
	std::ofstream f(filename, std::ios::binary);
	int nx = sim_data.get_num_x(), ny = sim_data.get_num_y();
	f.write(reinterpret_cast<const char*>(&nx), sizeof(int));
	f.write(reinterpret_cast<const char*>(&ny), sizeof(int));
	f.write(reinterpret_cast<const char*>(psi.psi), nx * ny * sizeof(MKL_Complex16));
	f.close();
	std::cout << "Ground state saved to " << filename << std::endl;
}

bool load_wavefunction_binary(WaveFunction &psi, SimulationData &sim_data, const char *filename) {
	std::ifstream f(filename, std::ios::binary);
	if (!f) return false;
	int nx, ny;
	f.read(reinterpret_cast<char*>(&nx), sizeof(int));
	f.read(reinterpret_cast<char*>(&ny), sizeof(int));
	if (nx != sim_data.get_num_x() || ny != sim_data.get_num_y()) {
		std::cerr << "Grid mismatch in " << filename << ": saved "
		          << nx << "x" << ny << ", expected "
		          << sim_data.get_num_x() << "x" << sim_data.get_num_y() << std::endl;
		return false;
	}
	f.read(reinterpret_cast<char*>(psi.psi), nx * ny * sizeof(MKL_Complex16));
	if (!f) {
		std::cerr << "Failed to read wavefunction data from " << filename << std::endl;
		return false;
	}
	psi.calc_abs(sim_data);
	std::cout << "Ground state loaded from " << filename << std::endl;
	return true;
}


