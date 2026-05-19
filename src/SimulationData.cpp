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
#include "../include/SimulationData.hpp"
//Standard includes
#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <string.h>
#include <iomanip>

#include "../include/SaveData.hpp"

// ── User-configurable parameters ─────────────────────────────────────────────
static const double TRAP_FREQ_HZ    = 70.0;   // trap frequency [Hz]
static const double MASS_AMU        = 87.0;    // atomic mass [amu]
static const double LASER_NM        = 790.0;   // Raman laser wavelength [nm]
static const double OMEGA_R_EREC    = 8.2;     // Rabi coupling [E_rec]
static const double GRAD_KHZ_PER_UM = -0.6;     // detuning gradient [kHz/μm]
static const double RAMP_TOTAL_MS   = 80;    // total ramp duration [ms]
static const double RAMP_WIDTH_MS   = 60;    // tanh ramp width [ms]
static const double RAMP_SHIFT_MS   = 5;    // hold before ramp starts [ms]
static const double BOX_UM          = 35;    // simulation box size [μm]
static const double DT_US           = 1.0;     // time step [μs]
static const double RUN_TIME_MS     = 800.0;   // real-time evolution duration [ms]
static const double SAVE_INTERVAL_MS = 0.5;    // wavefunction save interval [ms]
static const double N_ATOMS         = 1e5;     // number of atoms
static const double A_SCATT_BOHR    = 100.4;   // s-wave scattering length [a_0] (Rb-87 ≈ 100.4)
static const double OMEGA_Z_HZ      = 100.0;   // axial trap frequency [Hz] (sets ℓ_z)
static const double GAMMA_X         = 1.0;     // trap anisotropy x
static const double GAMMA_Y         = 1.0;     // trap anisotropy y
// ─────────────────────────────────────────────────────────────────────────────

//Constructor
SimulationData::SimulationData(int num_x, int num_y) {
	//Folder based on current time
	time(&this->current_time);
	this->mytime = localtime(&this->current_time);
	sprintf(this->command, "mkdir -p fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + mytime->tm_mon, mytime->tm_mday, mytime->tm_hour, mytime->tm_min);
	system(this->command);
	printf("Directory created\n");
	sprintf(this->folder, "fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + mytime->tm_mon, mytime->tm_mday, mytime->tm_hour, mytime->tm_min);
	//Physical constants and trap parameters (Rb-87)
	const double hbar       = 1.0546e-34;
	const double mass       = MASS_AMU * 1.66054e-27;
	const double omega_trap = 2.0 * M_PI * TRAP_FREQ_HZ;
	double lx = sqrt(hbar / (mass * omega_trap));   // oscillator length [m]

	//Simulation grid
	this->num_x = num_x;
	this->num_y = num_y;
	this->total_num_pts = num_x * num_y;
	this->length_x = BOX_UM * 1e-6 / lx;
	this->length_y = BOX_UM * 1e-6 / lx;
	//Number of steps
	this->num_r_steps    = (int)(RUN_TIME_MS      * 1000.0 / DT_US);
	this->save_interval  = (int)(SAVE_INTERVAL_MS * 1000.0 / DT_US);
	this->num_i_steps = 1000000;
	this->gamma_x = GAMMA_X;
	this->gamma_y = GAMMA_Y;
	const double a_bohr = 5.2918e-11;                         // Bohr radius [m]
	const double a_s    = A_SCATT_BOHR * a_bohr;              // scattering length [m]
	const double omega_z = 2.0 * M_PI * OMEGA_Z_HZ;
	const double ell_z  = sqrt(hbar / (mass * omega_z));      // axial osc. length [m]
	this->beta = sqrt(8.0 * M_PI) * N_ATOMS * a_s / ell_z;   // β = √(8π)·N·a/ℓ_z
	printf("Interaction: N=%.0f, a=%.1f a0, ell_z=%.3f um, beta=%.1f\n",
	       N_ATOMS, A_SCATT_BOHR, ell_z*1e6, this->beta);
	this->count = 0;
	//Gauge potential parameters (specify in physical units, converted to dimensionless trap units)
	double k_L              = 2.0 * M_PI / (LASER_NM * 1e-9);  // laser wavevector [m⁻¹]
	this->recoil_k          = k_L * lx;                         // k_L in units of 1/ℓ
	double E_rec            = 0.5 * pow(this->recoil_k, 2.0);  // E_rec in units of ℏω_trap
	this->omega_r           = OMEGA_R_EREC * E_rec;             // E_rec → ℏω_trap
	this->detuning_gradient = GRAD_KHZ_PER_UM * (2.0 * M_PI * 1e9 * lx / omega_trap); // kHz/μm → ℏω_trap/ℓ
	printf("Unit conversion: lx=%.3e m, k_rec=%.3f/ell, E_rec=%.2f hOmega, Omega_r=%.1f hOmega\n",
	       lx, this->recoil_k, E_rec, this->omega_r);
	//Memory allocation time
	this->x = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->y = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
	this->px = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->py = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->chemical_potential = (double*)mkl_malloc(this->num_i_steps * sizeof(double), 64);
	//Steps for filling the momentum arrays with correct values
	double ax = -0.5 * this->num_x;
	double bx = 0.5 * this->num_x - 1.0;
	double ay = -0.5 * this->num_y;
	double by = 0.5 * this->num_y - 1.0;
	double step_x = (2 * M_PI / this->length_x) * ((bx - ax) / (this->num_x - 1.0));
	double step_y = (2 * M_PI / this->length_y) * ((by - ay) / (this->num_y - 1.0));
	//Fill the newly allocated arrays
	for (int i = 0; i < this->num_x; ++i) {
		this->x[i] = -0.5 * this->length_x + i * this->length_x / ((double)this->num_x);
		this->px[i] = (2 * M_PI / this->length_x) * ax + i * step_x;
	}	
	for (int i = 0; i < this->num_y; ++i) {
		this->y[i] = -0.5 * this->length_y + i * this->length_y / ((double)this->num_y);
		this->py[i] = (2 * M_PI / this->length_y) * ay + i * step_y;
	}
	this->dx = this->x[2]-this->x[1];
	this->dy = this->y[2]-this->y[1];
	this->dt      = DT_US * 1e-6 * omega_trap;
	this->ell_um  = lx * 1e6;
	this->dt_us   = DT_US;
	this->trap_hz = TRAP_FREQ_HZ;
	//Perform FFT shift on momentum arrays to compensate for negative frequencies shifting
	double temp;
	int n[2];
	n[0] = this->get_num_x() / 2.0;
	n[1] = this->get_num_y() / 2.0;
	for (int i = 0; i < n[0]; ++i) {
		temp = this->px[i];
		this->px[i] = this->px[i + n[0]];
		this->px[i + n[0]] = temp;
	}
	for (int i = 0; i < n[1]; ++i) {
		temp = this->py[i];
		this->py[i] = this->py[i + n[1]];
		this->py[i + n[1]] = temp;
	}

	double ms_to_steps       = 1e-3 * omega_trap / this->dt;
	this->detuning_ramp_time = (int)(RAMP_TOTAL_MS * ms_to_steps);
	this->detuning_ramp_shape = (double*)mkl_malloc(this->detuning_ramp_time * sizeof(double), 64);
	double temp_val = 0;
	double ramp_width = RAMP_WIDTH_MS * ms_to_steps;
	double shift      = RAMP_SHIFT_MS * ms_to_steps;

	for (int i = 0; i < this->detuning_ramp_time; ++i) {

		temp_val = i;
		if (i < shift) {
			this->detuning_ramp_shape[i] = 0;
		}
		else {
			this->detuning_ramp_shape[i] = 0.5 * (1.0 + tanh(6.0 * (temp_val - shift - 0.5*ramp_width) / ramp_width));
		}
	}

	save_binary(this->detuning_ramp_shape, "detuning_ramp.txt", this->detuning_ramp_time);	

};

//Destructor for the class
SimulationData::~SimulationData() {
	mkl_free(x);
	mkl_free(y);
	mkl_free(px);
	mkl_free(py);
};
		


















