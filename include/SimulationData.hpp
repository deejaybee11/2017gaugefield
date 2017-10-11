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

#ifndef _GAUGE_SIM_DATA_H
#define _GAUGE_SIM_DATA_H

//Standard includes
#include <stdlib.h>
#include <iostream>
//MKL includes
#include "mkl.h"

/**
 * SimulationData class
 *
 * Stores all relevant variables and arrays used for initialization of simulation
 */

class SimulationData {
public:
	//Constructor and Destructor
	SimulationData(int num_x, int num_y);
	~SimulationData();
	//Folder for saving results
	time_t current_time;
	struct tm *mytime;
	char folder[80];
	char command[80];
	int count;
	int max_count = 200;
	//Position and Momentum
	double *x;
	double *y;
	double *px;
	double *py;
	double length_x;
	double length_y;
	//Discretization
	double dx;
	double dx2;
	double dy;
	double dy2;
	double dt;
	//Potential parameters
	double gamma_x;
	double gamma_y;
	//Energy things
	double beta;
	double *chemical_potential;
	//Gauge potential parameters
	double omega_r;
	double detuning_gradient;
	double recoil_k;
	double *detuning_ramp_shape;
	int detuning_ramp_time;
	//Simulation parameters
	int num_r_steps;
	int num_i_steps;
	const char *R = "REAL";
	const char *I = "IMAG";
	
	//Member functions
	int get_num_x() { return this->num_x; };
	int get_num_y() { return this->num_y; };
	int get_total_pts() { return this->total_num_pts; };
	
private:
	int num_x;
	int num_y;
	int total_num_pts;
};

#endif		//	_GAUGE_SIM_DATA_H
