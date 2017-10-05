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

//Constructor
SimulationData::SimulationData(int num_x, int num_y) {
	//Folder based on current time
	time(&this->current_time);
	this->mytime = localtime(&this->current_time);
	sprintf(this->command, "mkdir fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + mytime->tm_mon, mytime->tm_mday, mytime->tm_hour, mytime->tm_min);
	system(this->command);
	printf("Directory created\n");
	sprintf(this->folder, "fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + mytime->tm_mon, mytime->tm_mday, mytime->tm_hour, mytime->tm_min);
	//Simulation grid
	this->num_x = num_x;
	this->num_y = num_y;
	this->total_num_pts = num_x * num_y;
	this->length_x = 30;
	this->length_y = 30;
	//Number of steps
	this->num_r_steps = 100000;
	this->num_i_steps = 5;
	this->gamma_x = 1;
	this->gamma_y = 1.4;
	this->beta = 1;
	double dt = 0.01;
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
};

//Destructor for the class
SimulationData::~SimulationData() {
	mkl_free(x);
	mkl_free(y);
	mkl_free(px);
	mkl_free(py);
};
		


















