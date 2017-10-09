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

#ifndef _SAVE_DATA_H
#define _SAVE_DATA_H

#include <stdlib.h>

#include "mkl.h"
#include "fitsio.h"

#include "SimulationData.hpp"
#include "WaveFunction.hpp"
#include "PotentialData.hpp"

void save_fits_wavefunction(SimulationData &sim_data, WaveFunction &psi, const char * fits_file_name);
void save_fits_potential(SimulationData &sim_data, double *pot_data, const char * fits_file_name);
void save_binary(double *data, const char *filename, int length);
void save_binary(MKL_Complex16 *data, const char *filename, int length);

#endif    //    _SAVE_DATA_H
