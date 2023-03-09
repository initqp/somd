/*
 * SOMD is an ab-initio molecular dynamics package designed for the SIESTA code.
 * Copyright (C) 2023 github.com/initqp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <vector>
#include <cstdio>
#include <iostream>
#include "nhc.h"
#include "constants.h"

using namespace std;

// T : temperature, unit : K
// t : coupling time scale, unit : ps
// n_bead : chain length
// n_d : DOF
// n_r : number of RESPA iterations.
NHC::NHC(double T, double t, int n_bead, int n_d, int n_r)
{
	int i = 0;

	tau = t;
	n_dof = n_d;
	length = n_bead;
	n_respa = n_r;
	temperature = T;
	Q = vector<double>(length);
	p = vector<double>(length);
	q = vector<double>(length);
	for (i = 0; i < length; i++) {
		p[i] = 1.0 - 2.0 * (double)((i % 2) == 0);
		q[i] = 0.0;
	}
	set_Q();
}

NHC::~NHC(void)
{
	return;
}

// Return the chain length.
int const &NHC::get_length(void) const
{
	return length;
}

// Return the degree of freedom.
int const &NHC::get_n_dof(void) const
{
	return n_dof;
}

// Return number of RESPA.
int const &NHC::get_n_respa(void) const
{
	return n_respa;
}

// Return the coupling time scale.
double const &NHC::get_tau(void) const
{
	return tau;
}

// Return chain temperature.
double const &NHC::get_temperature(void) const
{
	return temperature;
}

// Return the chain masses.
std::vector<double> const &NHC::get_Q(void) const
{
	return Q;
}

// Return the chain positions.
std::vector<double> const &NHC::get_q(void) const
{
	return q;
}

// Return the chain momentums.
std::vector<double> const &NHC::get_p(void) const
{
	return p;
}

// Set the chain positions.
void NHC::set_q(vector<double> &v)
{
	int i = 0;
	int len = 0;

	if ((size_t)length >= v.size())
		len = v.size();
	else
		len = length;
	for (i = 0; i < len; i++)
		q[i] = v[i];
}

// Set the chain momentums.
void NHC::set_p(vector<double> &v)
{
	int i = 0;
	int len = 0;

	if ((size_t)length >= v.size())
		len = v.size();
	else
		len = length;
	for (i = 0; i < len; i++)
		p[i] = v[i];
}

// Set the degree of freedom.
void NHC::set_n_dof(int n_d)
{
	n_dof = n_d;
	set_Q();
};

// Set number of RESPA.
void NHC::set_n_respa(int n_r)
{
	n_respa = n_r;
};

// Set the coupling time scale.
void NHC::set_tau(double t)
{
	tau = t;
	set_Q();
};

// Set the chain temperature.
void NHC::set_temperature(double T)
{
	temperature = T;
	set_Q();
};

void NHC::set_Q(void)
{
	int i = 0;

	for (i = 0; i < length; i++)
		Q[i] = temperature * BOLTZCONST * tau * tau;
	Q[0] *= n_dof;
}

double NHC::propagate(double E_k, double dt)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double G = 0.0;
	double T = 0.0;
	double tmp = 0.0;
	double dt_1 = 0.0;
	double dt_2 = 0.0;
	double dt_4 = 0.0;
	double factor = 1.0;

	T = BOLTZCONST * temperature;
	// Loop over the seven Suzuki-Yoshida numbers.
	for (i = 0; i < 7; i++) {
		dt_1 = dt * w_sy[i] / n_respa;
		dt_2 = 0.5 * dt_1;
		dt_4 = 0.5 * dt_2;
		// Loop over the RESPA iterations.
		for (j = 0; j < n_respa; j++) {
			// update momentum of bead[length - 1]
			if (length == 1)
				G = E_k * 2.0 - n_dof * T;
			else
				G = p[length - 2] * p[length - 2] / \
					Q[length - 2] - T;
			p[length - 1] += G * dt_2;
			// update momentums of bead[length - 2:0]
			for (k = length - 2; k >= 0; k--) {
				tmp = exp(-1.0 * dt_4 * p[k + 1] / Q[k + 1]);
				p[k] *= tmp;
				if (k == 0)
					G = E_k * 2.0 - n_dof * T;
				else
					G = p[k - 1] * p[k - 1] / Q[k - 1] - T;
				p[k] += G * dt_2;
				p[k] *= tmp;
			}
			// update positions of bead[0:length - 1]
			for (k = 0; k < length; k++)
				q[k] += p[k] / Q[k] * dt_1;
			// update the scale factor
			tmp = exp(-1.0 * dt_1 * p[0] / Q[0]);
			factor *= tmp;
			// Since we did not really scale the velocities, we
			// should scale the E_k instead.
			E_k *= tmp * tmp;
			// update momentums of bead[0:length - 2]
			for (k = 0; k <= length - 2; k++) {
				tmp = exp(-1.0 * dt_4 * p[k + 1] / Q[k + 1]);
				p[k] *= tmp;
				if (k == 0)
					G = E_k * 2.0 - n_dof * T;
				else
					G = p[k - 1] * p[k - 1] / Q[k - 1] - T;
				p[k] += G * dt_2;
				p[k] *= tmp;
			}
			// update momentum of bead[length - 1]
			if (length == 1)
				G = E_k * 2.0 - n_dof * T;
			else
				G = p[length - 2] * p[length - 2] / \
					Q[length - 2] - T;
			p[length - 1] += G * dt_2;
		}
	}
	return factor;
}
