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

#include <omp.h>
#include <cmath>
#include <iostream>
#include "math_utils.h"
#include "constraints.h"

using namespace std;

// Return the distance between atoms.
static inline double _distance(double *positions, vector<int> indices)
{
	return vecdis(&(positions[indices[0] * 3]),
		&(positions[indices[1] * 3]), 3);
}

// Return the angle between three atoms.
static inline double _angle(double *positions, vector<int> indices)
{
	int idx_1 = 0;
	int idx_2 = 0;
	int idx_3 = 0;
	double a[3];
	double b[3];

	idx_1 = indices[0] * 3;
	idx_2 = indices[1] * 3;
	idx_3 = indices[2] * 3;
	vecsub(&(positions[idx_2]), &(positions[idx_1]), a, 3);
	vecsub(&(positions[idx_2]), &(positions[idx_3]), b, 3);
	return acos(vecdot(a, b, 3) / ( \
		vecdis(&(positions[idx_2]), &(positions[idx_1]), 3) * \
		vecdis(&(positions[idx_2]), &(positions[idx_3]), 3)));
}

// Return the torsion between four atoms
static inline double _torsion(double *positions, vector<int> indices)
{
	int idx_1 = 0;
	int idx_2 = 0;
	int idx_3 = 0;
	int idx_4 = 0;
	double n_1[3];
	double n_2[3];
	double m_1[3];
	double r_12[3];
	double r_23[3];
	double r_34[3];

	idx_1 = indices[0] * 3;
	idx_2 = indices[1] * 3;
	idx_3 = indices[2] * 3;
	idx_4 = indices[3] * 3;
	vecsub(&(positions[idx_2]), &(positions[idx_1]), r_12, 3);
	vecsub(&(positions[idx_3]), &(positions[idx_2]), r_23, 3);
	vecsub(&(positions[idx_4]), &(positions[idx_3]), r_34, 3);
	veccross3(r_12, r_23, n_1);
	veccross3(r_23, r_34, n_2);
	vecnorm(n_1,  3);
	vecnorm(n_2,  3);
	vecnorm(r_23, 3);
	veccross3(n_1, r_23, m_1);
	return (-1.0 * atan2(vecdot(m_1, n_2, 3), vecdot(n_1, n_2, 3)));
}

// Calculate derivatives: d_distance / d_r.
static inline void _d_distance_d_r(double *positions, vector<int> indices, \
	double value, double *deriv)
{
	int i = 0;
	int idx_1 = 0;
	int idx_2 = 0;

	idx_1 = indices[0] * 3;
	idx_2 = indices[1] * 3;
	for (i = 0; i < 3; i++) {
		deriv[i + 0] = \
			(positions[idx_1 + i] - positions[idx_2 + i]) / value;
		deriv[i + 3] = \
			(positions[idx_2 + i] - positions[idx_1 + i]) / value;
	}
}

// Calculate derivatives: d_angle / d_r.
static inline void _d_angle_d_r(double *positions, vector<int> indices, \
	double value, double *deriv)
{
	int i = 0;
	int idx_1 = 0;
	int idx_2 = 0;
	int idx_3 = 0;
	double r_12 = 0.0;
	double r_32 = 0.0;
	double cos_angle = 0.0;
	double sin_angle = 0.0;

	idx_1 = indices[0] * 3;
	idx_2 = indices[1] * 3;
	idx_3 = indices[2] * 3;
	cos_angle = cos(value);
	sin_angle = sin(value);
	r_12 = vecdis(&(positions[idx_1]), &(positions[idx_2]), 3);
	r_32 = vecdis(&(positions[idx_3]), &(positions[idx_2]), 3);
	for (i = 0; i < 3; i++) {
		deriv[i + 0] = \
			((positions[idx_1 + i] - positions[idx_2 + i]) / \
			r_12 * cos_angle - (positions[idx_3 + i] - \
			positions[idx_2 + i]) / r_32) / r_12 / sin_angle;
		deriv[i + 6] = \
			((positions[idx_3 + i] - positions[idx_2 + i]) / \
			r_32 * cos_angle - (positions[idx_1 + i] - \
			positions[idx_2 + i]) / r_12) / r_32 / sin_angle;
		deriv[i + 3] = -1.0 * deriv[i + 0] - deriv[i + 6];
	}
}

// Calculate derivatives: d_torsion / d_r.
static inline void _d_torsion_d_r(double *positions, vector<int> indices, \
	double value, double *deriv)
{
	int i = 0;
	int idx_1 = 0;
	int idx_2 = 0;
	int idx_3 = 0;
	int idx_4 = 0;
	double r_12[3];
	double r_34[3];
	double r_32[3];
	double r_32_len = 0.0;
	double r_12_cross_r_32[3];
	double r_32_cross_r_34[3];
	double r_12_dot_r_32 = 0.0;
	double r_34_dot_r_32 = 0.0;
	double r_12_cross_r_32_norm_2 = 0.0;
	double r_32_cross_r_34_norm_2 = 0.0;

	idx_1 = indices[0] * 3;
	idx_2 = indices[1] * 3;
	idx_3 = indices[2] * 3;
	idx_4 = indices[3] * 3;
	vecsub(&(positions[idx_2]), &(positions[idx_1]), r_12, 3);
	vecsub(&(positions[idx_4]), &(positions[idx_3]), r_34, 3);
	vecsub(&(positions[idx_2]), &(positions[idx_3]), r_32, 3);
	r_32_len = vecdis(&(positions[idx_2]), &(positions[idx_3]), 3);
	veccross3(r_12, r_32, r_12_cross_r_32);
	veccross3(r_32, r_34, r_32_cross_r_34);
	r_12_dot_r_32 = vecdot(r_12, r_32, 3);
	r_34_dot_r_32 = vecdot(r_34, r_32, 3);
	r_12_cross_r_32_norm_2 = veclen(r_12_cross_r_32, 3);
	r_12_cross_r_32_norm_2 *= r_12_cross_r_32_norm_2;
	r_32_cross_r_34_norm_2 = veclen(r_32_cross_r_34, 3);
	r_32_cross_r_34_norm_2 *= r_32_cross_r_34_norm_2;
	for (i = 0; i < 3; i++) {
		deriv[i + 0] = r_32_len / \
			r_12_cross_r_32_norm_2 * r_12_cross_r_32[i];
		deriv[i + 9] = r_32_len / \
			r_32_cross_r_34_norm_2 * r_32_cross_r_34[i];
		deriv[i + 9] = deriv[i + 9] * -1.0;
	}
	for (i = 0; i < 3; i++) {
		deriv[i + 3] = deriv[i + 0] * \
			(r_12_dot_r_32 / r_32_len / r_32_len - 1.0) - \
			(r_34_dot_r_32 / r_32_len / r_32_len) * deriv[i + 9];
		deriv[i + 6] = deriv[i + 9] * \
			(r_34_dot_r_32 / r_32_len / r_32_len - 1.0) - \
			(r_12_dot_r_32 / r_32_len / r_32_len) * deriv[i + 0];
	}
}

RATTLE::RATTLE(void)
{
	n_constraints = 0;
	die_on_fail = true;
	max_cycles = 300;
	positions_t1 = vector<double>(1);
	values.push_back(vector<double>(1));
	values.push_back(vector<double>(1));
	deriv.push_back(vector<vector<double>>(1));
	deriv.push_back(vector<vector<double>>(1));
	positions_t1.clear();
	values[0].clear();
	values[1].clear();
	deriv[0].clear();
	deriv[1].clear();
}

RATTLE::~RATTLE(void)
{
	clear();
}

int const &RATTLE::get_n_constraints(void) const
{
	return n_constraints;
}

int const &RATTLE::get_max_cycles(void) const
{
	return max_cycles;
}

bool const &RATTLE::get_die_on_fail(void) const
{
	return die_on_fail;
}

std::vector<int> const &RATTLE::get_types(void) const
{
	return types;
}

std::vector<double> const &RATTLE::get_targets(void) const
{
	return targets;
}

std::vector<double> const &RATTLE::get_tolerances(void) const
{
	return tolerances;
}

std::vector<std::vector<int>> const &RATTLE::get_indices(void) const
{
	return indices;
}

void RATTLE::set_max_cycles(int v)
{
	max_cycles = v;
}

void RATTLE::set_die_on_fail(bool v)
{
	die_on_fail = v;
}

void RATTLE::check(int tp, std::vector<int> &idx, double tg, double tol)
{
	int flag = 0;

	for (int i = 0; i < NCONSTTYPES; i++) {
		if (tp == i) {
			flag = 1;
			break;
		}
	}
	if (flag == 0) {
		throw out_of_range("Unknown constraint type: " + to_string(tp));
	} else if ((int)idx.size() != n_atoms_req[tp]) {
		throw length_error(to_string(n_atoms_req[tp]) + \
			" atoms are required to define a type " + \
			to_string(tp) + " constraint!");
	}
}

void RATTLE::append(int tp, std::vector<int> &idx, double tg, double tol)
{
	check(tp, idx, tg, tol);
	n_constraints++;
	types.push_back(tp);
	targets.push_back(tg);
	tolerances.push_back(tol);
	indices.push_back(vector<int>(idx.size()));
	indices[(indices.size() - 1)] = idx;
	values[0].push_back(0.0);
	values[1].push_back(0.0);
	deriv[0].push_back(vector<double>(idx.size() * 3));
	deriv[1].push_back(vector<double>(idx.size() * 3));
}

void RATTLE::pop(int idx)
{
	n_constraints--;
	types.erase(types.begin() + idx);
	targets.erase(targets.begin() + idx);
	tolerances.erase(tolerances.begin() + idx);
	indices.erase(indices.begin() + idx);
	values[0].erase(values[0].begin() + idx);
	values[1].erase(values[1].begin() + idx);
	deriv[0].erase(deriv[0].begin() + idx);
	deriv[1].erase(deriv[1].begin() + idx);
}

void RATTLE::clear(void)
{
	n_constraints = 0;
	types.clear();
	targets.clear();
	tolerances.clear();
	indices.clear();
	values[0].clear();
	values[1].clear();
	deriv[0].clear();
	deriv[1].clear();
}

/* calculate values of a series of geometry variables */
void RATTLE::calculate_geometry_variables(double *positions, int t)
{
	#pragma omp parallel for
	for (int i = 0; i < n_constraints; i++) {
		switch (types[i]) {
		case 0:
			values[t][i] = _distance(positions, indices[i]);
			break;
		case 1:
			values[t][i] = _angle(positions, indices[i]);
			break;
		case 2:
			values[t][i] = _torsion(positions, indices[i]);
			break;
		}
	}
}

/* calculate derivatives of a series of geometry variables */
void RATTLE::calculate_geometry_variables_derivatives(double *positions, int t)
{
	#pragma omp parallel for
	for (int i = 0; i < n_constraints; i++) {
		switch (types[i]) {
		case 0:
			_d_distance_d_r(positions, indices[i], values[t][i], \
				deriv[t][i].data());
			break;
		case 1:
			_d_angle_d_r(positions, indices[i], values[t][i], \
				deriv[t][i].data());
			break;
		case 2:
			_d_torsion_d_r(positions, indices[i], values[t][i], \
				deriv[t][i].data());
			break;
		}
	}
}

void RATTLE::rattle_constrain_q(double *positions, double *velocities, \
	double *mass, double dt, int n_atoms)
{
	int flag = 0;
	double gamma = 0.0;
	double epsilon = 0.0;

	if (omp_get_num_threads() > n_constraints)
		omp_set_num_threads(n_constraints);
	// Check position vector size.
	if (positions_t1.size() != (size_t)(n_atoms * 3))
		positions_t1.resize(n_atoms * 3);
	/* calculate values and derivatives of the geometry variables
	** at current time step. */
	calculate_geometry_variables(positions, 0);
	calculate_geometry_variables_derivatives(positions, 0);
	/* perform RATTLE loops part 2 */
	for (int j = 0; j < max_cycles; j++) {
		flag = 0;
		/* update positions by dt */
		#pragma omp parallel for
		for (int i = 0; i < n_constraints; i++) {
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l] * 3;
				positions_t1[idx + 0] = positions[idx + 0] + \
					velocities[idx + 0] * dt;
				positions_t1[idx + 1] = positions[idx + 1] + \
					velocities[idx + 1] * dt;
				positions_t1[idx + 2] = positions[idx + 2] + \
					velocities[idx + 2] * dt;
			}
		}
		/* calculate values and derivatives of the geometry variables
		** at next time step. */
		calculate_geometry_variables(positions_t1.data(),1);
		calculate_geometry_variables_derivatives(positions_t1.data(),1);
		/* loop over each contraint */
		for (int i = 0; i < n_constraints; i++) {
			double tmp = 0.0;
			/* if current contraint satisfies our tolerance,
			** go to next one */
			epsilon = values[1][i] - targets[i];
			if (fabs(epsilon) < tolerances[i])
				goto CONST_I_END;
			else
				flag++;
			/* calculate the Lagrange multiplier */
			gamma = epsilon / dt;
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l];
				tmp += vecdot(&(deriv[0][i].data()[l * 3]), \
					&(deriv[1][i].data()[l * 3]), 3) / \
					(2.0 * mass[idx]);
			}
			gamma /= tmp;
			/* scale the velocities */
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l];
				velocities[idx * 3 + 0] -= 0.5 * gamma * \
					deriv[0][i][l * 3 + 0] / mass[idx];
				velocities[idx * 3 + 1] -= 0.5 * gamma * \
					deriv[0][i][l * 3 + 1] / mass[idx];
				velocities[idx * 3 + 2] -= 0.5 * gamma * \
					deriv[0][i][l * 3 + 2] / mass[idx];
			}
			CONST_I_END:;
		}
		if (flag == 0) {
			break;
		} else if (j == (max_cycles - 1)) {
			string msg = "RATTLE PART 1 DID NOT CONVERGE AFTER " + \
				to_string(max_cycles) + " STEPS!";
			cerr << msg << endl;
			if (die_on_fail)
				throw(runtime_error(msg));
		}
	}
}

void RATTLE::rattle_constrain_p(double *positions, double *velocities, \
	double *mass, double dt, int n_atoms)
{
	int flag = 0;
	double eta = 0.0;

	if (omp_get_num_threads() > n_constraints)
		omp_set_num_threads(n_constraints);
	/* calculate values and derivatives of the geometry variables
	** at next time step. */
	calculate_geometry_variables(positions, 1);
	calculate_geometry_variables_derivatives(positions, 1);
	/* perform RATTLE loops part 2 */
	for (int j = 0; j < max_cycles; j++) {
		flag = 0;
		/* loop over each contraint */
		for (int i = 0; i < n_constraints; i++) {
			eta = 0.0;
			double tmp = 0.0;
			/* if the velocities of current contraint satisfies our
			** tolerance, go to next one */
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l];
				eta += vecdot(&(velocities[idx * 3]), \
					&(deriv[1][i].data()[l * 3]), 3);
			}
			if (fabs(eta) < tolerances[i] / fabs(dt))
				goto CONST_I_END;
			else
				flag++;
			/* calculate the Lagrange multiplier */
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l];
				tmp += vecdot(&(deriv[1][i].data()[l * 3]), \
					&(deriv[1][i].data()[l * 3]), 3) / \
					(2.0 * mass[idx]);
			}
			eta /= tmp;
			/* scale the velocities */
			for (int l = 0; l < (int)(indices[i].size()); l++) {
				int idx = indices[i][l];
				velocities[idx * 3 + 0] -= 0.5 * eta * \
					deriv[1][i][l * 3 + 0] / mass[idx];
				velocities[idx * 3 + 1] -= 0.5 * eta * \
					deriv[1][i][l * 3 + 1] / mass[idx];
				velocities[idx * 3 + 2] -= 0.5 * eta * \
					deriv[1][i][l * 3 + 2] / mass[idx];
			}
			CONST_I_END:;
		}
		if (flag == 0) {
			break;
		} else if (j == (max_cycles - 1)) {
			string msg = "RATTLE PART 2 DID NOT CONVERGE AFTER " + \
				to_string(max_cycles) + " STEPS!";
			cerr << msg << endl;
			if (die_on_fail)
				throw(runtime_error(msg));
		}
	}
}
