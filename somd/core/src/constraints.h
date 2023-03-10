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

#ifndef _RATTLE_H_
#define _RATTLE_H_

#include <vector>

// Number of valid types.
#define NCONSTTYPES 3

class RATTLE
{
protected:
	// Number of constraints.
	int n_constraints;
	// Max number of RATTLE iterations.
	int max_cycles;
	// If exits when fails.
	bool die_on_fail;
	// Type of constraints.
	std::vector<int> types;
	// Target value of each constraint.
	std::vector<double> targets;
	// Tolerances of the RATTLE iterations.
	std::vector<double> tolerances;
	// Indices of atoms in each constraint.
	std::vector<std::vector<int>> indices;

	// Name of the constraints.
	const char type_str[3][10] = {"distance", "angle", "torsion"};
	// Required atom numbers of each type of constraints.
	int n_atoms_req[NCONSTTYPES] = {2, 3, 4};
	// Positions at next timestep.
	//    size : <n_atom * 3>
	std::vector<double> positions_t1;
	// Values of CVs at this and next timestep.
	//    size : <2 * <n_constraints>>
	std::vector<std::vector<double>> values;
	// Derivatives about positions of CVs at this and next timestep.
	//    size : <2 * <n_constraints * <n_atom_per_constraint * 3>>>
	std::vector<std::vector<std::vector<double>>> deriv;
	// Calculate CV values at this and next timestep.
	void calculate_geometry_variables(double *positions, int t);
	// Calculate derivatives about positions of CVs at this and next
	// timestep.
	void calculate_geometry_variables_derivatives(double *positions, int t);
public:
	// The constructor.
	RATTLE(void);
	// The destructor.
	~RATTLE(void);
	// Get number of constraints.
	int const &get_n_constraints(void) const;
	// Get max_cycles.
	int const &get_max_cycles(void) const;
	// Get die on file.
	bool const &get_die_on_fail(void) const;
	// Get type of a constraint.
	std::vector<int> const &get_types(void) const;
	// Get target value of a constraint.
	std::vector<double> const &get_targets(void) const;
	// Get RATTLE tolerance of a constraint.
	std::vector<double> const &get_tolerances(void) const;
	// Get atom indices of a constraint.
	std::vector<std::vector<int>> const &get_indices(void) const;
	// Set max_cycles.
	void set_max_cycles(int v);
	// Set die_on_fail
	void set_die_on_fail(bool v);
	// Check if a constraint is valid.
	void check(int tp, std::vector<int> &idx, double tg, double tol);
	// Append one constraint.
	void append(int tp, std::vector<int> &idx, double tg, double tol);
	// Delete one constraint.
	void pop(int idx);
	// Delete all constraints.
	void clear(void);

	// Perform the upper part of RATTLE, with the timestep of dt.
	void rattle_constrain_q(double *positions, double *velo, \
		double *mass, double dt, int n_atoms);
	// Perform the lower part of RATTLE, with the timestep of dt.
	void rattle_constrain_p(double *positions, double *velo, \
		double *mass, double dt, int n_atoms);
}; //class RATTLE

#endif /* _RATTLE_H_ */
