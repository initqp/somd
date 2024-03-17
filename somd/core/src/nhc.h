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

#ifndef _NHC_H_
#define _NHC_H_

#include <vector>

class NHC
{
protected:
    // Chain length.
    int length;
    // Degree of freedom.
    int n_dof;
    // Number of RESPA iterations.
    int n_respa;
    // Chain temperature.
    //   unit : (K)
    double temperature;
    // Coupling time scale.
    //   unit : (ps)
    double tau;
    // The chain "mass".
    //   unit : (kJ/mol*ps*ps)
    //   size : length
    std::vector<double> Q;
    // The "momentum" of each bead.
    //   unit : (kJ/mol*ps)
    //   size : length
    std::vector<double> p;
    // The "position" of each bead.
    //   unit : (1)
    //   size : length
    std::vector<double> q;
    // The sixth-order Suzuki-Yoshida parameters. Comes from Tuckerman.
    const double w_sy[7] = {
        0.784513610477560,
        0.235573213359357,
        -1.17767998417887,
        1.315186320683907,
        -1.17767998417887,
        0.235573213359357,
        0.784513610477560,
    };
    // Set the Q vector from temperature, dof and tau.
    void set_Q(void);
public:
    // The constructors.
    NHC(double T, double t, int n_bead, int n_d, int n_r);
    // The destructor.
    ~NHC(void);
    // Return the chain length.
    int const &get_length(void) const;
    // Return the degree of freedom.
    int const &get_n_dof(void) const;
    // Return number of RESPA.
    int const &get_n_respa(void) const;
    // Return the coupling time scale.
    double const &get_tau(void) const;
    // Return chain temperature.
    double const &get_temperature(void) const;
    // Return chain masses.
    std::vector<double> const &get_Q(void) const;
    // Return chain positions.
    std::vector<double> const &get_q(void) const;
    // Return chain momentums.
    std::vector<double> const &get_p(void) const;
    // Set the degree of freedom.
    void set_n_dof(int n_d);
    // Set number of RESPA.
    void set_n_respa(int n_r);
    // Set the coupling time scale.
    void set_tau(double t);
    // Set the chain temperature.
    void set_temperature(double T);
    // Set the chain positions.
    void set_q(std::vector<double> &v);
    // Set the chain momentums.
    void set_p(std::vector<double> &v);
    // Propagate the chains by dt / 2.
    double propagate(double E_k, double dt);
}; //class NHC

#endif // _NHC_H_
