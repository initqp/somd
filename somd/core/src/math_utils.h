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

double veclen(double *vec, int n_elem);
double veclen2(double *vec, int n_elem);
double vecdot(double *vec_1, double *vec_2, int n_elem);
double vecdis(double *vec_1, double *vec_2, int n_elem);
void vecnorm(double *vec, int n_elem);
void vecscale(double *vec, double factor, int n_elem);
void vecadd(double *vec_1, double *vec_2, double *out, int n_elem);
void vecsub(double *vec_1, double *vec_2, double *out, int n_elem);
void veccross3(double vec_1[3], double vec_2[3], double out[3]);
