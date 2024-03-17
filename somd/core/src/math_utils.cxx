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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* return length of a vector */
double veclen(double *vec, int n_elem)
{
    int i = 0;
    double result = 0.0;

    for (i = 0; i < n_elem; i++)
        result += vec[i] * vec[i];
    return sqrt(result);
}

/* return length^2 of a vector */
double veclen2(double *vec, int n_elem)
{
    int i = 0;
    double result = 0.0;

    for (i = 0; i < n_elem; i++)
        result += vec[i] * vec[i];
    return result;
}

/* normalize a vector */
void vecnorm(double *vec, int n_elem)
{
    int i = 0;
    double len = 0.0;

    len = veclen(vec, n_elem);
    for (i = 0; i < n_elem; i++)
        vec[i] /= len;
}

/* multiply each element of a vector with factor. */
void vecscale(double *vec, double factor, int n_elem)
{
    int i = 0;

    for (i = 0; i < n_elem; i++)
        vec[i] *= factor;
}

/* return v1 dot v2 */
double vecdot(double *vec_1, double *vec_2, int n_elem)
{
    int i = 0;
    double result = 0.0;

    for (i = 0; i < n_elem; i++)
        result += vec_1[i] * vec_2[i];
    return result;
}

/* return ||v2 - v1|| */
double vecdis(double *vec_1, double *vec_2, int n_elem)
{
    int i = 0;
    double result = 0.0;

    for (i = 0; i < n_elem; i++)
        result += (vec_1[i] - vec_2[i]) * (vec_1[i] - vec_2[i]);
    return sqrt(result);
}

/* write (vec_1 + vec_2) to out */
void vecadd(double *vec_1, double *vec_2, double *out, int n_elem)
{
    int i = 0;

    for (i = 0; i < n_elem; i++)
        out[i] = vec_1[i] + vec_2[i];
}

/* write (vec_1 - vec_2) to out */
void vecsub(double *vec_1, double *vec_2, double *out, int n_elem)
{
    int i = 0;

    for (i = 0; i < n_elem; i++)
        out[i] = vec_1[i] - vec_2[i];
}

/* write (v1 x v2) to out */
void veccross3(double vec_1[3], double vec_2[3], double out[3])
{
    out[0] = vec_1[1] * vec_2[2] - vec_2[1] * vec_1[2];
    out[1] = vec_1[2] * vec_2[0] - vec_2[2] * vec_1[0];
    out[2] = vec_1[0] * vec_2[1] - vec_2[0] * vec_1[1];
}
