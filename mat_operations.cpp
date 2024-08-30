#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mat_operations.h"
#include "double_arithmetics.h"


bool mat_ctor (matrix* a, const size_t size_x, const size_t size_y)
{
    mat_dtor(a);

    a->arr = (double*)calloc((size_x*size_y), sizeof(*(a->arr)));

    if (a->arr == NULL)
        return 0;

    a->size_x = size_x;
    a->size_y = size_y;

    return 1;
}

void mat_dtor (matrix* a)
{
    free(a->arr);
    a->arr = NULL;
    a->size_x = 0;
    a->size_y = 0;
}

void print_mat (const matrix* a)
{
    for (size_t i = 0; i < a->size_x; i++)
    {
        printf("|");
        for (size_t j = 0; j < a->size_y; j++)
            printf("%*.*lg ", ELEMENTS_LENGTH, ELEMENTS_LENGTH, *(a->arr + i*a->size_y + j));

        printf("|\n");
    }
}

bool mat_resize(matrix* a, const size_t size_x, const size_t size_y)
{
    if (a->size_x == size_x && a->size_y == size_y)
        return 1;

    double* temp = a->arr;

    a->arr = (double*)realloc((void*)a->arr, (size_x*size_y)*sizeof(*(a->arr)));

    if (a->arr == NULL)
    {
        free(temp);
        return 0;
    }

    for (size_t i = 0; i < size_x; i++)
    {
        for (size_t j = 0; j < size_y; j++)
        {
            a->arr[i*size_x + j] = 0;
        }
    }

	a->size_x = size_x;
	a->size_y = size_y;

    return 1;
}

int mat_const_mult (matrix* a, const double c)
{
    for (size_t i = 0; i < a->size_x; i++)
        for (size_t j = 0; j < a->size_y; j++)
            a->arr[i*a->size_y + j] *= c;

    return 1;
}

bool mat_swap_rows (matrix* a, const size_t x1, const size_t x2)
{
    if (x1 >= a->size_x || x2 >= a->size_x)
        return 0;

    for (size_t j = 0; j < a->size_y; j++)
    {
        double temp = 0;

        temp = a->arr[x1*a->size_y + j];
        a->arr[x1*a->size_y + j] = a->arr[x2*a->size_y + j];
        a->arr[x2*a->size_y + j] = temp;
    }
    return 1;
}

bool mat_swap_cols (matrix* a, const size_t y1, const size_t y2)
{
    if (y1 >= a->size_y || y2 >= a->size_y)
        return 0;

    for (size_t i = 0; i < a->size_y; i++)
    {
        double temp = 0;

        temp = a->arr[y1*a->size_y + i];
        a->arr[i*a->size_y + y1] = a->arr[i*a->size_y + y2];
        a->arr[i*a->size_y + y2] = temp;
    }
    return 1;
}

int mat_transp(matrix* a)
{
    matrix buff = {};

    if (!mat_ctor(&buff, a->size_y, a->size_x))
        return -1;

    for (size_t i = 0; i < a->size_x; i++)
    {
        for (size_t j = 0; j < a->size_y; j++)
        {
            buff.arr[j*buff.size_y + i] = a->arr[i*a->size_y + j];
        }
    }

    mat_swap(a, &buff);
    mat_dtor(&buff);

    return 1;
}

int mat_swap (matrix* dst, matrix* src)
{
    size_t temp = 0;

    temp = dst->size_x;
    dst->size_x = src->size_x;
    src->size_x = temp;

    temp = dst->size_y;
    dst->size_y = src->size_y;
    src->size_y = temp;

    double* buff = NULL;

    buff = dst->arr;
    dst->arr = src->arr;
    src->arr = buff;

    return 1;
}

int mat_copy (matrix* dst, const matrix* src)
{
    size_t size_x = src->size_x;
    size_t size_y = src->size_y;

    if (!mat_resize(dst, size_x, size_y))
        return -1;

    for (size_t i = 0; i < size_x; i++)
        for (size_t j = 0; j < size_y; j++)
            dst->arr[i*size_y + j] = src->arr[i*size_y + j];

    return 1;
}

int mat_minmat (matrix* dst, const matrix* a, const size_t i, const size_t j)
{
    if (!mat_resize(dst, a->size_x - 1, a->size_y - 1))
        return -1;

    for (size_t ic = 0; ic < dst->size_x; ic++)
        for (size_t jc = 0; jc < dst->size_y; jc++)
        {
            dst->arr[ic*dst->size_y + jc] = a->arr[a->size_y * (ic >= i ? ic+1 : ic) +
                                                           1 * (jc >= j ? jc+1 : jc)];
        }

    return 1;
}

int mat_adjugate (matrix* a)
{
    if (a->size_x != a->size_y)
        return 0;

    matrix adj = {};
    mat_ctor(&adj, a->size_x, a->size_y);

    for (size_t i = 0; i < a->size_x; i++)
        for (size_t j = 0; j < a->size_y; j++)
            adj.arr[i*adj.size_y + j] = mat_minor(a, j, i); // Adj_ij = M_ji

    mat_swap(a, &adj);

    mat_dtor(&adj);

    return 1;
}

int mat_invert (matrix* a)
{
    double det = mat_det(a);

    if (isnan(det) || fabs(det) < 1e-30 || mat_adjugate(a) == 0)
        return 0;

    mat_const_mult(a, 1/det);

    return 1;
}

double mat_det (const matrix* a)
{
    if  (a->size_x != a->size_y)
        return NAN;

    const size_t n = a->size_x;

    matrix d = {};
    mat_ctor(&d, n, n);

    mat_copy(&d, a);

    double det = 1;
    bool sign = true;

    for(size_t i = 0; i < n; i++)
    {
        size_t j = i;
        for (; j < n; j++)
        {
            if (!is_zero(d.arr[i*n + j]))
                break;
        }

        if (j == n)
        {
            if (i == 0)
                return 0;

            continue;
        }
        if (j != i)
        {
            mat_swap_rows(&d, i, j);
            sign = !sign;
        }

        for (size_t u = i+1; u < n; u++)
        {
            double norm_coeff = d.arr[u*n + i] / d.arr[i*n + i];

            for (size_t v = 0; v < n; v++)
            {
                d.arr[u*n + v] -= d.arr[i*n + v] * norm_coeff;
            }
        }
        det *= d.arr[i*n + i];
    }

    det = (sign ? 1 : -1) * det;

    mat_dtor(&d);
    return det;
}

double mat_minor(const matrix* a, const size_t i, const size_t j)
{
    if  (i >= a->size_x || j >= a->size_y || a->size_x != a->size_y)
        return 0;

    matrix com = {};
    mat_ctor(&com, a->size_x, a->size_y);

    mat_minmat(&com, a, i, j);

    double res = mat_det(&com);

    mat_dtor(&com);

    return res;
}

int mat_sum(const matrix* a, const matrix* b, matrix* c)
{
    if (a->size_x != b->size_x || a->size_y != b->size_y)
        return 0;

    size_t size_x = a->size_x;
    size_t size_y = a->size_y;

    if (!mat_resize(c, size_x, size_y))
    {
        return -1;
    }

	for (size_t i = 0; i < size_x; i++)
	{
		for (size_t j = 0; j < size_y; j++)
        {
			c->arr[i*c->size_y + j] = a->arr[i*size_y + j] + b->arr[i*size_y + j];
        }
	}

    return 1;
}

int mat_mult(const matrix* a, const matrix* b, matrix* c)
{
    if (a->size_y != b->size_x)
        return 0;

    // TODO: add size_t
    size_t size_a = a->size_x,
           size_b = b->size_y,
           size_n = a->size_y;

    if (!mat_resize(c, size_a, size_b))
    {
        return -1;
    }

    for (size_t ia = 0; ia < size_a; ia++)
    {
        for (size_t ib = 0; ib < size_b; ib++)
        {
            double value = 0;

            for (size_t n = 0; n < size_n; n++)
            {
                value += a->arr[ia*size_n + n] * b->arr[n*size_b + ib];
            }
            *(c->arr + ia*c->size_y + ib) = value;
        }
    }

    return 1;
}

