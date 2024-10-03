#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mat_operations.h"

int main()
{
    matrix a = {}, b = {}, res = {};
    const size_t size_n = 7;

    mat_ctor(&b, size_n, size_n);
    mat_ctor(&a, size_n, size_n);
    mat_ctor(&res, size_n, size_n);

    double d = 1;
    for(size_t i = 0; i < size_n*size_n; i++, d++)
        a.arr[i] = rand() / d;

    mat_copy(&b, &a);
    mat_invert(&b);

    mat_print(&a);
    mat_print(&b);

    mat_mult(&a, &b, &res);
    mat_print(&res);

    mat_dtor(&a);
    mat_dtor(&b);
    mat_dtor(&res);
    return 0;
}
