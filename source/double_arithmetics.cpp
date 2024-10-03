#include <math.h>
#include "double_arithmetics.h"

bool is_zero(const double a)
{
    return (fabs(a) < DOUBLE_DEVIATION);
}

bool are_equal(const double a, const double b)
{
    return is_zero(a-b);
}

double min(const double a, const double b)
{
    return ((a < b)? a : b);
}

double max(const double a, const double b)
{
    return ((a > b)? a : b);
}
