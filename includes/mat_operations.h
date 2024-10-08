#ifndef MAT_OPERATIONS_H
#define MAT_OPERATIONS_H

struct matrix {
    size_t size_x, size_y;
    double* arr;
};

int mat_ctor (matrix* a, const size_t size_x, const size_t size_y);
void mat_dtor (matrix* a);

void mat_print  (const matrix* a);

int mat_resize (matrix* a, const size_t size_x, const size_t size_y);

int mat_swap_cols (matrix* a, const size_t y1, const size_t y2);
int mat_swap_rows (matrix* a, const size_t x1, const size_t x2);
int mat_swap (matrix* dst, matrix* src);
int mat_copy (matrix* dst, const matrix* src);


int mat_const_mult (matrix* a, double c);

int mat_submat (matrix* com, const matrix* a, const size_t i, const size_t j);

int mat_adjugate_mat (matrix* m);
int mat_transp (matrix* a);
int mat_invert (matrix* a);

double mat_minor (const matrix* a, const size_t i, const size_t j);
double mat_adjugate (const matrix* a, const size_t i, const size_t j);

double mat_det (const matrix* a);

int mat_mult (const matrix* a, const matrix* b, matrix* c);
int mat_sum (const matrix* a, const matrix* b, matrix* c);

#endif // MAT_OPERATIONS_H
