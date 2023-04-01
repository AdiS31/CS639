#include "PointwiseOps.h"

// #define DO_NOT_USE_MKL
#ifndef DO_NOT_USE_MKL
#include <mkl.h>
#endif

void Copy(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM])
{
#pragma omp parallel for    
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        y[i][j][k] = x[i][j][k];
}

void CopyMKL(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM])
{
    cblas_scopy (
       (int) XDIM*YDIM*ZDIM,
       x[0][0],
       1,
       y[0][0],
       1);
}
void Saxpy(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const float scale)
{
#pragma omp parallel for
    for (int i = 0; i < XDIM; i++)
    for (int j = 0; j < YDIM; j++)
    for (int k = 0; k < ZDIM; k++)
        z[i][j][k] = x[i][j][k] * scale + y[i][j][k];
}

void SaxpyFull(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const float scale)
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];
    float *tempRaw = new float [XDIM*YDIM*ZDIM];
    array_t temp = reinterpret_cast<array_t>(*tempRaw);

    cblas_scopy (
       XDIM*YDIM*ZDIM,
       y[0][0],
       1,
       temp[0][0],
       1);
    
    cblas_saxpy(
        XDIM * YDIM * ZDIM,
        scale,
        &x[0][0][0],
        1,
        &temp[0][0][0],
        1);
    
    const float *t = temp;

    cblas_scopy(
        XDIM * YDIM * ZDIM,
        t[0][0],
        1,
        z[0][0],
        1);
}
void SaxpyX(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM], float (&copyX)[XDIM][YDIM][ZDIM], const float scale)
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];
    float *tempRaw = new float [XDIM*YDIM*ZDIM];
    array_t temp = reinterpret_cast<array_t>(*tempRaw);

    cblas_scopy (
       XDIM*YDIM*ZDIM,
       y[0][0],
       1,
       temp[0][0],
       1);
    
    cblas_saxpy(
        XDIM * YDIM * ZDIM,
        scale,
        &x[0][0][0],
        1,
        &temp[0][0][0],
        1);
    
    cblas_scopy(
        XDIM * YDIM * ZDIM,
        temp[0][0],
        1,
        copyX[0][0],
        1);
    
    
}

void Saxpy(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM],
    const float scale)
{

#ifdef DO_NOT_USE_MKL    
    // Just for reference -- implementation without MKL
#pragma omp parallel for
    for (int i = 0; i < XDIM; i++)
    for (int j = 0; j < YDIM; j++)
    for (int k = 0; k < ZDIM; k++)
        y[i][j][k] += x[i][j][k] * scale;
#else
    cblas_saxpy(
        XDIM * YDIM * ZDIM, // Length of vectors
        scale,              // Scale factor
        &x[0][0][0],        // Input vector x, in operation y := x * scale + y
        1,                  // Use step 1 for x
        &y[0][0][0],        // Input/output vector y, in operation y := x * scale + y
        1                   // Use step 2 for y
    );
#endif
}
