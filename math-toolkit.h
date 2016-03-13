#ifndef __RAY_MATH_TOOLKIT_H
#define __RAY_MATH_TOOLKIT_H

#include <math.h>
#include <stdio.h>
#include <assert.h>

#ifdef _SIMD_
#include <immintrin.h>  // AVX
#define ALIGN __attribute__ ((aligned (32)))
#endif

#ifdef _MACRO_
static inline
void normalize(double *v)
{
    double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    assert(d != 0.0 && "Error calculating normal");

    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}

#define length(v)   \
    sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

#define add_vector(a, b, out) \
{   \
    out[0] = a[0] + b[0];    \
    out[1] = a[1] + b[1];    \
    out[2] = a[2] + b[2];   \
}

#define subtract_vector(a, b, out) \
{   \
    out[0] = a[0] - b[0];    \
    out[1] = a[1] - b[1];    \
    out[2] = a[2] - b[2];   \
}

#define multiply_vectors(a, b, out) \
{   \
    out[0] = a[0] * b[0];    \
    out[1] = a[1] * b[1];    \
    out[2] = a[2] * b[2];   \
}

#define multiply_vector(a, b, out) \
{   \
    out[0] = a[0] * ((double)b);    \
    out[1] = a[1] * ((double)b);    \
    out[2] = a[2] * ((double)b);   \
}

#define cross_product(v1, v2, out) \
{   \
    out[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
    out[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
    out[2] = v1[0] * v2[1] - v1[1] * v2[0]; \
}

#define dot_product(v1, v2) \
    ((v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]))
#else
static inline
void normalize(double *v)
{
    double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    assert(d != 0.0 && "Error calculating normal");

    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}

static inline
double length(const double *v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline
void add_vector(const double *a, const double *b, double *out)
{
#ifdef _UNROLLING_
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
#elif _SIMD_
    double ALIGN c[4];
    __m256d a1 = _mm256_set_pd(a[0], a[1], a[2], 0.0);
    __m256d b1 = _mm256_set_pd(b[0], b[1], b[2], 0.0);
    __m256d t1 = _mm256_add_pd(a1, b1);
    _mm256_store_pd(c, t1);
    out[0]=c[3];
    out[1]=c[2];
    out[2]=c[1];
#else
    for (int i = 0; i < 3; i++)
        out[i] = a[i] + b[i];
#endif
}

static inline
void subtract_vector(const double *a, const double *b, double *out)
{
#ifdef _UNROLLING_
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
#elif _SIMD_
    double ALIGN c[4];
    __m256d a1 = _mm256_set_pd(a[0], a[1], a[2], 0.0);
    __m256d b1 = _mm256_set_pd(b[0], b[1], b[2], 0.0);
    __m256d t1 = _mm256_sub_pd(a1, b1);
    _mm256_store_pd(c, t1);
    out[0]=c[3];
    out[1]=c[2];
    out[2]=c[1];
#else
    for (int i = 0; i < 3; i++)
        out[i] = a[i] - b[i];
#endif
}

static inline
void multiply_vectors(const double *a, const double *b, double *out)
{
#ifdef _UNROLLING_
    out[0] = a[0] * b[0];
    out[1] = a[1] * b[1];
    out[2] = a[2] * b[2];
#elif _SIMD_
    double ALIGN c[4];
    __m256d a1 = _mm256_set_pd(a[0], a[1], a[2], 0.0);
    __m256d b1 = _mm256_set_pd(b[0], b[1], b[2], 0.0);
    __m256d t1 = _mm256_mul_pd(a1, b1);
    _mm256_store_pd(c, t1);
    out[0]=c[3];
    out[1]=c[2];
    out[2]=c[1];
#else
    for (int i = 0; i < 3; i++)
        out[i] = a[i] * b[i];
#endif
}

static inline
void multiply_vector(const double *a, double b, double *out)
{
#ifdef _UNROLLING_
    out[0] = a[0] * b;
    out[1] = a[1] * b;
    out[2] = a[2] * b;
#elif _SIMD_
    double ALIGN c[4];
    __m256d a1 = _mm256_set_pd(a[0], a[1], a[2], 0.0);
    __m256d b1 = _mm256_set1_pd(b);
    __m256d t1 = _mm256_mul_pd(a1, b1);
    _mm256_store_pd(c, t1);
    out[0]=c[3];
    out[1]=c[2];
    out[2]=c[1];
#else
    for (int i = 0; i < 3; i++)
        out[i] = a[i] * b;
#endif
}

static inline
void cross_product(const double *v1, const double *v2, double *out)
{
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static inline
double dot_product(const double *v1, const double *v2)
{
#if _UNROLLING_
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
#elif _SIMD_
//    double ALIGN out[4];
    __m256d a1 = _mm256_set_pd(v1[0], v1[1], v1[2], 0.0);
    __m256d b1 = _mm256_set_pd(v2[0], v2[1], v2[2], 0.0);
    __m256d t1 = _mm256_mul_pd(a1, b1);
//    _mm256_store_pd(out, t1);
//    return out[1]+out[2]+out[3];
    return t1[1]+t1[2]+t1[3];
#else
    double dp = 0.0;

    for (int i = 0; i < 3; i++)
        dp += v1[i] * v2[i];
    return dp;
#endif
}
#endif

static inline
void scalar_triple_product(const double *u, const double *v, const double *w,
                           double *out)
{
    cross_product(v, w, out);
    multiply_vectors(u, out, out);
}

static inline
double scalar_triple(const double *u, const double *v, const double *w)
{
    double tmp[3];
    cross_product(w, u, tmp);
    return dot_product(v, tmp);
}

#endif
