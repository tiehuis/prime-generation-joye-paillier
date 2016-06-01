/**
 * bench.c - benchmarking routines for functions
 *
 * This provides implementation of naive methods and also modified
 * methods of the joye-paillier algorithms which measure various iteration
 * parameters to get a better idea of the overall performance.
 *
 * These could have beebn hidden behind a DEBUG preprocessor macro, but that
 * is a bit messy, so we simply repeat the implementation with minor
 * alterations here.
 */

#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "chrono.h"

/**
 * Unit Generation comparision:
 *
 * For Unit Generation, we are interested in the number of iterations required
 * before we reach a unit value.
 */

int bench_ug_n(mpz_t rop, gmp_randstate_t rs, const mpz_t m)
{
    int iterations = 0;
    mpz_t z1;
    mpz_init(z1);

    do {
        mpz_urandomm(rop, rs, m);
        mpz_gcd(z1, rop, m);
        iterations++;
    } while (mpz_cmp_ui(z1, 1) != 0);

    mpz_clear(z1);
    return iterations;
}

int bench_ug(mpz_t rop, gmp_randstate_t rs, const mpz_t m, const mpz_t cm)
{
    int iterations = 0;
    mpz_t u, r;
    mpz_inits(u, r, 0);

    // Generate in range [0, m - 1].
    mpz_urandomm(rop, rs, m);
    // Transform to [1, m]
    mpz_add_ui(rop, rop, 1);

    while (1) {
        // TODO: Investigate repeated k values being used, is this just
        // a repurcussion of a bad r value being chosen?
        // (1 - k^cm) mod m == (1 - (k^cm mod m)) mod m
        mpz_powm(u, rop, cm, m);
        mpz_ui_sub(u, 1, u);
        mpz_mod(u, u, m);
        iterations++;

        // If u == 0 we are done!
        if (mpz_cmp_ui(u, 0) == 0) {
            break;
        }

        // r in [1, m]
        mpz_urandomm(r, rs, m);
        mpz_add_ui(r, r, 1);

        // k = k + ru mod m
        mpz_mul(r, r, u);
        mpz_add(rop, rop, r);
        mpz_mod(rop, rop, m);
    }

    mpz_clears(u, r, 0);
    return iterations;
}

void cf(mpz_t rop, const mpz_t c)
{
    mpz_t r, z1;
    mpz_inits(r, z1, 0);

    mpz_set(r, c);
    mpz_set_ui(z1, 2);
    mpz_set_ui(rop, 1); // carmichael(p) running product

    while (mpz_cmp_ui(r, 1) != 0) {
        if (mpz_divisible_p(r, z1)) {
            // Reduce r by z1. We do not have any prime powers.
            mpz_divexact(r, r, z1);

            // Found a prime factor, compute carmichael of the factor.
            // This is always (p - 1).
            mpz_sub_ui(z1, z1, 1);
            mpz_lcm(rop, rop, z1);
            mpz_add_ui(z1, z1, 1);
        }

        mpz_nextprime(z1, z1);
    }

    mpz_clears(r, z1, 0);
}

/**
 * Prime Generation comparison:
 *
 * We compare three different methods for their approximate performance and
 * how many values are tested per prime.
 */
int bench_pg_g(mpz_t p, mpz_t v, mpz_t w, mpz_t t,                      // output
         const mpz_t qmin, const mpz_t qmax, const mpz_t epsilon) // input
{
    int iterations = 0;
    mpz_sub(t, qmax, qmin);
    mpz_add_ui(t, t, 1);

    mpz_set_ui(p, 2); // current primorial
    mpz_set_ui(v, 2); // current prime

    while (1) {
        // terminating condition: minimal remainder is beyond epsilon
        {
            // Can we fix w to be fairly smooth?
            mpz_tdiv_r(w, t, p);
            if (mpz_cmp(w, epsilon) >= 0) {
                // remainder exceeds bound, revert last prime multiply
                mpz_divexact(p, p, v);
                mpz_tdiv_q(w, t, p);
                break;
            }
        }

        iterations++;
        mpz_nextprime(v, v);
        mpz_mul(p, p, v);
    }

    // p and w are now fixed
    mpz_tdiv_qr(v, t, qmin, p);
    return iterations;
}

/*
int bench_pg_f(mpz_t p, mpz_t bmin, mpz_t bmax, mpz_t v,          // output
         const mpz_t qmin, const mpz_t qmax, const mpz_t epsilon) // input
{
    mpz_t w, z1, z2;
    mpz_init(w);
    mpz_init(z1);
    mpz_init(z2);

    int iterations = 0;
    mpz_sub(t, qmax, qmin);
    mpz_add_ui(t, t, 1);

    mpz_set_ui(p, 3); // current primorial
    mpz_set_ui(v, 3); // current prime

    while (1) {
        // terminating condition: minimal remainder is beyond epsilon
        {
            mpz_tdiv_r(w, t, p);
            if (mpz_cmp(w, epsilon) >= 0) {
                // remainder exceeds bound, revert last prime multiply
                mpz_divexact(p, p, v);
                mpz_tdiv_q(w, t, p);
                break;
            }
        }

        iterations++;
        mpz_nextprime(v, v);
        mpz_mul(p, p, v);
    }

    // We have calculated w the same method as the naive method, but now vary
    // such that w = bmax - bmin and the values chosen are optimal in regards
    // to (P2) and (P3).

    // (P2) v*p + bmin*p >= qmin
    // (P3) (v + 1)*p + bmax*p - 1 <= qmax

    // Determine value of bmax if v is 0 as a upper bound.
    mpz_cdiv_q(z1, qmin, p);

    // Determine value of bmin if v is 0 as a lower bound.
    mpz_add_ui(z2, qmax, 1);
    mpz_tdiv_q(z2, z2, p);

    // bmax in [z2 + w, z1]
    mpz_add_ui(z2, z2, w);
    mpz_sub(z1, z1, z2);
    mpz_urandomm(bmax, rs, z1);
    mpz_add(bmax, bmax, z2);

    // generate optimal v which satisfies (P3)
    mpz_add_ui(z1, qmax, 1);
    mpz_tdiv_q(z1, z1, p);
    mpz_mul(z2, bmax, p);
    mpz_sub(z1, z1, z2);
    mpz_tdiv_q(z1, z1, p);
    mpz_sub_ui(v, z1, p);

    // Confirm other bound

    // Vary v to satisfy constraints

    mpz_clear(w);
    mpz_clear(z1);
    mpz_clear(z2);
    return iterations;
}
*/

int bench_gpg(mpz_t rop, gmp_randstate_t rs,
              const mpz_t p, const mpz_t v, const mpz_t w, const mpz_t t)
{
    int iterations = 0;
    mpz_t l, m, k, a;
    mpz_init(l);
    mpz_init(m);
    mpz_init(k);
    mpz_init(a);

    mpz_mul(l, v, p);
    mpz_mul(m, w, p);

    // a in (Z/mZ)* \ {1}
    do { bench_ug_n(a, rs, m); } while (mpz_cmp_ui(a, 1) == 0);

    bench_ug_n(k, rs, m);

    while (1) {
        // q <- [(k - t) mod m] + t + l
        mpz_sub(rop, k, t);
        mpz_mod(rop, rop, m);
        mpz_add(rop, rop, t);
        mpz_add(rop, rop, l);
        iterations++;

        // Found a prime
        if (mpz_probab_prime_p(rop, 25) != 0) {
            break;
        }

        mpz_mul(k, a, k);
        mpz_mod(k, k, m);
    }

    mpz_clear(l);
    mpz_clear(m);
    mpz_clear(k);
    mpz_clear(a);
    return iterations;
}

int bench_gpg_f(mpz_t rop, gmp_randstate_t rs,
                const mpz_t p, const mpz_t bmin, const mpz_t bmax, const mpz_t v)
{
    int iterations = 0;
    mpz_t l, b, k, t;
    mpz_init(l);
    mpz_init(b);
    mpz_init(k);
    mpz_init(t);

    mpz_mul(l, v, p);

    // k in (Z/pZ)*
    bench_ug_n(k, rs, p);

    // b in [bmin, bmax]
    mpz_sub(b, bmax, bmin);
    mpz_urandomm(b, rs, b);

    // t = b*p
    mpz_mul(t, b, p);

    while (1) {
        // q <- k + t + l
        mpz_add(rop, k, t);
        mpz_add(rop, rop, l);

        // q <- p - k + t + l
        if (mpz_even_p(rop)) {
            mpz_sub(rop, p, rop);
        }

        iterations++;

        // Found a prime
        if (mpz_probab_prime_p(rop, 25) != 0) {
            break;
        }

        // TODO: mpn_lshift instead
        mpz_mul_ui(k, k, 2);
        mpz_mod(k, k, p);
    }

    mpz_clear(l);
    mpz_clear(b);
    mpz_clear(k);
    mpz_clear(t);
    return iterations;
}

int bench_gpg_n(mpz_t rop, gmp_randstate_t rs, const mpz_t qmin, const mpz_t qmax)
{
    int iterations = 0;
    mpz_t z1;
    mpz_init(z1);
    mpz_sub(z1, qmax, qmin);

    do {
        mpz_urandomm(rop, rs, z1);
        mpz_add(rop, rop, qmin);
        iterations++;
    } while (mpz_probab_prime_p(rop, 25) == 0);

    mpz_clear(z1);
    return iterations;
}

int bench_gspg_n(mpz_t rop, mpz_t qop, gmp_randstate_t rs, const mpz_t qmin, const mpz_t qmax)
{
    int iterations = 0;
    mpz_t z1;
    mpz_init(z1);
    mpz_sub(z1, qmax, qmin);

    while (1) {
        mpz_urandomm(rop, rs, z1);
        mpz_add(rop, rop, qmin);

        // Ensure p % 4 == 3
        mpz_setbit(rop, 0);
        mpz_setbit(rop, 1);
        iterations++;

        if (mpz_probab_prime_p(rop, 25) != 0) {
            mpz_sub_ui(qop, rop, 1);
            mpz_divexact_ui(qop, qop, 2);

            // Both p and q are prime
            if (mpz_probab_prime_p(qop, 25) != 0) {
                break;
            }
        }
    }

    mpz_clear(z1);
    return iterations;
}

int bench_gspg_f(mpz_t rop, mpz_t qop, gmp_randstate_t rs,
                 const mpz_t p, const mpz_t v, const mpz_t w, const mpz_t t)
{
    int iterations = 0;
    /*
    mpz_t l, m, md, a, u, xi;
    mpz_inits(l, m, md, a, u, xi, NULL);

    // Generate initial values
    mpz_mul(l, v, p);
    mpz_mul(m, w, p);
    mpz_set(md, m);     // How is rho set?

    // a in QR(m)
    // m = w*p, we need to be able to factor w fairly quickly for this to be viable
    bench_qr(a, rs, m);

    // xi in (Z/mZ)*
    bench_ug_n(xi, rs, m);

    // k <- 4*u*xi^2 + 3*md mod m
    mpz_mul(xi, xi, xi);
    mpz_mul(xi, xi, u);
    mpz_mul_ui(xi, xi, 4);
    mpz_mul_ui(k, md, 3);
    mpz_add(k, k, xi);
    mpz_mod(k, k, m);

    while (1) {
        // q <- [(k - t) mod m] + t + l
        mpz_sub(rop, k, t);
        mpz_mod(rop, rop, m);
        mpz_add(rop, rop, t);
        mpz_add(rop, rop, l);
        iterations++;

        // Found a prime
        if (mpz_probab_prime_p(rop, 25) != 0) {
            mpz_sub_ui(qop, rop, 1);
            mpz_divexact(qop, qop, 2);

            // Found a safe prime
            if (mpz_probab_prime_p(qop, 25) != 0) {
                break;
            }
        }

        mpz_mul(k, a, k);
        mpz_mod(k, k, m);
    }

    mpz_clears(l, m, md, a, u, xi, NULL);
    */
    return iterations;
}

/**
 * Reporting code:
 *
 * These simply use the above functions, relaying the output to stdout in a
 * fairly consistent manner.
 */
void bench_unit_generation(const int runs)
{
    chrono tm_naive, tm_fast;
    int it_naive = 0, it_fast = 0;

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, rand());

    mpz_t r, m, cm;
    mpz_inits(r, m, cm, NULL);

    // Generate a value m that is our group. This is a primorial since it
    // usually is expected to be.
    const int primorial_n = rand() % 100 + 30;
    mpz_primorial_ui(m, primorial_n);

    // Naive method
    chrono_start(&tm_naive);
    {
        for (int i = 0; i < runs; ++i) {
            it_naive += bench_ug_n(r, rs, m);
        }
    }
    chrono_end(&tm_naive);

    // Generate the carmichael value of m
    cf(cm, m);

    // Fast method
    chrono_start(&tm_fast);
    {
        for (int i = 0; i < runs; ++i) {
            it_fast += bench_ug(r, rs, m, cm);
        }
    }
    chrono_end(&tm_fast);

    // Report statistics
    printf("# Unit Generation\n");
    printf("  Primorial Value: %d\n", primorial_n);
    printf("  Runs: %d\n\n", runs);
    printf("Naive:\n  %lfms total\n  %lfns per run\n"
           "  total iterations: %d\n  iterations per run: %f\n\n",
           chrono_get_msec(&tm_naive), chrono_get_nsec(&tm_naive) / runs,
           it_naive, (double) it_naive / runs);

    printf("Fast:\n  %lfms total\n  %lfns per run\n"
           "  total iterations: %d\n  iterations per run: %f\n\n",
           chrono_get_msec(&tm_fast), chrono_get_nsec(&tm_fast) / runs,
           it_fast, (double) it_fast / runs);

    gmp_randclear(rs);
    mpz_clears(r, m, cm, NULL);
}

void bench_prime_generation(const int runs, const mp_bitcnt_t bs)
{
    chrono tm_naive, tm_fast;
    int it_naive = 0, it_fast = 0;

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, rand());

    mpz_t p, v, w, t, epsilon; // Fast generation specific
    mpz_t qmin, qmax;          // Generation range
    mpz_t r;                   // result values
    mpz_inits(p, v, w, t, epsilon, qmin, qmax, r, NULL);

    // Compute [qmin, qmax] from bitsize
    mpz_set_ui(r, 1);
    mpz_mul_2exp(qmax, r, bs + 1);
    mpz_mul_2exp(qmin, r, bs);

    // Compute fast generation parameters
    mpz_tdiv_q_ui(epsilon, qmin, 1000); // 10^-3 is a good choice
    bench_pg_g(p, v, w, t, qmin, qmax, epsilon);

    // Naive method
    chrono_start(&tm_naive);
    {
        for (int i = 0; i < runs; ++i) {
            it_naive += bench_gpg_n(r, rs, qmin, qmax);
        }
    }
    chrono_end(&tm_naive);

    // Fast method
    chrono_start(&tm_fast);
    {
        for (int i = 0; i < runs; ++i) {
            it_fast += bench_gpg(r, rs, p, v, w, t);
        }
    }
    chrono_end(&tm_fast);

    // Report statistics
    printf("# Prime Generation\n");
    printf("  Bitsize: %lu\n", bs);
    printf("  Runs: %d\n\n", runs);
    printf("Naive:\n  %lfms total\n  %lfns per run\n"
           "  total tested: %d\n  tested per run: %f\n\n",
           chrono_get_msec(&tm_naive), chrono_get_nsec(&tm_naive) / runs,
           it_naive, (double) it_naive / runs);

    printf("Fast:\n  %lfms total\n  %lfns per run\n"
           "  total tested: %d\n  tested per run: %f\n\n",
           chrono_get_msec(&tm_fast), chrono_get_nsec(&tm_fast) / runs,
           it_fast, (double) it_fast / runs);

    gmp_randclear(rs);
    mpz_clears(p, v, w, t, epsilon, qmin, qmax, r, NULL);
}

void bench_safe_prime_generation(const int runs, const mp_bitcnt_t bs)
{
    chrono tm_naive, tm_fast;
    int it_naive = 0, it_fast = 0;

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, rand());

    mpz_t qmin, qmax;          // Generation range
    mpz_t r1, r2;              // result values
    mpz_inits(qmin, qmax, r1, r2, NULL);

    // Compute [qmin, qmax] from bitsize
    mpz_set_ui(r1, 1);
    mpz_mul_2exp(qmax, r1, bs + 1);
    mpz_mul_2exp(qmin, r1, bs);

    chrono_start(&tm_naive);
    {
        for (int i = 0; i < runs; ++i) {
            it_naive += bench_gspg_n(r1, r2, rs, qmin, qmax);
        }
    }
    chrono_end(&tm_naive);

    // Report statistics
    printf("# Safe Prime Generation\n");
    printf("  Bitsize: %lu\n", bs);
    printf("  Runs: %d\n\n", runs);
    printf("Naive:\n  %lfms total\n  %lfns per run\n"
           "  total tested: %d\n  tested per run: %f\n\n",
           chrono_get_msec(&tm_naive), chrono_get_nsec(&tm_naive) / runs,
           it_naive, (double) it_naive / runs);

    gmp_randclear(rs);
    mpz_clears(qmin, qmax, r1, r2, NULL);
}

int main(void)
{
    srand(time(NULL));
    bench_unit_generation(1000);
    bench_prime_generation(10, 512);
    bench_safe_prime_generation(1, 256);
}
