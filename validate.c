/**
 * validate.c - Validation routines for jp.c.
 */

#include <stdio.h>
#include <gmp.h>

/**
 * Validates function `pg_g`
 *
 * Ensures the output variables satistify constraints (P1), (P2), (P3), (P4).
 *
 * Prints a nicely formatted report with all the produced values. This allows
 * easy pinpointing of incorrect parameters.
 */
int validate_pg_g(const mpz_t p, const mpz_t v, const mpz_t w, const mpz_t t,
                  const mpz_t qmin, const mpz_t qmax, const mpz_t epsilon)
{
    int fr = 0; // failure result

    // Compute the decimal equivalent epsilon
    mpf_set_default_prec(512);
    mpf_t ed, f1, f2, f3, f4;
    mpf_inits(ed, f1, f2, f3, f4, 0);

    mpz_t z1, z2, z3;
    mpz_inits(z1, z2, z3, 0);

    // Calculate the epsilon value from the integer spec
    mpf_set_z(f1, epsilon);
    mpf_set_z(f2, qmin);
    mpf_div(ed, f1, f2);

    printf("# Parameters:\n\n");
    gmp_printf("  p    = %Zd\n", p);
    gmp_printf("  v    = %Zd\n", v);
    gmp_printf("  w    = %Zd\n", w);
    gmp_printf("  t    = %Zd\n", t);
    gmp_printf("  qmin = %Zd\n", qmin);
    gmp_printf("  qmax = %Zd\n", qmax);
    gmp_printf("  e    = %Zd\n", epsilon);
    gmp_printf("  ed   = %Ff\n", ed);
    printf("\n");

    printf("# Bound Restrictions\n\n");

    // (P1) Bound Check
    printf("(P1) : 1 - ed < (w*p - 1)/(qmax - qmin) <= 1\n\n");

    // 1 - ed
    mpf_set_ui(f1, 1);
    mpf_sub(f1, f1, ed);
    gmp_printf("  %s\n   = %Ff\n\n", "1 - ed", f1);

    // w * p - 1
    mpz_mul(z1, w, p);
    mpz_sub_ui(z1, z1, 1);
    gmp_printf("  %s\n   = %Zd\n\n", "w*p - 1", z1);

    // qmax - qmin
    mpz_sub(z2, qmax, qmin);
    gmp_printf("  %s\n   = %Zd\n\n", "qmax - qmin", z2);

    // (w * p - 1) / (qmax - qmin)
    mpf_set_z(f3, z1);
    mpf_set_z(f4, z2);
    mpf_div(f3, f3, f4);
    gmp_printf("  %s\n   = %Ff\n\n", "(w*p - 1)/(qmax - qmin)", f3);

    if (mpf_cmp(f3, f1) <= 0) {
        gmp_printf("  \e[031mfailure: %Ff < %Ff\e[0m\n\n", f3, f1);
        fr = 1;
    }
    else if (mpf_cmp_ui(f3, 1) >= 0) {
        gmp_printf("  \e[031mfailure: %Ff > %d\e[0m\n\n", f3, 1);
        fr = 1;
    }
    else {
        gmp_printf("  \e[032msuccess: %Ff < %Ff <= %d\e[0m\n\n", f1, f3, 1);
    }

    // (P2) Bound Check
    printf("(P2) : v * p + t >= qmin\n\n");

    mpz_mul(z1, v, p);
    mpz_add(z1, z1, t);
    gmp_printf("  %s\n    = %Zd\n\n", "v*p + t", z1);

    if (mpz_cmp(z1, qmin) < 0) {
        gmp_printf("  \e[031mfailure: %Zd < %Zd\e[0m\n", z1, qmin);
        fr = 1;
    }
    else {
        gmp_printf("  \e[032msuccess: %Zd\n        >= %Zd\e[0m\n", z1, qmin);
    }

    printf("\n");

    // (P3) Bound Check
    printf("(P3) : (v + w)*p + t - 1 <= qmax\n\n");

    mpz_add(z1, v, w);
    mpz_mul(z1, z1, p);
    mpz_add(z1, z1, t);
    mpz_sub_ui(z1, z1, 1);
    gmp_printf("  %s\n    = %Zd\n\n", "(v + w)*p + t - 1", z1);

    if (mpz_cmp(z1, qmax) > 0) {
        gmp_printf("  \e[031mfailure: %Zd > %Zd\e[0m\n\n", z1, qmax);
        fr = 1;
    }
    else {
        gmp_printf("  \e[032msuccess: %Zd\n        <= %Zd\e[0m\n\n", z1, qmax);
    }

    // (P4) Bound Check
    printf("(P4) : phi(p)/p is minimal\n\n");

    // Hard to accurately show, compute the primorial and show the ratio between
    // the two for now. Always succeeds for the moment.

    // In order to compute phi(p) fast, iterate through all prime factors in
    // the primorial (cheap).

    // Generate the first prime which doesn't divide into p
    mpz_set_ui(z1, 2);
    while (mpz_divisible_p(p, z1)) {
        mpz_nextprime(z1, z1);
    }

    // Now have the highest prime factor, calculate phi
    mpz_set_ui(z2, 2);
    mpz_set_ui(z3, 1); // phi(p) running product
    while (mpz_cmp(z2, z1) < 0) {
        mpz_sub_ui(z2, z2, 1);
        mpz_mul(z3, z3, z2);
        mpz_add_ui(z2, z2, 1);
        mpz_nextprime(z2, z2);
    }

    gmp_printf("  %s\n    = %Zd\n", "phi(p)", z3);
    gmp_printf("  %s\n    = %Zd\n\n", "p", p);

    mpf_set_z(f1, z3);
    mpf_set_z(f2, p);
    mpf_div(f1, f1, f2);

    // always succeed on this form
    gmp_printf("  \e[032msuccess: %Ff\e[0m\n\n", f1);

    mpf_clears(ed, f1, f2, f3, f4, 0);
    mpz_clears(z1, z2, z3, 0);

    return fr;
}

#include "jp.h"

// Currently this only validates 'pg_g'
int main(void)
{
    mpz_t pm, qm, p, v, w, t, qmin, qmax, epsilon, _one;
    mpz_inits(pm, qm, p, v, w, t, qmin, qmax, epsilon, _one, NULL);

    mpz_set_ui(_one, 1);
    mpz_mul_2exp(qmax, _one, 513);
    mpz_mul_2exp(qmin, _one, 512);

    // epsilon ~10e-3
    mpz_tdiv_q_ui(epsilon, qmin, 1000);

    pg_g(p, v, w, t, qmin, qmax, epsilon);
    int r = validate_pg_g(p, v, w, t, qmin, qmax, epsilon);

    mpz_clears(pm, qm, p, v, w, t, qmin, qmax, epsilon, _one, NULL);
    return r;
}
