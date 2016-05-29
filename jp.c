/**
 * jp.c - joye-paillier prime generation implementation.
 */

#include "jp.h"

/**
 * pg_g - parameter generation (general)
 *
 * Produce the following parameters (p, v, w, t) given the inputs (bit-length)
 *
 * This is a slightly modified version of the usual algorithm. It assumes that
 * qmin and qmax are are single bit-range (i.e. all 512-bit numbers). This
 * allows a trivial set of near optimal parameters.
 *
 * p, v, w, t satisify the following:
 *
 * (P1) 1 - epsilon < (w*p - 1)/(qmax - qmin) <= 1
 *
 * (P2) v*p + t >= qmin
 *
 * (P3) (v + w)*p + t - 1 <= qmax
 *
 * (P4) phi(p) / p is minimized
 *      where p is a product of primes
 */
void pg_g(mpz_t p, mpz_t v, mpz_t w, mpz_t t,                      // output
          const mpz_t qmin, const mpz_t qmax, const mpz_t epsilon) // input
{
    // (P4)
    //
    // phi(p) / p is smaller if the prime p is smaller. i.e.
    // 2/3 < 10/11 clearly. Therefore, a product of primes of length
    // < n is minimized if it is constructed of the n smallest primes.
    //
    // This product is also known as a primorial. We place this initial
    // restriction on p's form to satisfy (P4).
    //
    // Note: We could possibly choose a better value of p if we didn't
    // take a consecutive run of primes, however this requires more heuristics
    // and would be a fair bit more complicated possibly.
    //
    // (P1)
    //
    // In an effort for simplicity, we simply iterate through successive
    // primorials until we reach one which exceeds the given epsilon, or
    // we reach a certain upper threshold?
    //
    // (P2), (P3)
    //
    // These follow trivially. Given p which satisfies (P1):
    //
    // let v = qmin / p
    // let t = qmin % p

    // Calculate (qmin - qmax + 1). We work with epsilons signifying the
    // difference between w*p and (qmax - qmin + 1) to avoid working with
    // floating points.

    mpz_sub(t, qmax, qmin);
    mpz_add_ui(t, t, 1);

    mpz_set_ui(p, 2); // current primorial
    mpz_set_ui(v, 2); // current prime

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

        mpz_nextprime(v, v);
        mpz_mul(p, p, v);
    }

    // p and w are now fixed
    mpz_tdiv_qr(v, t, qmin, p);
}

/**
 * cf - carmichael function
 *
 * Calculate the carcmichael function of number. This is equivalent to
 * factoring, so consider the smoothness of the value `c`.
 */
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
 * ug - unit generation
 *
 * A relatively fast function for generating units (relatively prime values
 * mod m).
 *
 * This requires computation of carmichael value of a number. This is slightly
 * expensive and as a result, this function needs to be called a fair amount
 * (~100-200 times for a ~200-smooth number) before this is quicker than the
 * naive random search for units.
 */
void ug(mpz_t rop, gmp_randstate_t rs, const mpz_t m, const mpz_t cm)
{
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
}

/**
 * ug_n - unit generation (naive)
 *
 * A random search through units. This is sometimes quicker than the
 * faster variant if the group is hard to factor and/or we only require
 * a small number of units.
 */
void ug_n(mpz_t rop, gmp_randstate_t rs, const mpz_t m)
{
    mpz_t z1;
    mpz_init(z1);

    do {
        mpz_urandomm(rop, rs, m);
        mpz_gcd(z1, rop, m);
    } while (mpz_cmp_ui(z1, 1) != 0);

    mpz_clear(z1);
}

/**
 * gpg - generic prime generation
 *
 * The generic prime generation algorithm as described by joye-paillier.
 *
 * This requires precomputed values from `pg_g`.
 */
void gpg(mpz_t rop, gmp_randstate_t rs,
         const mpz_t p, const mpz_t v, const mpz_t w, const mpz_t t,
         const mpz_t qmin, const mpz_t qmax)
{
    mpz_t l, m, k, a;
    mpz_init(l);
    mpz_init(m);
    mpz_init(k);
    mpz_init(a);

    mpz_mul(l, v, p);
    mpz_mul(m, w, p);

    // a in (Z/mZ)* \ {1}
    do { ug_n(a, rs, m); } while (mpz_cmp_ui(a, 1) == 0);

    ug_n(k, rs, m);

    while (1) {
        // q <- [(k - t) mod m] + t + l
        mpz_sub(rop, k, t);
        mpz_mod(rop, rop, m);
        mpz_add(rop, rop, t);
        mpz_add(rop, rop, l);

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
}
