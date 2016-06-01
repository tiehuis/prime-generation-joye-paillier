/**
 * jp.h - joye-paillier prime generation declarations.
 */

#ifndef JPN_H
#define JPN_H

#include <gmp.h>

/**
 * gpg - generic prime generation
 *
 * The generic prime generation algorithm as described by joye-paillier.
 *
 * This requires precomputed values from `pg_g`.
 */
void gpg(mpz_t rop, gmp_randstate_t rs,
         const mpz_t p, const mpz_t v, const mpz_t w, const mpz_t t,
         const mpz_t qmin, const mpz_t qmax);

/**
 * pg_g - parameter generation (general)
 *
 * Generate the parameters required for the general purpose prime
 * generation algorithm.
 */
void pg_g(mpz_t p, mpz_t v, mpz_t w, mpz_t t,
          const mpz_t qmin, const mpz_t qmax, const mpz_t epsilon);

/**
 * cf - carmichael function
 *
 * Computes the carmichael function for relatively smooth numbers.
 * This is equivalent to factoring, so be wary of inputs.
 */
void cf(mpz_t rop, const mpz_t c);

/**
 * ug - unit generation
 *
 * Produce units in (Z/mZ)* relatively quickly.
 */
void ug(mpz_t rop, gmp_randstate_t rs, const mpz_t m, const mpz_t cm);

/**
 * ug_n - unit generation (naive)
 *
 * This is slower, but does not require the computation of the carmichael
 * number of m.
 */
void ug_n(mpz_t rop, gmp_randstate_t rs, const mpz_t m);

#endif
