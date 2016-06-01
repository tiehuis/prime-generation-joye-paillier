#include <time.h>
#include <gmp.h>
#include "jp.h"

int main(void)
{
    // Use the standard prime generation algorithm
    mpz_t pm, p, v, w, t, qmin, qmax, epsilon, _one;
    mpz_inits(pm, p, v, w, t, qmin, qmax, epsilon, _one, NULL);

    mpz_set_ui(_one, 1);
    mpz_mul_2exp(qmax, _one, 1025);
    mpz_mul_2exp(qmin, _one, 1024);

    // epsilon ~10e-3
    mpz_tdiv_q_ui(epsilon, qmin, 1000);

    pg_g(p, v, w, t, qmin, qmax, epsilon);

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, time(NULL));

    gpg(pm, rs, p, v, w, t, qmin, qmax);

    gmp_printf("%Zd\n", pm);

    gmp_randclear(rs);
    mpz_clears(pm, p, v, w, t, qmin, qmax, epsilon, _one, NULL);
}
