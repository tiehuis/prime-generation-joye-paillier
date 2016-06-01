#ifndef CHRONO_H
#define CHRONO_H

/* An external library (e.g. GMP) may define this themselves, so do not
 * redefine if already exists, to avoid warnings. */
#ifndef _POSIX_C_SOURCE
#   define _POSIX_C_SOURCE 199309L
#endif

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define DEFAULT_CLOCK CLOCK_MONOTONIC_RAW

#define _1E3 1000ULL
#define _1E6 1000000ULL
#define _1E9 1000000000ULL

typedef struct chrono_struct__ {
    struct timespec begin;
    struct timespec end;
    struct timespec result;
    int clock_type;
} chrono;

static void chrono_diff(struct timespec *r, struct timespec *t, struct timespec *u)
{
    if (u->tv_nsec - t->tv_nsec < 0) {
        r->tv_sec = u->tv_sec - t->tv_sec - 1;
        r->tv_nsec = _1E9 + u->tv_nsec - t->tv_nsec;
    } else {
        r->tv_sec = u->tv_sec - t->tv_sec;
        r->tv_nsec = u->tv_nsec - t->tv_nsec;
    }
}

static inline void chrono_start_with_clock(chrono *t, int clock_type)
{
    t->clock_type = clock_type;
    clock_gettime(clock_type, &t->begin);
}

static inline void chrono_start(chrono *t)
{
    chrono_start_with_clock(t, DEFAULT_CLOCK);
}

static inline void chrono_end(chrono *t)
{
    clock_gettime(t->clock_type, &t->end);
    chrono_diff(&t->result, &t->begin, &t->end);
}

static double chrono_get__(chrono *t, const unsigned int o1, const unsigned int o2)
{
    return (double) t->result.tv_sec * o1 + (double) t->result.tv_nsec / o2;
}

#define chrono_get_sec(t)  chrono_get__(t,    1, _1E9)
#define chrono_get_msec(t) chrono_get__(t, _1E3, _1E6)
#define chrono_get_usec(t) chrono_get__(t, _1E6, _1E3)
#define chrono_get_nsec(t) chrono_get__(t, _1E9, 1)

#endif
