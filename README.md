# Fast Prime Generation (Joye-Paillier)

The provided code implements the algorithms as described Joye, Paillier in the
paper [Fast Generation of Prime Numbers on Portable Devices: An Update (2006)]
(https://www.iacr.org/archive/ches2006/13/13.pdf).

# Documentation

See the function declarations/implementations for a description of the code and
a rough overview of how it works/design choices.

# Performance

The following is some output produced by `bench.c`. Not all routines are complete,
but we note the following:

### Reduced Iteration Count on Prime Generation
The fast prime generation algorithm only requires 10% of the primality tests on average
compared to the naive method. This is a good sign, however, we note that the performance is
still slightly slower, largely due to the more intensive operations in the inner loop.

This is addressed in the faster variant that is yet ot be implemented, which reduces
the modular exponentiations with simple shifts and the like.

### Unit Generation Performance
This is slightly underwhelming in most cases, and the implementation should be
looked over to verify if it is actually correct. The naive method in most cases is
quicker and it does not require precomputations of carmichael values.

# Benchmark

```
# Unit Generation
  Primorial Value: 46
  Runs: 1000

Naive:
  1.453659ms total
  1453.659000ns per run
  total iterations: 7033
  iterations per run: 7.033000

Fast:
  1.947018ms total
  1947.018000ns per run
  total iterations: 2489
  iterations per run: 2.489000

# Prime Generation
  Bitsize: 1560
  Runs: 100

Naive:
  27814.689353ms total
  278146893.530000ns per run
  total tested: 104549
  tested per run: 1045.490000

Fast:
  30277.630804ms total
  302776308.040000ns per run
  total tested: 9429
  tested per run: 94.290000

# Safe Prime Generation
  Bitsize: 256
  Runs: 1

Naive:
  267.250007ms total
  267250007.000000ns per run
  total tested: 21062
  tested per run: 21062.000000
```

# Todo

This currently only implements the general prime generation algorithm listed in
the paper. There remains a fast variant, and also a safe-prime variant that
would be very useful to implement.

The safe prime generation algorithm is supposedly ~25 times faster than a
naive search. I cannot find any other implementations around, so it could have
some useful practical use if implemented.

Benchmarks should be added to show whether these are indeed worth using in
practice.

# License

All the provided code is MIT licensed.
