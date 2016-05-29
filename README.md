# Fast Prime Generation (Joye-Paillier)

The provided code implements the algorithms as described Joye, Paillier in the
paper [Fast Generation of Prime Numbers on Portable Devices: An Update (2006)]
(https://www.iacr.org/archive/ches2006/13/13.pdf).

# Documentation

See the function declarations/implementations for a description of the code and
a rough overview of how it works/design choices.

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
