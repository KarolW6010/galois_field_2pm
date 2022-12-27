# galois_field_2pm
A Rust library for representing and performing arithmetic of elements of Galois Fields of size 2<sup>M</sup>.

This library suppoorts addition, subtraction, multiplication, division, and inversion of elements in GF(2<sup>M</sup>) for 2 ≤ M ≤ 127.
Two implementations are supported. The implementations affect how multiplication, division, and inversion are computed.
The first implementation uses look up tables while the second implementation uses the Extended Euclidean Algorithm.

## Representing Galois Fields
GF(2<sup>M</sup>) is isomorphic to ${GF(2)[x] \over p(x)}$ where p(x) is an irreducible polynomial over GF(2) of degree M.
Elements of GF(2<sup>M</sup>) can therfore be uniquely be mapped to polynomials over GF(2) with degree less than M.
In order to represent all the polynomials over GF(2) with degree less than M, M bits are needed.

## Implementations
The look up table implementation can only be used when p(x) is a primitive polynomial with degree less than or equal to 16.
The other implementation will work with any irreducible polynomial up to degree 127.
