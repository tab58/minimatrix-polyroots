# cubic-roots
Implementation of Kahan's cubic root finder.

## Motivation
According to Kahan, standard cubic root finders based on classical methods are prone to numerical inaccuracies. This is an implementation of a numerically stable method illustrated in Kahan's notes titled "To Solve a Real Cubic Equation".

The notes can be found here: http://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf.

## Usage
The quadratic equation solved is Ax<sup>2</sup> + Bx + C = 0.
The cubic equation solved is Ax<sup>3</sup> + Bx<sup>2</sup> + Cx + D = 0.

- getQuadraticRoots(A, B, C)
  - returns 4 real numbers: A1, B1, A2, B2
  - this corresponds to the (possibly) complex roots A1 + B1 * i and A2 + B2 * i.
- getCubicRoots(A, B, C, D)
  - returns 5 real numbers: A0, A1, B1, A2, B2
  - this corresponds to the 1 real and 2 (possibly) complex roots A0, A1 + B1 * i, and A2 + B2 * i.

## Author
Algorithm and pseudocode by William Kahan.
Implementation in Javascript by Tim Bright

## License
(c) Tim Bright, 2017. MIT License.