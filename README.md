# minimatrix-polyroots
Implementation of Kahan's polynomial root finders for polynomials up to degree 4.

## Motivation
According to Kahan, standard cubic root finders based on classical methods are prone to numerical inaccuracies. This is an implementation of a numerically stable method illustrated in Kahan's notes titled "To Solve a Real Cubic Equation".

The notes can be found here: http://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf.

Kahan's cubic root finder is based on a quadratic root finder. In like fashion, the quartic root finder is based on Kahan's cubic root finder. The linear root finder is easy to calculate and is added for completeness.

## Usage
The linear equation solved is Ax + B = 0.

The quadratic equation solved is Ax<sup>2</sup> + Bx + C = 0.

The cubic equation solved is Ax<sup>3</sup> + Bx<sup>2</sup> + Cx + D = 0.

The quartic equation solved is Ax<sup>4</sup> + Bx<sup>3</sup> + Cx<sup>2</sup> + Dx + E = 0.

- getLinearRoot(A, B)
- getQuadraticRoots(A, B, C)
- getCubicRoots(A, B, C, D)
- getQuarticRoots(A, B, C, D, E)

The roots returned are in an unsorted array. Since roots can be complex, they are given in a complex form:
```
[
  {
    real: <number>,
    imag: <number>
  }
  ...
]
```

## Author
Algorithm and pseudocode by William Kahan.
Implementation in Javascript by Tim Bright.

## License
(c) Tim Bright, 2017. MIT License.