# Filter / Polynomial design using SDP

This MATLAB software package designs one-dimensional filters/polynomials via a semi-definite program formulation, following Bumitrescu's book [1].

The SDP formulation solves the problem without discretizing the frequency / polynomial.
This allows us to specify arbitrarily fine starting points and end points of the bands.

The SDP is solved using CVX. 

## fd
Design a zero-phase order-N filter that minimizes maximum absolute deviation from the specified magnitude response

Inputs:

      N - Filter order.
      a - Frequency band start points. From -pi to pi.
      b - Frequency band end points. From -pi to pi.
      m - Magnitude response at frequency bands.

      a, b, m are length-B vectors where B is the number of bands.
Output:

      x - Order-N filter


## pd
Design a degree-N polynomial that minimizes maximum absolute deviation from the specified magnitude response

Inputs:

      N - Polynomial degree.
      a - Bands start points. From 0 to 1.
      b - Bands end points. From 0 to 1.
      m - Magnitude response at intevals.

      a, b, m are length-B vectors where B is the number of bands.
Output:

      x - Order-N filter


## References
[1]	B. Dumitrescu, Positive Trigonometric Polynomials and Signal Processing Applications. Springer Publishing Company, Incorporated, 2007.