## Filter design using SDP

Design a minimum energy zero-phase order-N filter subject to bounded constraints

Inputs:

      N - Filter order
      a - Frequency bands start point
      b - Frequency bands end point
      l - Lower bounds
      u - Upper bounds
Output:

      x - Designed order-N filter

Example: Low pass filter with cutoff pi/2 and 1.0% passband ripple

      x = filt_design(20, -pi/2, pi/2, 0.99, 1.01);


## References
[1]	B. Dumitrescu, Positive trigonometric polynomials and signal processing applications. 2007.
