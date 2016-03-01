## Filter design using SDP

Design a minimum energy zero-phase order-N filter
such that l < |X(exp(jw))| < u for a < w < b

Inputs:
      N - Filter order
      a - Frequency bands start point
      b - Frequency bands end point
      l - Lower bound
      u - Upper bound
Output:
      x - Designed order-N filter

Example:
      Low pass filter with cutoff pi/2 and 1.0% passband ripple
      x = filt_design(20, -pi/2, pi/2, 0.99, 1.01);


## References
[1]	B. Dumitrescu, Positive trigonometric polynomials and signal processing applications. 2007.