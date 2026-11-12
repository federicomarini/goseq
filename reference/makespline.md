# Monotonic Spline

Fits a monotonic cubic spline to the data provided, using the penalized
constrained least squares method from the `mgcv` package.

## Usage

``` r
makespline(x, y, newX = NULL, nKnots = 6, lower_bound = 10^-3)
```

## Arguments

- x:

  The predictor variable.

- y:

  The response variable. Must be the same length as `x`.

- newX:

  The points at which to return the value on the fitted spline. If not
  specified `x` is used.

- nKnots:

  The number of knots to use in fitting the spline.

- lower_bound:

  The spline cannot drop below this value.

## Value

Returns a vector of values containing the value of the fit at each point
`newX`.

## Details

This uses the `pcls` function from the mgcv package to produce the fit.
The monotonicity constraint is enforced using `mono.con` from the same
package. The `lower_bound` argument is only used on the rare occasions
when the fitting function becomes negative or arbitrarily close to zero.
If this does occur `lower_bound` is added everywhere to ensure that no
one length is given essentially infinite weighting.

## References

Package mgcv. In particular this function is a modification of an
example given in the man page for `pcls`.

## Author

Matthew D. Young <myoung@wehi.edu.au>.

## Examples

``` r
y <- c( rbinom(50,p=0.4,size=1), rbinom(50,p=0.6,size=1) )
x <- 1:100
plot(x,y)
p <- makespline(x,y)
#> Warning: initial point very close to some inequality constraints
lines(x,p)

```
