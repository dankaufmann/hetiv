# Contributing to hetiv

Thanks for improving `hetiv`. The package is research software, so the
most helpful contributions are reproducible bug reports, small
regression tests, and documentation fixes.

## Development workflow

1.  Install development dependencies with
    `devtools::install_deps(dependencies = TRUE)`.
2.  Load the package with `devtools::load_all()`.
3.  Run tests with `devtools::test()`.
4.  Check test coverage with
    [`covr::package_coverage()`](http://covr.r-lib.org/reference/package_coverage.md).
5.  Regenerate documentation after roxygen changes with
    `devtools::document()`.
6.  Run `devtools::check()` before release-oriented changes.

Please keep examples deterministic, offline-safe, and reasonably fast.
Use small simulated data in tests unless a larger sample is needed to
reproduce a numerical issue.

## Style

The public API uses mathematical names such as `Psi`, `Phi`, `SigE`, and
`H`. Those names are intentionally preserved. New internal helpers
should prefer clear lower-case names, focused functions, and explicit
error messages.

## Bug reports

For numerical or econometric bugs, include:

- the exact function call;
- the dimensions of all input matrices;
- [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html);
- whether the issue depends on the random seed;
- the expected behavior and the observed behavior.
