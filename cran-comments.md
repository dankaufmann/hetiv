## Test environments

* Local macOS 26.5.1, R 4.4.0, Pandoc 3.4 from Quarto
* GitHub Actions macOS, R release: pending branch/PR run
* GitHub Actions Windows, R release: pending branch/PR run
* GitHub Actions Ubuntu, R release: pending branch/PR run

## R CMD check results

Local full check:

* 0 errors | 0 warnings | 3 notes

The local notes are:

* New submission.
* Unable to verify current time.
* HTML manual validation warnings from the local HTML validator.

The full local build rebuilt vignettes successfully and `R CMD check --as-cran`
reported `inst/doc`, vignette files, package vignettes, and vignette output
rebuilding as OK.

## Notes

This is the first CRAN submission of hetiv.
