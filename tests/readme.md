# Test `aries`

In R, source `run-test.r` to run the tests.
The results of the test will be shown
the file `test.html` (and `test.md`)
in the `test-output` folder.

```
aries.dir <- "path/to/aries"
alspac.dir <- "path/to/alspac"
source("run-test.r")
```

*The tests do load and make use of ALSPAC
   data using the `alspac` R package.
   Consequently, only direct ALSPAC users will be
   able to run the tests and will need to
   set directories containing ALSPAC and ARIES data.*

This test requires 5-10 minutes to run.

Testing the `aries.copy.release()` function requires at least 3-4 hours to run
and therefore uses a separate script.

```
aries.dir <- "path/to/aries"
alspac.dir <- "path/to/alspac"
new.dir <- "path/to/new/aries"
source("run-copy-release-test.r")
```




