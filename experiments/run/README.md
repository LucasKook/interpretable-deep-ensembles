
Training models
===============

Models can be run from the command line using
```bash
Rscript -e <path-to-script.R> <loss> <formula> <model>
```
For instance,
```bash
Rscript -e run-utkface.R "rps" "age_group ~ 1" "ci"
```
fits the complex intercept model for the UTKFace data using the RPS loss.

