# Contributing to fractal-noise

Thanks for your interest. Bug reports, questions and pull requests are all welcome.

## Reporting a problem or asking a question

Open an issue at <https://github.com/dagmarfraser/fractal-noise/issues>. Please include:

- your MATLAB release and operating system, and which toolboxes are installed
  (`ver` output is ideal);
- a minimal, runnable example (fix the seed with `rng(...)`), what you expected,
  and what happened instead, including the full error message.

Questions about which estimator or generator suits your data are fine as issues too.

## Proposing a change

1. Fork the repository and create a branch from `main`.
2. Make the change. Keep functions self-contained (the toolbox has no
   dependency on outside research code) and keep name-value options backward
   compatible: a new behaviour should be opt-in, with the old behaviour as default.
3. Add or update tests in `tests/` and run the full suite:

   ```matlab
   addpath("toolbox"); results = runtests("tests"); table(results)
   ```

   Tests that need optional software (FieldTrip, DSP System Toolbox, the
   original research code) skip themselves when it is absent.
4. Update the help text, `toolbox/functionSignatures.json` if a signature
   changed, and `CHANGELOG.md`.
5. Open a pull request describing what changed and why. CI runs the test suite
   on every pull request.

## Error-handling convention

Fail loudly rather than return a plausible-looking result: invalid input, an
optimiser stuck on its bound or an empty series should raise an error with an
identifier (`functionName:reason`). Any fallback must be visible (a warning or a
flag in the output), never silent.

## Releases

Releases are tagged `vX.Y.Z` on `main` (semantic versioning). Note that every
push to `main` also publishes a new version on MATLAB File Exchange, so changes
land on a branch first and reach `main` through a pull request.

## Support and governance

fractal-noise is maintained by Dagmar Scott Fraser (University of Birmingham).
Issues are usually answered within two weeks. Decisions on scope and API are
made by the maintainer, in the open on issues and pull requests.
