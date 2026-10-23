# Contributing

Contributions are welcome.

- Report bugs and questions on the issue tracker. Include a small data excerpt and the call that failed.
- Pull requests should describe the change and keep `tests/run_checks.R` passing. Add a check for any new behaviour.
- Run the checks locally from the repository root before opening a pull request:

  ```
  Rscript tests/run_checks.R
  ```

- Style: plain base R, namespaced calls to suggested packages (`pkg::fun`) and roxygen comments for exported functions.
- Record user-visible changes in `NEWS.md`.
