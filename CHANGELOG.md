# Change Log

## [Unreleased](https://github.com/Leeds-monc/monc/tree/HEAD)

[Full Changelog](https://github.com/Leeds-monc/monc/compare/v0.9.0...HEAD)

### New features

- Archer2 compatibility changes (#45) by [Chris Symonds](https://github.com/cemac-css), based on changes by Adrian Hines. Includes changes to the following:
    - Casim core - adds dummy modules for compatibility with newer versions of casim
    - IO bridge adds a new field to configuration data structure which replaces `command_to_send` variable in mpi comms
    - Added socrates locations to keyword fcm conf file
    - Added job handling and compilation scripts for archer2
    - halo_swap_neighbors process modified
    - modified bomex mcf file to use 64x64 grids and l_constant_forcing_theta_height variable

- ARC4 compilation modifications and instructions (#28) by [Chris Symonds](https://github.com/cemac-css)

- Continuous integration testing with Travis (using simplified version of Straka test case) (#9 and #10) by [Leif Denby](https://github.com/leifdenby)


### Bug fixes

**Bugfixes**:

- Fix for piecewise linear 1D interpolation (#35) by [Steven
  Boeing](https://github.com/sjboeing)


## [v0.9.0](https://github.com/Leeds-monc/monc/tree/v0.9.0) (2020-04-29)

First tagged version representing fork from MOSRS repository at revision
7765.
