# ALGCM1D

[![Build Status](https://travis-ci.com/udistr/ALGCM1D.jl.svg?branch=master)](https://travis-ci.com/udistr/ALGCM1D.jl)
[![Coverage](https://codecov.io/gh/udistr/ALGCM1D.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/udistr/ALGCM1D.jl)

This code is under construction

A code for a simple 1D coupled land-atmosphere model.

Prognostic variables:

Atmosphere: Heat, horizontal momentum (1 direction) and moisture.
Ocean: Heat and momentum (salinity not yet implemented)
Grid:

height coordinates
Processess included:

* Vertical advaction, explicit integration
* Vetical diffustion:
  * Holtslag and Boville (1993) for the atmosphere
* Bulk formulea for the land-sea interface
* Solar radiation
* Condensation and deposition of super-saturated water vapor

