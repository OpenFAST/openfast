# TurbSim Module
The legacy version of this module and additional documentation are available
at the [NWTC Software Portal](https://nwtc.nrel.gov/TurbSim/).

## Overview
TurbSim is a stochastic, full-field, turbulent-wind simulator primarialy for
use with [InflowWind](https://nwtc.nrel.gov/InflowWind "InflowWind")-based
simulation tools. It uses a statistical model (as opposed to a physics-based
model) to numerically simulate time series of three-component wind-speed
vectors at points in a two-dimensional vertical rectangular grid that is fixed
in space. 

Spectra of velocity components and spatial coherence are defined in the
frequency domain, and an inverse Fourier transform produces time series. The
underlying theory behind this method of simulating time series assumes a
stationary process.
