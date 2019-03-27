# TurbSim
A stochastic, full-field, turbulence simulator, primarialy for use with [InflowWind](https://nwtc.nrel.gov/InflowWind "InflowWind")-based simulation tools 

**Authors**: Neil Kelley and [Bonnie Jonkman](mailto:bonnie.jonkman@nrel.gov), NREL

TurbSim is a stochastic, full-field, turbulent-wind simulator. It uses a statistical model (as opposed to a physics-based model) to 
numerically simulate time series of three-component wind-speed vectors at points in a two-dimensional vertical rectangular 
grid that is fixed in space. 

Spectra of velocity components and spatial coherence are defined in the frequency domain, and 
an inverse Fourier transform produces time series. The underlying theory behind this method of 
simulating time series assumes a stationary process.

For more information, please refer to documentation on the [TurbSim web site](https://nwtc.nrel.gov/TurbSim).

