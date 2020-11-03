# CoDiPack extension

Extension for the [CoDiPack - Code Differentiation Package](https://www.scicomp.uni-kl.de/software/codi/).

|CMake flags|```-DGISMO_WITH_CODIPACK=ON``` (default ```OFF```)|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Build status| [CDash](link) |
|Repository|[gismo/gismo](https://github.com/gismo/gismo)|
|Status|completed|
|Developer|Matthias Möller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|03-11-2020|

***

The CoDiPack extension provides two new scalar types
`codi::RealForward` and `codi::RealReverse` that can be used in place
of `real_t` to enable algorithmic differentiation (AD) in forward and
reverse mode, respectively. This extension requires C++11 or better.