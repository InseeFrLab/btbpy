## Beyond the Border (Python version) : Kernel Density Estimation for Urban Geography

### btbpy (Beyond the Border for Python Users)

`btbpy` is a partial transposition of R `btb` package, available on the [CRAN](https://cran.r-project.org/web/packages/btb/index.html).

### Contributions

Developer and maintainer of `btbpy` package
Julien Jamme, <julien.jamme@protonmail.com>

Authors and Contributors of R `btb` package:
Arlindo Dos Santos [cre],
Francois Semecurbe [drt, aut],
Auriane Renaud [ctb],
Farida Marouchi [ctb]
Joachim Timoteo [ctb]

#### What do `btbpy` and `btb` do ?

The kernelSmoothing() function allows you to square and smooth geolocated data. It calculates a classical kernel smoothing (conservative) or a geographically weighted median. There are only two major call modes of the function. The smoothing with quantiles method is not available on the btbpy package.
The first call mode is kernelSmoothing(obs, epsg, cellsize, bandwith) for a classical kernel smoothing and automatic grid.
The second call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, centroids) for a classical kernel smoothing and user grid.
        
Geographically weighted summary statistics : a framework for localised exploratory data analysis, C.Brunsdon & al., in Computers, Environment and Urban Systems C.Brunsdon & al. (2002) <doi:10.1016/S0198-9715(01)00009-6>, 
Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, pp. 83-86, (2003) <doi:10.1080/13658816.2014.937718>.
