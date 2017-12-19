# extremes precip scaling

Translation of Paul O'Gorman's matlab code (http://www.mit.edu/~pog/src/precip_extremes_scaling.m) into python, following O'Gorman and Schneider, PNAS, 106, 14773-14777, 2009. Takes vertical profile of vertical velocity; vertical profile of temperature; vertical profile of pressure levels; and surface pressure, and outputs rainfall in kg/m^2/s. 

### Example
```
import numpy as np
from precip_extremes_scaling import scaling

# generate dummy data for this example
omega = np.linspace(-0.3, 0.2, 10)
temp = np.linspace(270, 350, 10)
plev = np.linspace(100000, 10000, 10)
ps = 100000

precip = scaling(omega, temp, plev, ps)  # multiply by 86400 for mm/day
```
