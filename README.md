# superkerr
The superkerr XSPEC model: for continuum fitting in the general Kerr metric.   

## More detailed README
See the pdf file README_superkerr for more detailed instructions. 

## Paper
This repository contains the XSPEC model that was developed in [Mummery, Balbus and Ingram 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3373M/abstract). 

Please cite this paper if you use this model. 

## Code 
* The main program is a fortran90 file ``superkerr.f90``
* Raytracing modules are included in the fortran90 file ``amodules.f90``, and associated ``X.mod`` files. 
* Files to load this into XSPEC ``load.xcm`` and ``lmodel.dat`` is also included. 


## Full model input parameter description

All of the model parameters are listed below, along with their units and limits of validity. 

| Parameter | Units | range |
 --- | --- | ---
| $`M `$ |  Solar masses  | $`0 < M `$ |
| $`a_\star `$ | Dimensionless | $`-\infty < a_\star < \infty`$ |
| $`i `$  | Degrees | $`0 < i < 90`$ |
| $`\dot M`$  | $`L_{\rm edd}/c^2 `$  | $`0 < \dot M`$ |
| $`\delta_{\cal J}`$  | Dimensionless  |  $`0 < \delta_{\cal J} < 1`$ |
| $`{\tt norm}`$  |  1/kpc$`^2`$   | $`0 < {\tt norm}`$ |

## Loading into XSPEC 
* For use in the XSPEC software (https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
* XSPEC requires the files: ``lmodel.dat``, ``load.xcm``, ``superkerr.f90``, ``amodules.f90`` and all ``.mod`` files
* ``superkerr.f90`` contains the actual fortran implementation of the model 
* See https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html for more information on how to load the model into XSPEC
* For any questions about the model, please email andrew.mummery@physics.ox.ac.uk 
