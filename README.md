# ir3
MATLAB codes for performing iterative refinement in up to 3 
different precisions. 

## Related publications
* E. Carson and N. J. Higham. [A new analysis of iterative refinement and its application 
to accurate solution of ill-conditioned sparse linear systems](http://eprints.ma.man.ac.uk/2537/01/covered/MIMS_ep2017_12.pdf). MIMS EPrint 2017.12.
* E. Carson and N. J. Higham. [Accelerating the solution of linear systems by 
iterative refinement in three precisions](http://eprints.ma.man.ac.uk/2562/01/paper.pdf). MIMS EPrint 2017.24. 

## Included MATLAB files
* **_sir3.m_** is a function that performs LU-based iterative refinement with three precisions.

* **_gmresir3.m_** is a function that performs GMRES-based iterative refinement in three precisions.

* **_gmres_hs.m, gmres_sd.m, and gmres_dq.m_** are functions that run left-preconditioned GMRES using precisions half/single, single/double, and double/quad, resp. Application of the preconditioned coefficient matrix to a vector and the preconditioner to the right-hand-side vector are performed in the higher precision; other computations performed all use the lower precision.  

* **_gmresir_example.m_** is an example script for comparing LU-IR and GMRES-IR (with 2 precisions). The test problems correspond to those used in Figures 5.1 and 5.2 of MIMS EPrint 2017.12.

* **_ir3_example.m_** is an example script for running iterative refinement with 3 precisions. The test problems correspond to those used in Figure 10.1 of MIMS EPrint 2017.#.

* The **_vpa/_** directory contains equivalent versions of these files that work with MATLAB's built-in vpa instead of the Advanpix toolkit. 


## Requirements
* The codes have been developed and tested with MATLAB 2015b.
* This code requires Cleve Laboratory to perform half precision computations and 
Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.
    * If you prefer to not use the Advanpix library, use the code within the vpa directory, which uses MATLAB's vpa to perform higher precision computations. Note that this is much slower than Advanpix. 

## License
See license.txt for licensing information.
