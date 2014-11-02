#molden2gamess#

This python script is intented to convert cartesian orbitals and basis set specified in [Molden Format](http://www.cmbi.ru.nl/molden/molden_format.html).

As of right now this script works for
* RHF and UHF methods
* Basis sets up to f functions
* Basis sets with SP functions

If the Molden file was generated by Turbomole, add the `-t` flag to scale the d and f orbitals to match what is expected in GAMESS.

Example Molden format files are provided in the `tests` directory.

<b>Reminder:</b> This script cannot set `ISPHER` or `SCFTYP`, you will need to fix those values after conversion. Also, ECP inputs are not included in the Molden format and must be added manually.

Usage:
```
> ./molden2gamess.py [-t | --turbomole] <molden-file> > gamess.inp
```
