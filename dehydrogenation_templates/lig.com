%nprocshared=12
%mem=24GB
# opt int=grid=ultrafine M06 gen pseudo=read

Title

1 1
Ir    0.000000000000      0.475900000000      0.000000000000
H     0.000000000000     -1.591100000000      0.000000000000
H     2.294500000000      0.163800000000     -0.085400000000
H    -2.294500000000      0.163700000000      0.085400000000

t[corelesslist] 0
6-31G**
******
Ir 0
SDD
******

Ir 0
SDD



