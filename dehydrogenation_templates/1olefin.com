%nprocshared=12
%mem=24GB
# opt int=grid=ultrafine M06 gen pseudo=read

Title

0 1
Ir    0.064100000000     -0.000200000000      0.283500000000
C    -0.006900000000      1.322800000000      2.076800000000
C    -0.602400000000      2.054900000000      1.037800000000
H    -1.693100000000      2.033100000000      0.986600000000
C    -0.023900000000      3.359300000000      0.547800000000
C    -0.589100000000      3.896500000000     -0.762300000000
H     1.066100000000      3.285400000000      0.491100000000
H    -0.222200000000      4.111500000000      1.332900000000
H     0.100700000000      4.646900000000     -1.180200000000
H    -0.634600000000      3.080500000000     -1.502100000000
H    -0.115200000000      0.805500000000     -1.180800000000
H     0.296300000000     -1.116100000000      1.503300000000
H     0.987100000000      1.616700000000      2.418400000000
H    -0.611600000000      0.857500000000      2.851800000000
C    -1.963300000000      4.538300000000     -0.614200000000
H    -2.651200000000      3.839400000000     -0.114800000000
H    -1.881400000000      5.406200000000      0.058000000000
C    -2.549900000000      4.967000000000     -1.948300000000
H    -3.527400000000      5.448300000000     -1.834400000000
H    -1.886300000000      5.676300000000     -2.459500000000
H    -2.677400000000      4.103800000000     -2.614600000000
H    -2.154100000000     -0.823900000000      0.254200000000
H     2.431400000000      0.058300000000     -0.015300000000
H     0.338600000000     -1.709800000000     -0.938600000000

t[corelesslist] 0
6-31G**
******
Ir 0
SDD
******

Ir 0
SDD



