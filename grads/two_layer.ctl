dset ^two_layer.dat
options little_endian
undef 9.999E+20
title /Data/data5/ncep-ncar/prs/prs.grib.mean.y58
xdef 144 linear   0 2.5
ydef  1 linear 0  2.5
tdef 744 linear 1jan2000 1hr
zdef 1 linear 1 1
vars 5
vort1 0 99 vorticity
vort3 0 99 vorticity
psi1 0 99 psi
psi3 0 99 psi
w 0 99 w
ENDVARS
