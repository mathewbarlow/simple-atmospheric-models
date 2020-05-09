dset ^forecast_1d.dat
options little_endian
undef 9.999E+20
title CDAS data
xdef 144 linear   0  2.5
ydef  1 levels 45
tdef 4225 linear 01jan2007 1dy
zdef 1 levels 500
vars 4
vort00 0 99 initial vorticity
vortfor 0 99 forecast 24hr vorticity
vort24 0 99 obs 24hr vorticity
vortadvec 0 99 vorticity forecast from pure advection
endvars
