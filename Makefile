

all:
	f2py --f77exec=gfortran -c -m dispar3_sub dispar3_sub.f

