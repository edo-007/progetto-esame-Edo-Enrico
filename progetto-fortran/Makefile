FLAGS = -g -Wall 
DIR = heat-dat-files


ifeq ($(out), y)
	FLAGS += -D_OUTPUT=1
	
endif

clear-out:
	rm heat-dat-files/Heat2D*

all:
	gfortran $(FLAGS) -c -o lib/constants-module.o lib/constants-module.f90 
	gfortran $(FLAGS) -c -o lib/utils-module.o lib/utils-module.f90 
	gfortran $(FLAGS) -c -o lib/hd2-serial.o -cpp hd2-serial.f90 
	gfortran $(FLAGS) -o hd2-serial lib/hd2-serial.o lib/constants-module.o lib/utils-module.o 