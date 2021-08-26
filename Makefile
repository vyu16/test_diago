#FC = mpif90
#FFLAGS = -fast -tp=pwr9 -cuda -gpu=cc70,cuda11.0 -cudalib=cusolver,cublas
#ELPA_INC = -I/ccs/home/vwzyu/opt/elsi/build/include
#ELPA_LIB = -L/ccs/home/vwzyu/opt/elsi/build/lib -lelpa
#MATH_LIB = -L/ccs/home/vwzyu/opt/scalapack-2.1.0 -lscalapack -L$(OLCF_ESSL_ROOT)/lib64 -lesslsmp -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64 -llapack

SRC = test_elpa.f90 test_scalapack.f90 test_lapack.f90 test_cusolver.f90
EXE = $(SRC:.f90=.x)

all: $(EXE)

%.x: %.f90
	$(FC) $(FFLAGS) -o $@ $< $(ELPA_INC) $(ELPA_LIB) $(MATH_LIB)

clean:
	rm -f *.[ox]
