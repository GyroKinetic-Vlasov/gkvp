### Fujitsu Fortran Compiler ###
FC = mpifrtpx
FFLAGS = -Kfast,parallel # Optimization
FFLAGS += -X9 # Fortran95
FFLAGS += -Koptmsg=2 -Nlst=t # Optimization report
FFLAGS += -fw # Suppress message
FFLAGS += -Kopenmp #-Nfjomplib # OpenMP
FFLAGS += -mcmodel=large # Static memory larger than 2GB
#FFLAGS += -Haefosux -NRtrap #-O0 # Debug
OPTRPT = 'lst'
#FFLAGS += -Nfjprof # Fujitsu profiler fapp
FFLAGS += -Ksimd_nouse_multiple_structures # Specific option for compiler tcs1.2.26 to avoid slowing down GKV
FFLAGS += -Knosch_pre_ra # Specific option for compiler tcs1.2.26 to avoid slowing down GKV
INC = 
LIB =


PROG = 'gkvp.exe'

SRC = ../src/
MYL = ../lib/

MATH = gkvp_math_portable

FFT = gkvp_fft_fftw
### Usage of FFTW (module load fftw-tune)
ifeq ($(FFT),gkvp_fft_fftw)
  #INC += -I$(FFTW_DIR)/include
  #LIB += -L$(FFTW_DIR)/lib -lfftw3 -lm
  LIB += -lfftw3 -lm
endif

FILEIO=gkvp_fileio_fortran
#FILEIO=gkvp_fileio_netcdf
### Usage of NetCDF (module load netcdf-fortran netcdf-c phdf5)
### NetCDF does not work on the FLOW supercomputer for now, Jan 17 2021
ifeq ($(FILEIO),gkvp_fileio_netcdf)
  #INC += -I$(NETCDF_FORTRAN_DIR)/include -I$(NETCDF_DIR)/include -I$(PHDF5_DIR)/include
  #LIB += -L$(NETCDF_FORTRAN_DIR)/lib -L$(NETCDF_DIR)/lib -L$(PHDF5_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
  LIB += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
endif


gkvp:	$(SRC)gkvp_header.f90\
	$(SRC)gkvp_mpienv.f90\
	$(MYL)$(MATH).f90\
	$(SRC)gkvp_clock.f90\
	$(SRC)$(FILEIO).f90\
	$(SRC)gkvp_intgrl.f90\
	$(SRC)gkvp_tips.f90\
	$(SRC)gkvp_vmecbzx.f90\
	$(SRC)gkvp_igs.f90\
	$(SRC)gkvp_bndry.f90\
	$(SRC)gkvp_colli.f90\
	$(SRC)$(FFT).f90\
	$(SRC)gkvp_fld.f90\
	$(SRC)gkvp_colliimp.f90\
	$(SRC)gkvp_freq.f90\
	$(SRC)gkvp_zfilter.f90\
	$(SRC)gkvp_exb.f90\
	$(SRC)gkvp_trans.f90\
	$(SRC)gkvp_shearflow.f90\
	$(SRC)gkvp_advnc.f90\
	$(SRC)gkvp_dtc.f90\
	$(SRC)gkvp_out.f90\
	$(SRC)gkvp_set.f90\
	$(SRC)gkvp_main.f90

	$(FC) $(FFLAGS) -c $(SRC)gkvp_header.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_mpienv.f90
	$(FC) $(FFLAGS) -c $(MYL)$(MATH).f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_clock.f90
	$(FC) $(FFLAGS) -c $(SRC)$(FILEIO).f90 $(INC)
	$(FC) $(FFLAGS) -c $(SRC)gkvp_intgrl.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_tips.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_vmecbzx.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_igs.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_bndry.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_colli.f90
	$(FC) $(FFLAGS) -c $(SRC)$(FFT).f90 $(INC)
	$(FC) $(FFLAGS) -c $(SRC)gkvp_fld.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_colliimp.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_freq.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_zfilter.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_exb.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_trans.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_shearflow.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_advnc.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_dtc.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_out.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_set.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_main.f90

	$(FC) $(FFLAGS)	\
			gkvp_header.o\
			gkvp_mpienv.o\
			$(MATH).o\
			gkvp_clock.o\
			$(FILEIO).o\
			gkvp_intgrl.o\
			gkvp_tips.o\
			gkvp_vmecbzx.o\
			gkvp_igs.o\
			gkvp_bndry.o\
			gkvp_colli.o\
			$(FFT).o\
			gkvp_fld.o\
			gkvp_colliimp.o\
			gkvp_freq.o\
			gkvp_zfilter.o\
			gkvp_exb.o\
			gkvp_trans.o\
			gkvp_shearflow.o\
			gkvp_advnc.o\
			gkvp_dtc.o\
			gkvp_out.o\
			gkvp_set.o\
                        gkvp_main.o\
			-o $(PROG) $(LIB)

	cp *.o *.mod *.$(OPTRPT) ../src/
	rm -f *.o *.mod *.$(OPTRPT)

clean:
	rm -f ../src/*.o ../src/*.mod ../src/*.$(OPTRPT) ./*.exe ./sub.q.*.o* \
	      ./*.o ./*.mod ./*.$(OPTRPT) ./*namelist.* ./sub.q.*

clear:
	rm -f ./*.o ./*.mod ./*.$(OPTRPT) ./*namelist.* ./sub.q.*
 
