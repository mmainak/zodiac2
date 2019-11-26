#
# This is the makefile for zodiac.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#

COMPILER = mpif90 #-openmp #-vec-report0 -mcmodel=medium -shared-intel 
COMPOPTS = #-I/u2/wes/bgayen/fftw_2/include -I/opt/intel/Compiler/11.0/083/include
LINKOPTS = -lrfftw -lfftw #-L/u2/wes/bgayen/fftw_2/lib/ -lrfftw -lfftw
DLDFLAGS = #-lm /home/156/mxm156/tecplot10/lib/tecio64.a -lstdc++

#COMPILER = ftn
#COMPOPTS = -I/opt/fftw/2.1.5.1/cnos/include -I/usr/local/usp/PETtools/CE/pkgs/netcdf-3.6.2/include
#LINKOPTS = -L/opt/fftw/2.1.5.1/cnos/lib/ -ldrfftw -ldfftw

PARALLEL = TRUE
LES  = FALSE
DUCT = TRUE
CHAN = FALSE


ifeq ($(LES),TRUE)
LES_CHAN = les_chan.o les_chan_th.o
else
LES_CHAN = no_les.o
endif

#ifeq ($(DUCT),TRUE)
#DUCT_CASE = duct.o
#else
#CHAN_CASE = chan_baines.o
#endif


ifeq ($(PARALLEL),TRUE)
zodiac: zodiac.F90 modules.o ALLOCATION.o periodic.o $(LES_CHAN) \
	duct.o cavity.o fft.o mpi_duct.o wall_model.o dstretch.o dmgd9v.o flow_statistics.o flow_output.o boundary.o \
	grid_def
	$(COMPILER) $(COMPOPTS) zodiac.F90 -o zodiac \
	periodic.o $(LES_CHAN) \
	duct.o cavity.o fft.o mpi_duct.o wall_model.o dstretch.o dmgd9v.o modules.o ALLOCATION.o flow_statistics.o flow_output.o boundary.o $(LINKOPTS) ${DLDFLAGS}
else
zodiac: zodiac.F90 modules.o ALLOCATION.o periodic.o $(LES_CHAN) \
        duct.o cavity.o fft.o mpi_chan_serial.o wall_model.o dstretch.o dmgd9v.o\
         grid_def 
	$(COMPILER) $(COMPOPTS) zodiac.F90 -o zodiac \
        periodic.o $(LES_CHAN) \
	duct.o cavity.o fft.o mpi_chan_serial.o wall_model.o dstretch.o dmgd9v.o modules.o ALLOCATION.o $(LINKOPTS) ${DLDFLAGS}
endif

periodic.o: periodic.f fft.o  grid_def
	$(COMPILER) $(COMPOPTS) -c periodic.f

#channel_baines.o: channel_baines.f fft.o mpi_chan_serial.o wall_model.o header_duct grid_def
#	$(COMPILER) $(COMPOPTS) -c channel_baines.f

ifeq ($(LES),TRUE) 
les_chan.o: les_chan.F90 fft.o  grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan.F90

les_chan_th.o: les_chan_th.F90 fft.o  grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan_th.F90
else
no_les.o: no_les.f
	$(COMPILER) $(COMPOPTS) -c no_les.f
endif

ifeq ($(PARALLEL),TRUE)
mpi_duct.o: mpi_duct.F90  grid_def 
	$(COMPILER) $(COMPOPTS) -c mpi_duct.F90
else
mpi_chan_serial.o: mpi_chan_serial.f 
	$(COMPILER) $(COMPOPTS) -c mpi_chan_serial.f
endif

duct.o: duct.F90 fft.o mpi_duct.o wall_model.o 
	$(COMPILER) $(COMPOPTS) -c duct.F90

cavity.o: cavity.f 
	$(COMPILER) $(COMPOPTS) -c cavity.f

fft.o:  fft.F90 
	$(COMPILER) $(COMPOPTS) -c fft.F90

wall_model.o:   wall_model.f  
	$(COMPILER) $(COMPOPTS) -c wall_model.f

dmgd9v.o: dmgd9v.f  
	$(COMPILER) $(COMPOPTS) -c dmgd9v.f

dstretch.o: dstretch.F90  
	$(COMPILER) $(COMPOPTS) -c dstretch.F90


ALLOCATION.o:   ALLOCATION.F90  modules.o
	$(COMPILER) $(COMPOPTS) -c ALLOCATION.F90

modules.o:      modules.F90
	$(COMPILER) $(COMPOPTS) -c modules.F90

boundary.o:      boundary.F90
	$(COMPILER) $(COMPOPTS) -c boundary.F90

flow_statistics.o:      flow_statistics.F90
	$(COMPILER) $(COMPOPTS) -c flow_statistics.F90
flow_output.o:      flow_output.F90
	$(COMPILER) $(COMPOPTS) -c flow_output.F90

clean:
	rm -f *.o* fort.* *~ zodiac output.txt out_screen.txt  core *.mod

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.

