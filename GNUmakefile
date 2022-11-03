#============================================================================== 
# Copyright (c) Wenjun Ying, School of Mathematical Sciences & Institute of
# Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.
#
# GNU makefile for C/C++/FORTRAN Programming.
#
#  $@ = name of current target
#  $? = list of dependencies newer than the target
#  $< = name of the dependency file, as if selected by make 
#       for use with an implicit rule
#  $^ = the dependency list
#  $* = the basename of the current target
#  $% = the name of the library member being processed
#
#  $(@F) = filename part of a target
#  $(@D) = directory part of a target (= . if no / in value)
#============================================================================== 

PROGRAMS = 
PROGRAMS += solveFitzHughNagumoModel2d 
PROGRAMS += solveGrayScottEqn2d 
PROGRAMS += solveWaveEqn2d 
PROGRAMS += solveFisherEqn2d 
PROGRAMS += solveDiffusionEqn2d2 

OBJECTS = 
#OBJECTS += NewMathGL.o MathGLUT2d.o gl2ps.o 
OBJECTS += QuadraticInterpolation.o Interpolation.o 
OBJECTS += ParametricCurve.o 
OBJECTS += ClosedSplineCurve.o 
OBJECTS += CartesianGrid.o 
OBJECTS += DiffusionEqn.o FisherEqn.o GrayScottEqn.o FitzHughNagumo.o 
OBJECTS += WaveEqn.o 
OBJECTS += DiffusionSolver1d-001-O4.o 
OBJECTS += ReacDiffSolver1d3-2.o 
OBJECTS += ReacDiffSolver1dN-D.o 
OBJECTS += ReacDiffSolver1dD-D-D.o 

DIRS = tool

#============================================================================== 

CC = g++ 

INCFLAGS := $(foreach DIR_NAME, $(DIRS), -I$(DIR_NAME))

CPPFLAGS = -O3 -MMD -DNDEBUG
CPPFLAGS += -fopenmp -DUSE_OPENMP
#CPPFLAGS += -DUSE_OPENGL -DUSE_GL2PS

LIBS =
LIBS += -fopenmp
#LIBS += -lGL -lGLU -lglut 
#LIBS += $(HOME)/lib/libeps3d.so

#============================================================================== 

ifneq ($(PROGRAMS), )
programs:
	@for PROGRAM_NAME in $(PROGRAMS) ; do \
	   $(MAKE) $$PROGRAM_NAME PROGRAM_NAME=$$PROGRAM_NAME; \
	 done
endif

ifneq ($(PROGRAM_NAME), )
$(PROGRAM_NAME) : $(PROGRAM_NAME).o $(OBJECTS) 
	    $(CC) -o $(PROGRAM_NAME) $(PROGRAM_NAME).o $(OBJECTS) $(LIBS)
endif

#============================================================================== 

%.o : %.c 
	$(CC) $(CPPFLAGS) $(INCFLAGS) -c -o $*.o $< 

%.O : %.C 
	$(CC) $(CPPFLAGS) $(INCFLAGS) -c -o $*.o $<

#============================================================================== 

BASE_NAME = $(shell basename $(shell pwd))

clean:
	@rm -f *~ 

save: clean
	@if test ! -e bak; then mkdir bak; fi; 
	@cp -f GNUmakefile *.[HIChc] bak 

zip: cleanall
	@rm -f $(BASE_NAME).zip 
	#@zip $(BASE_NAME).zip -r *
	@zip $(BASE_NAME).zip -r README GNUmakefile *.[HIChcm] $(DIRS)

cleanpics:
	@rm -f *.eps *.eps.gz

cleanobjs:
	@rm -f *~ *.o $(OBJECTS) $(PROGRAMS) 

cleanall:
	@rm -f *~ *.o $(OBJECTS) $(PROGRAMS) 
	@find . -name '*.d' -type f -exec rm -rf {} \;

restart: cleanall
	@$(MAKE) 


