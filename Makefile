#################################################################################
#                                                                               #
#        Makefile for Reinitz lab transcription model                           #
#                                                                               #
#        The options below can be different depending on your personal          #
#        installation. To use this file hassle free I reccomend setting         #
#        these options as environmental variable in your own .bashrc            #
#        so that these do not need to be updated every time you update code     #
#                                                                               #
#################################################################################

# commands to link to the xml2 library, used in the current neoParSA
XML_CFLAGS ?= `xml2-config --cflags`
XML_LIBS   ?= `xml2-config --libs`

# directories that must be present for all compiled files. We need the Boost
# libraries, uncompiled, available online, then links to neoParSA

BOOST_DIR  ?=/users/kenneth/Boost
PARSA_ROOT ?=/users/kenneth/Sandbox/git/neoParSA

PARSA_DIR  = $(PARSA_ROOT)/parsa
LDLIBS     += -L$(PARSA_ROOT)/build/lib -lparsa
LIBPARSA   = $(PARSA_ROOT)/build/lib/libparsa.a

# the user can specify which compiler to use by changing CXX
CXX ?= g++

ifdef PARALLEL
	ifeq ($(CXX),icpc)
		PFLAGS = -openmp -DPARALLEL
	else
		PFLAGS = -fopenmp -DPARALLEL
	endif
endif

ifdef DEBUG
  FLAGS  = $(USR_FLAGS) -g -O2 -Wall -Wstrict-aliasing=0 -I$(BOOST_DIR) -I$(PARSA_DIR) $(XML_CFLAGS) $(PFLAGS)
  RFLAGS = -c -DNDEBUG -fPIC -O2 -g
else
	FLAGS  = $(USR_FLAGS) -O3 -I$(BOOST_DIR) -I$(PARSA_DIR) $(XML_CFLAGS) $(PFLAGS)
	RFLAGS = -c -fPIC -O3
endif

SOURCE = src/subgroup.cpp src/competition.cpp src/sequence.cpp src/score.cpp src/coeffects.cpp \
src/cooperativity.cpp src/scalefactor.cpp src/fasta.cpp src/mode.cpp src/pwm.cpp \
src/TF.cpp src/gene.cpp src/nuclei.cpp src/datatable.cpp src/twobit.cpp \
src/parameter.cpp src/organism.cpp src/bindings.cpp \
src/bindingsite.cpp src/distance.cpp src/promoter.cpp \
src/quenching.cpp


OBJECT=$(SOURCE:.cpp=.o)

HEADER=$(SOURCE:.cpp=.h)

all: transcpp scramble unfold test_moves

# the basic binaries required to run and probe the transcription model
transcpp: $(OBJECT:.o=.$(CXX).o) $(LIBPARSA) src/utils.$(CXX).o src/transcpp.o
	$(CXX) $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/transcpp.o $(XML_LIBS) $(PFLAGS) -o transcpp $(LDLIBS) 
	
ptranscpp: $(OBJECT:.o=.$(CXX).o) $(LIBPARSA) src/utils.$(CXX).o src/ptranscpp.o
	$(CXX) $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/ptranscpp.o $(XML_LIBS) $(PFLAGS) -DUSE_BOOST -o ptranscpp $(LDLIBS) 
	
scramble: $(OBJECT:.o=.$(CXX).o) $(LIBPARSA) src/utils.$(CXX).o src/scramble.o
	$(CXX) $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/scramble.o $(XML_LIBS) $(PFLAGS) -o scramble $(LDLIBS) 
	
unfold: $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/unfold.o
	$(CXX) $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/unfold.o $(XML_LIBS) $(PFLAGS) -o unfold $(LDLIBS) 

test_moves: $(OBJECT:.o=.$(CXX).o) $(LIBPARSA) src/utils.$(CXX).o src/test_moves.o
	$(CXX) $(OBJECT:.o=.$(CXX).o) src/utils.$(CXX).o src/test_moves.o $(XML_LIBS) $(PFLAGS) -o test_moves $(LDLIBS) 

src/transcpp.o: src/transcpp.cpp
	$(CXX) -c $(FLAGS) src/transcpp.cpp -o src/transcpp.o
	
src/ptranscpp.o: src/ptranscpp.cpp
	$(CXX) -c $(FLAGS) src/ptranscpp.cpp -DUSE_BOOST -o src/ptranscpp.o
	
src/scramble.o: src/scramble.cpp
	$(CXX) -c $(FLAGS) src/scramble.cpp -o src/scramble.o
	
src/unfold.o: src/unfold.cpp
	$(CXX) -c $(FLAGS) src/unfold.cpp -o src/unfold.o
	
src/test_moves.o: src/test_moves.cpp
	$(CXX) -c $(FLAGS) src/test_moves.cpp -o src/test_moves.o

# The default way to compile object files
$(OBJECT:.o=.$(CXX).o) : $(SOURCE) $(HEADER)
	$(CXX) -c -fPIC $(FLAGS) $(@:.$(CXX).o=.cpp) -o $@
	
# the utils file needs to be compiled separately so that error and print
# messages can be passed to R or matlab if used
src/utils.$(CXX).o: src/utils.cpp $(HEADER)
	$(CXX) -c -fPIC $(FLAGS) src/utils.cpp -o src/utils.$(CXX).o
	

#################################################################################

# compilation of the R interface

#################################################################################

# if compiling the R interface you must have Rcpp installed. You can chose to
# use a different compiler than CXX by specifying R_COMPILER. You may also chose
# a different directory than the default for R_HOME and R_LIBS, otherwise
# this will use the R command to find them

R_COMPILER ?=$(CXX)
R_LIB_DIR  ?=$(shell R --slave --vanilla -e "writeLines(.libPaths()[1])")
R_HOME_DIR ?=$(shell R --slave --vanilla -e "writeLines(paste(strsplit(Sys.getenv('R_HOME'),':')[[1]][1],sep=''))")

R_SOURCE = Rtranscpp/src/r_datatable.cpp Rtranscpp/src/r_defaults.cpp \
Rtranscpp/src/r_gene.cpp Rtranscpp/src/r_mode.cpp Rtranscpp/src/r_organism.cpp \
Rtranscpp/src/r_pwm.cpp Rtranscpp/src/r_tf.cpp

R_OBJECT=$(R_SOURCE:.cpp=.o)
R_HEADER=$(R_SOURCE:.cpp=.h)

$(R_OBJECT) : $(R_SOURCE) $(R_HEADER)
	$(R_COMPILER) $(RFLAGS) -I$(R_HOME_DIR)/include -Isrc/ -I$(R_LIB_DIR)/Rcpp/include $(@:.o=.cpp) -o $@
	
$(OBJECT:.o=.R.o) : $(SOURCE) $(HEADER)
	$(CXX) $(RFLAGS) -I$(BOOST_DIR) $(@:.R.o=.cpp) -o $@
	
src/utils.R.o: src/utils.cpp $(HEADER)
	$(R_COMPILER) $(RFLAGS) -I$(BOOST_DIR) -D R_LIB -I$(R_HOME_DIR)/include -I$(R_LIB_DIR)/Rcpp/include src/utils.cpp -o src/utils.R.o
	
Rtranscpp: $(R_OBJECT) $(OBJECT:.o=.R.o) src/utils.R.o
	ar cr Rtranscpp/src/liborganism.a $(OBJECT:.o=.R.o) src/utils.R.o
	$(R_COMPILER) -shared -O3 -o Rtranscpp/src/Rtranscpp.so $(R_OBJECT) -LRtranscpp/src -lorganism
	R CMD INSTALL Rtranscpp
	
	
#################################################################################

# compilation of the Matlab interface

#################################################################################

	
# Matlab is VERY picky about which compiler version it uses. I highly reccomend
# using the version that was used to compile matlab and using an identical
# version. Change this by setting the environmental variable MATLAB_COMPILER.
# Make sure that the matlab bin/ directory is in your PATH. 

MATLAB_COMPILER ?=$(CXX)
TMP_DIR         = $(shell which matlab)
MATLAB_INCLUDE  = $(TMP_DIR:bin/matlab=extern/include)

# compiles separate object file for matlab use
$(OBJECT:.o=.matlab.o) : $(SOURCE) $(HEADER)
	$(MATLAB_COMPILER) -c -fPIC -O3 -I$(BOOST_DIR) $(@:.matlab.o=.cpp) -o $@
	
src/utils.matlab.o: src/utils.cpp $(HEADER)
	$(R_COMPILER) -c -fPIC -O3 -I$(BOOST_DIR) -I$(MATLAB_INCLUDE) -D MEX src/utils.cpp -o src/utils.matlab.o

matlab: $(OBJECT:.o=.matlab.o) src/utils.matlab.o matlab/organism_interface_mex.cpp
	mex -g matlab/organism_interface_mex.cpp $(OBJECT:.o=.matlab.o) src/utils.matlab.o -o matlab/organism_interface_mex.mexa64
	
#################################################################################

# Cleanup

#################################################################################

# make clean removes everything
clean:
	rm -f src/*.o
	rm -f Rtranscpp/src/*.o
	rm -f Rtranscpp/src/*.so
	rm -f Rtranscpp/src/*.a
	rm -f matlab/*.mex64
	
# just cleanup R
cleanR:
	rm -f src/*.R.o
	rm -f Rtranscpp/src/*.o
	rm -f Rtranscpp/src/*.so
	rm -f Rtranscpp/src/*.a
	
# just cleanup matlab
cleanmatlab:
	rm -f src/*.matlab.o
	rm -f matlab/*.mex64
	

	


