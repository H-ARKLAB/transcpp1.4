XML_CFLAGS = `xml2-config --cflags`
XML_LIBS =  `xml2-config --libs`

CC  ?= gcc
CXX ?= g++

MATLAB_DIR ?=/usr/local/matlab-R2013b
RCPP_DIR   ?=/users/kenneth/RLibs/Rcpp/include
R_DIR      ?=/usr/include/R
BOOST_DIR  ?=/users/kenneth/Boost1
PARSA_ROOT ?=/users/kenneth/Annealer/neoParSA-1
PARSA_DIR = $(PARSA_ROOT)/parsa

LDLIBS += -L$(PARSA_ROOT)/build/lib -lparsa
LIBPARSA = $(PARSA_ROOT)/build/lib/libparsa.a

ifdef PARALLEL
	ifeq ($(CXX),icpc)
		PFLAGS = -openmp -DPARALLEL
	else
		PFLAGS = -fopenmp -DPARALLEL
	endif
endif

ifdef DEBUG
  FLAGS = -g -Wall -Wstrict-aliasing=0 -O2 -I$(BOOST_DIR) -I$(PARSA_DIR) $(XML_CFLAGS)
else
  ifeq ($(CXX),icpc)
		FLAGS = -O3 -I$(BOOST_DIR) -I$(PARSA_DIR) $(XML_CFLAGS) $(PFLAGS)
	else
		FLAGS = -O3 -I$(BOOST_DIR) -I$(PARSA_DIR) $(XML_CFLAGS) $(PFLAGS)
	endif
endif

ifdef MEX
	CXX += -fPIC
endif

ifdef R
  CXX += $(CXX) -O3 -fPIC
  FLAGS = -m64 -I/usr/include/R -DNDEBUG  -I/usr/local/include -I"/users/kenneth/RLibs/Rcpp/include"
endif

TRANSC     = competition.o sequence.o score.o coeffects.o cooperativity.o scalefactor.o fasta.o mode.o pwm.o TF.o gene.o nuclei.o datatable.o twobit.o parameter.o organism.o utils.o subgroup.o bindings.o bindingsite.o distance.o promoter.o quenching.o transcpp.o
PTRANSC    = competition.o sequence.o score.o coeffects.o cooperativity.o scalefactor.o fasta.o mode.o pwm.o TF.o gene.o nuclei.o datatable.o twobit.o parameter.o organism.o utils.o subgroup.o bindings.o bindingsite.o distance.o promoter.o quenching.o ptranscpp.o
SCRAMBLE   = competition.o sequence.o score.o coeffects.o cooperativity.o scalefactor.o fasta.o mode.o pwm.o TF.o gene.o nuclei.o datatable.o twobit.o parameter.o organism.o utils.o subgroup.o bindings.o bindingsite.o distance.o promoter.o quenching.o scramble.o
UNFOLD     = competition.o sequence.o score.o coeffects.o cooperativity.o scalefactor.o fasta.o mode.o pwm.o TF.o gene.o nuclei.o datatable.o twobit.o parameter.o organism.o utils.o subgroup.o bindings.o bindingsite.o distance.o promoter.o quenching.o unfold.o
PARSEFASTA = fasta.o parsefasta.o
CALCSCORE  = competition.o sequence.o fasta.o gene.o pwm.o TF.o coeffects.o cooperativity.o scalefactor.o parameter.o utils.o twobit.o distance.o
MATLAB     = competition.o sequence.o score.o coeffects.o cooperativity.o scalefactor.o fasta.o mode.o pwm.o TF.o gene.o nuclei.o datatable.o twobit.o parameter.o organism.o utils.o subgroup.o bindings.o bindingsite.o distance.o promoter.o quenching.o

test: $(MATLAB) test.o
	$(CXX) $(MATLAB) test.o -o tes $(LDLIBS) 
	
all: transcpp scramble unfold test_moves

mat: $(MATLAB)
ifdef MEXOPT_DIR
	$(MATLAB_DIR)/bin/mex -f $(MEXOPT_DIR) -g matlab/organism_interface_mex.cpp $(MATLAB) -o matlab/organism_interface_mex.mexa64
else
	$(MATLAB_DIR)/bin/mex -g matlab/organism_interface_mex.cpp $(MATLAB) -o matlab/organism_interface_mex.mexa64
endif

rlib: $(MATLAB)
	rm -rf Rtranscpp/src/*.
	rm -rf Rtranscpp/src/*.o
	rm -rf Rtranscpp/src/*.so
	ar cr Rtranscpp/src/liborganism.a $(MATLAB) 
	R CMD INSTALL Rtranscpp
	

transcpp: $(TRANSC) $(LIBPARSA)
	$(CXX) $(TRANSC) $(XML_LIBS) $(PFLAGS) -o transcpp $(LDLIBS) 
	rm transcpp.o
	
ptranscpp: $(PTRANSC) $(LIBPARSA)
	$(CXX) $(PTRANSC) $(XML_LIBS) $(PFLAGS) -m64 -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic -fPIC -Wl,-z,noexecstack -I/usr/include/mpich-x86_64 -L/usr/lib64/mpich/lib -lmpichcxx -Wl,-rpath -Wl,/usr/lib64/mpich/lib -lmpich -lopa -lmpl -lrt -lpthread -o ptranscpp $(LDLIBS) 
	rm ptranscpp.o
	
test_moves: $(MATLAB) test_moves.o $(LIBPARSA)
	$(CXX) $(MATLAB) test_moves.o $(XML_LIBS) $(PFLAGS) -o test_moves $(LDLIBS) 
	rm test_moves.o

scramble: $(SCRAMBLE) $(LIBPARSA)
	$(CXX) $(SCRAMBLE) $(XML_LIBS) $(PFLAGS) -o scramble $(LDLIBS) 
	rm scramble.o

unfold: $(UNFOLD) $(LIBPARSA)
	$(CXX) $(UNFOLD) $(XML_LIBS) $(PFLAGS) -o unfold $(LDLIBS) 
	rm unfold.o

getsequence: twobit.o utils.o getsequence.o
	$(CXX) -O3 twobit.o utils.o getsequence.o -o getseq
	
parsefasta: fasta.o parsefasta.o
	$(CXX) -O3 fasta.o parsefasta.o -o parsefasta
	
calcscores: $(CALCSCORE) calc_scores.o
	$(CXX) -O3 $(CALCSCORE) calc_scores.o -o calcscores

test_moves.o: test_moves.cpp
	$(CXX) -c $(FLAGS) test_moves.cpp
	
competition.o: competition.cpp
	$(CXX) -c $(FLAGS) competition.cpp
	
calc_scores.o: calc_scores.cpp
	$(CXX) -c $(FLAGS) calc_scores.cpp
	
fasta.o: fasta.cpp
	$(CXX) -c $(FLAGS) fasta.cpp
	
scramble.o: scramble.cpp
	$(CXX) -c $(FLAGS) scramble.cpp

sequence.o: sequence.cpp
	$(CXX) -c $(FLAGS) sequence.cpp
	
score.o: score.cpp
	$(CXX) -c $(FLAGS) score.cpp
	
cooperativity.o: cooperativity.cpp
	$(CXX) -c $(FLAGS) cooperativity.cpp
	
scalefactor.o: scalefactor.cpp
	$(CXX) -c $(FLAGS) scalefactor.cpp
	
parsefasta.o: parsefasta.cpp
	$(CXX) -c $(FLAGS) parsefasta.cpp
	
pwm.o: pwm.cpp
	$(CXX) -c $(FLAGS) pwm.cpp
	
nuclei.o: nuclei.cpp
	$(CXX) -c $(FLAGS) nuclei.cpp

TF.o: TF.cpp
	$(CXX) -c $(FLAGS) TF.cpp
	
gene.o: gene.cpp
	$(CXX) -c $(FLAGS) gene.cpp
	
datatable.o: datatable.cpp
	$(CXX) -c $(FLAGS) datatable.cpp

ifdef MEX
utils.o: utils.cpp
	$(CXX) -c -D MEX -I$(MATLAB_DIR)/extern/include $(FLAGS) utils.cpp
else 
ifdef R
utils.o: utils.cpp
	$(CXX) -c -D R_LIB -I$(RCPP_DIR) -I$(R_DIR) $(FLAGS) utils.cpp
else 
utils.o: utils.cpp
	$(CXX) -c $(FLAGS) utils.cpp
endif
endif
	
twobit.o: twobit.cpp
	$(CXX) -c $(FLAGS) twobit.cpp
	
organism.o: organism.cpp
	$(CXX) -c $(FLAGS) organism.cpp

parameter.o: parameter.cpp
	$(CXX) -c $(FLAGS) parameter.cpp
	
transcpp.o: transcpp.cpp
	$(CXX) -c $(FLAGS) transcpp.cpp
	
ptranscpp.o: ptranscpp.cpp
	$(CXX) -c $(FLAGS) -DUSE_BOOST -m64 -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic -fPIC -Wl,-z,noexecstack -I/usr/include/mpich-x86_64 -L/usr/lib64/mpich/lib -lmpichcxx -Wl,-rpath -Wl,/usr/lib64/mpich/lib -lmpich -lopa -lmpl -lrt -lpthread ptranscpp.cpp
	
getsequence.o: getsequence.cpp
	$(CXX) -c $(FLAGS) getsequence.cpp
	
subgroup.o: subgroup.cpp
	$(CXX) -c $(FLAGS) subgroup.cpp
	
bindingsite.o: bindingsite.cpp
	$(CXX) -c $(FLAGS) bindingsite.cpp
	
bindings.o: bindings.cpp
	$(CXX) -c $(FLAGS) bindings.cpp
	
quenching.o: quenching.cpp
	$(CXX) -c $(FLAGS) quenching.cpp
	
distance.o: distance.cpp
	$(CXX) -c $(FLAGS) distance.cpp
	
promoter.o: promoter.cpp
	$(CXX) -c $(FLAGS) promoter.cpp
	
mode.o: mode.cpp
	$(CXX) -c $(FLAGS) mode.cpp
	
unfold.o: unfold.cpp
	$(CXX) -c $(FLAGS) unfold.cpp
	
coeffects.o: coeffects.cpp
	$(CXX) -c $(FLAGS) coeffects.cpp
	
clean:
	rm -rf *.o 
