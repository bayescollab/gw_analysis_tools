####################################################
EXE:=bin/test_fishers.exe bin/test_SNR.exe bin/test_util.exe bin/test_waveform.exe bin/test_waveform_util.exe bin/test_cosmology.exe bin/test_response_functions.exe
# Not compiling bin/test_mcmc.exe program because deprecated. Now contained in Bayesship
#EXE:=bin/test_waveform.exe 
####################################################
# Use the default CXX compiler and linker
#
CCLOCAL:=$(CXX)
LINKERLOCAL:=$(CXX)
#CFLAGS:=-fopenmp -O3  -Wall -I../include -I/usr/include/hdf5/serial
#CFLAGS:=-fopenmp  -Wall -I../include -I/home/sperkins/git-repos/PTRJMCMC/build/install/include -g 
#LFLAGS:= -fopenmp -L../lib  -L/home/sperkins/git-repos/gw_analysis_tools/build/lib -L/home/sperkins/git-repos/PTRJMCMC/build/install/lib
CFLAGS:=-fopenmp  -Wall -I../include  -g 
LFLAGS:= -fopenmp -L../lib  -L/opt/gw_analysis_tools/build/lib
#LIBS:=-lm -lgsl -lgwat -llal -llalsimulation
#LIBS:=-lm -lgsl -lgwat -llal -llalsimulation -lhdf5 -lhdf5_cpp
#LIBS:=-lm -lgsl -lgwat -llal -llalsimulation 
#LIBS:=-lm -lgsl -lgwat -llal -llalsimulation
#LIBS:=-lm -lgsl -lgwat -llal -llalsimulation 
#LIBS:=-lm -lgsl -lgwat -lptrjmcmc
LIBS:=-lm -lgsl -lgwat -llal -llalsimulation 

SRCDIR:=src
ODIR:=build
DATA_DIR:=data
PLOTS_DIR:=plots
EXEDIR:=bin
LIBDIR:=lib
IDIR:=include

SRCEXT:=cpp
IEXT:=h
EXEEXT:=exe

SRCTOT:=$(wildcard $(SRCDIR)/*.$(SRCEXT))
DEPS:=$(patsubst $(SRCDIR)/%.$(SRCEXT), $(IDIR)/%.$(IEXT), $(SRC))
OBJ:=$(patsubst $(SRCDIR)/%.$(SRCEXT), $(ODIR)/%.o,$(SRCTOT))

.PHONY: all

all: $(EXE) $(DATA_DIR) $(PLOTS_DIR)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS)| $(ODIR)
	@echo COMPILING
	$(CCLOCAL) $(CFLAGS) -I$(IDIR)/ -c -o $@ $<

$(ODIR):
	@echo "Making build directory"
	mkdir $(ODIR)
$(DATA_DIR):
	@echo "Making data directory"
	mkdir $(DATA_DIR)
$(PLOTS_DIR):
	@echo "Making plots directory"
	mkdir $(PLOTS_DIR)
$(EXEDIR):
	@echo "Making bin directory"
	mkdir ./bin

.PHONY: clean

clean: 
	rm $(ODIR)/*.o 
	rmdir $(ODIR)
.PHONY: remove
remove: 
	rm $(ODIR)/*.o 
	rmdir $(ODIR)
	rm $(EXE) 
	rmdir $(EXEDIR) 
#########################################
# Executable number dependent parameters
#########################################
EXE0:=bin/test_fishers.exe
EXE1:=bin/test_SNR.exe
EXE2:=bin/test_util.exe
EXE3:=bin/test_waveform.exe
EXE4:=bin/test_waveform_util.exe
EXE5:=bin/test_mcmc.exe
EXE6:=bin/test_cosmology.exe
EXE7:=bin/test_response_functions.exe
OBJEXE0:=build/test_fishers.o
OBJEXE1:=build/test_SNR.o
OBJEXE2:=build/test_util.o
OBJEXE3:=build/test_waveform.o
OBJEXE4:=build/test_waveform_util.o
OBJEXE5:=build/test_mcmc.o
OBJEXE6:=build/test_cosmology.o
OBJEXE7:=build/test_response_functions.o
OBJDEP:=$(OBJ)
OBJDEP:=$(filter-out $(OBJEXE0), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE1), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE2), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE3), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE4), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE5), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE6), $(OBJDEP))
OBJDEP:=$(filter-out $(OBJEXE7), $(OBJDEP))
$(EXE0): $(OBJDEP) $(OBJEXE0) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE1): $(OBJDEP) $(OBJEXE1) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE2): $(OBJDEP) $(OBJEXE2) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE3): $(OBJDEP) $(OBJEXE3) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE4): $(OBJDEP) $(OBJEXE4) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE5): $(OBJDEP) $(OBJEXE5) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE6): $(OBJDEP) $(OBJEXE6) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
$(EXE7): $(OBJDEP) $(OBJEXE7) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $< $(LIBS)
