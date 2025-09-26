####################################################
EXE:=bin/injection.exe
####################################################
# Use the default CXX compiler and linker
CCLOCAL:=$(CXX)
LINKERLOCAL:=$(CXX)
CFLAGS:=-fopenmp -O2  -Wall -g
LFLAGS:= -fopenmp
LIBS:=-lm -lgsl -lgwat -lbayesship

SRCDIR:=src
ODIR:=build
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

all: $(EXE)

$(ODIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(DEPS) | $(ODIR)
	@echo COMPILING
	$(CCLOCAL) $(CFLAGS) -I$(IDIR)/ -c -o $@ $<

$(ODIR):
	@echo "Making build directory"
	mkdir ./build
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
EXE0:=bin/injection.exe
OBJEXE0:=build/injection.o
OBJDEP:=$(OBJ)
OBJDEP:=$(filter-out $(OBJEXE0), $(OBJDEP))
$(EXE0): $(OBJEXE0) $(OBJDEP) | $(EXEDIR)
	@echo LINKING
	$(LINKERLOCAL) $(LFLAGS) -o $@ $^ $(LIBS)
