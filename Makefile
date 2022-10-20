SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic


#------------------------------------------------------------
#
# When you adapt this makefile to compile LayersC
# please copy this makefile and set HVC, EAF, GHSS to
# the directories where these codes are installed.
#
#------------------------------------------------------------

HVC	= ../HVC
EAF	= ../EAF3D
GHSS	= ../gHSS-master
DEP 	= $(HVC) $(EAF) $(GHSS)

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

# CPLEXDIR      = ../../..
# CONCERTDIR    = ../../../../concert

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio201/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio201/concert


# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CC  = gcc

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

COPT  = -m64 -fPIC -O3

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# For dynamic linking
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIB      = cplex$(dynamic:yes=2010)
run           = $(dynamic:yes=LD_LIBRARY_PATH=$(CPLEXBINDIR))

CLNDIRS   = -L$(CPLEXLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CLNFLAGS  = -l$(CPLEXLIB) -lm -lpthread -ldl


all:
	make -C $(HVC)
	make -C $(EAF)
	make -C $(GHSS)
	make all_c

	
	
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include


CFLAGS  = $(COPT) -I$(CPLEXINCDIR)

#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------

C_EX = layers


all_c: $(C_EX)


execute_c: $(C_EX)
	
	$(run) ./$(C_EX)

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp


# ------------------------------------------------------------
#

layers: layers.o io.o timer.o main-layers.c
	$(CC) $(CFLAGS) $(CLNDIRS) -o $(C_EX) main-layers.c io.o timer.o layers.o $(HVC)/hvc.o  $(HVC)/hv-plus.o $(EAF)/eaf3d.o  $(HVC)/avl.o $(GHSS)/gHSS.o $(CLNFLAGS) $(VFLAGS)


%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

example: $(C_EX)
	$(run) ./$(C_EX) test.inp -k 2 -r "4 4 4"
	

