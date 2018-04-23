# Operating system
OS = WINDOWS

# Compiler
CC = g++

# Compiler standard
#		Windows: issues with c++11 on 32 bit g++, used gnu++11 instead
#		Linux: g++ compiles with c++11
# You can play around with this if problems occur (i.e. gnu+11 not supported 
# on your windows compiler etc.)
ifeq ($(OS),WINDOWS)
	CSTD = c++11
else
	CSTD = c++11
endif

# Optimization flag
OFLAG = -O3

# Openmp flag
OMPFLAG = -fopenmp

# Static option flag
STATICFLAG = -static -static-libgcc -static-libstdc++

# Project file directories
INCLDDIR = include
SRCDIR = src
BUILDDIR = build
BINDIR = bin
TESTDIR = test

# include folders
EIGEN = $(INCLDDIR)/eigen__3_3_4
BOOST = $(INCLDDIR)/boost__1_66_0
AMGCL = $(INCLDDIR)/amgcl__801a2fb

# Compiler and flags
# -DEIGEN_DONT_ALIGN used for Windows 32 bit exe.
CFLAG = -w -std=$(CSTD) $(OFLAG) $(OMPFLAG) $(STATICFLAG) -I$(EIGEN) -I$(BOOST) -I$(AMGCL)

# Target and dependencies
ifeq ($(OS),WINDOWS)
	LINK_TARGET = $(BINDIR)/BES.exe
else
	LINK_TARGET = $(BINDIR)/BES
endif
OBJS = $(BUILDDIR)/BES.o \
	$(BUILDDIR)/Options.o \
  $(BUILDDIR)/FileReader.o \
	$(BUILDDIR)/Config.o \
	$(BUILDDIR)/MeshCommon.o \
	$(BUILDDIR)/Vertex.o \
	$(BUILDDIR)/Face.o \
	$(BUILDDIR)/Element.o \
	$(BUILDDIR)/Mesh.o \
	$(BUILDDIR)/MeshGmsh.o \
  $(BUILDDIR)/MatrixSparse.o \
  $(BUILDDIR)/CoefficientMatrix.o \
	$(BUILDDIR)/SolutionSnapshot.o \
	$(BUILDDIR)/Solution.o \
	$(BUILDDIR)/Gradient.o \
	$(BUILDDIR)/GradientGreenGauss.o \
	$(BUILDDIR)/Limiter.o \
	$(BUILDDIR)/LimiterVenkatakrishnan.o \
  $(BUILDDIR)/FluxJacobian.o \
	$(BUILDDIR)/TimeStepCalculator.o \
	$(BUILDDIR)/TimeStepCalculatorLocal.o \
  $(BUILDDIR)/TimeStepCalculatorGlobal.o \
	$(BUILDDIR)/Solver.o \
	$(BUILDDIR)/SolverExplicit.o \
  $(BUILDDIR)/SolverImplicit.o

# Make commands
all: $(LINK_TARGET)

$(LINK_TARGET): $(OBJS)
	$(CC) $(CFLAG) $(OBJS) -o $(LINK_TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAG) -c $^ -o $@

ifeq ($(OS),WINDOWS)
clean:
	del /q $(BUILDDIR)\*.o $(BINDIR)\BES.exe $(BINDIR)\BES.out
else
clean:
	rm -rf $(BUILDDIR)/*.o
	rm -rf $(BINDIR)/BES.exe
	rm -rf $(BINDIR)/BES
endif 
