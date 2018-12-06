# CCFDV3.0
CCFDV3.0 is an open-source written in fortran for the analysis of partial differential equations(PDEs) problems on structured meshes (CGNS format).

CCFDV3.0
in build/CCFD_MPI/Makefile or build/CCFD_SEQ/Makefile:

OBJS_DIR = obj/
EXE_DIR  = bin/

EXE      = CCFD_MPI
FC       = mpif90
IDIR     = -I/home/ccfd/cgns/include -I/home/ccfd/metis/include 
CFLAGS   = -r8 -Wall -O2 -w  -DDEBUG 
LFLAGS   = 
LIBS     = -L/home/ccfd/cgns/lib -lcgns -L/home/ccfd/metis/lib -lmetis
        
VPATH    = $(SRC_DIR_F90d1):$(OBJS_DIR):$(SRC_DIR_F90d2):$(OBJS_DIR):$(SRC_DIR_F90d3):$(OBJS_DIR):$(SRC_DIR_F90d4):$(OBJS_DIR):$(SRC_DIR_F90d5):$(OBJS_DIR):$(SRC_DIR_F90d6):$(OBJS_DIR):$(SRC_DIR_F90d7):$(OBJS_DIR):$(SRC_DIR_F90d8):$(OBJS_DIR):$(SRC_DIR_F90d9):$(OBJS_DIR):$(SRC_DIR_F90d10):$(OBJS_DIR):$(SRC_DIR_F90d11):$(OBJS_DIR)
OBJS     = $(addprefix $(OBJS_DIR), $(OBJS_F90d1) $(OBJS_F90d2) $(OBJS_F90d3) $(OBJS_F90d4) $(OBJS_F90d5) $(OBJS_F90d6) $(OBJS_F90d7) $(OBJS_F90d8) $(OBJS_F90d9) $(OBJS_F90d10) $(OBJS_F90d11))
    @ echo "========================CCFDV3.0 complier done ========================================="
    @ echo "================You can see the application at CCFDV3.0/build/bin/CCFDV3_SEQ ==========="
    @ echo "==========YOU CAN COPY THE APPLICATION TO THE COMPUTING FILE OR LINE TO THE FILES  ====="
all : $(EXE)

#we need cgns and Metis lib and -lcgns -lmetis
#so you can edit the path at Makefile 

#make the ccfdv3.0 install

make all 
#or
make
#make the ccfdv3.0 uninstall
make clean

#you can delete the -DPMPI when you want complier the ccfd_seq in a host processor computing!
