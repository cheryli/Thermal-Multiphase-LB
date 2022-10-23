.SUFFIXES:
# FC = mpiifort
FC = mpif90 
OBJDIR=objs
SRCDIR=src

OPT=-acc -fast -Minline=collision,collisionT

SRC=commondata.F90 mpi_starts.F90 main.F90 initial.F90 flow_update.F90  \
	message_exc.F90 \
	# output.F90
objects=$(addprefix $(OBJDIR)/, $(subst .F90,.o,$(SRC)))

macros= \
	# -D__OUT_PUT__

#define steadyFlow  
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
#define __OUT_PUT__


default: mpi_test

mpi_test : dirs $(objects)
	$(FC)  $(macros) $(OPT) -o $@ $(objects)

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(FC) $< $(macros) $(OPT) -c -o $@


.PHONY : clean  clean_all dirs
clean :
	rm -f *.mod *.plt *.dat mpi_test nohup $(objects)

# clean_all: 
# 	rm -f *.plt *.dat mpi_test  $(objects) 

dirs:
	mkdir -p $(OBJDIR)/