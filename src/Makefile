include ../Makefile.in

program = test_unsteady_transport

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

#=====================================================================#
# Targets
#=====================================================================#

all: $(program)

lib: $(objects) 
	ar rvs $@ $(objects) lib${TARGET}.a
	@echo " ----------- ${TARGET} created ----------- "

$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)

distclean: clean

clean:
	@rm -f *.o *.mod $(program)

neat:
	@rm -f *.o *.mod

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean neat distclean 

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
