include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CXXFLAG_DEF =

CXXFLAGS += -std=c++20 -Ofast ${CXXFLAG_DEF} -pedantic -march=native

std: std.o
	-${CLINKER} -o std std.o  ${PETSC_LIB}
	${RM} std.o

clean::
	@rm -f *~ std vtest *tmp
	@rm data/*
