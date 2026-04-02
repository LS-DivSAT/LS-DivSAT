EXEC      = minisat
DEPDIR    = mtl utils

MROOT= ../minisat

patch-minisat: minisat.patch
	@cd ../minisat && git restore . && git apply ../minisat_ext/minisat.patch && cd ..

include $(MROOT)/mtl/template.mk

extminisat: patch-minisat Solver.o Ext.o BlackBoxSolver.o

.PHONY: extminisat patch-minisat