# implicit pattern matching rules for generating .o files from .f and .c files

#  %.o: %.f
#	$(COMPILE.f) $< ; mv $(notdir $@) $@ 
#  %.o: %.c
#	$(COMPILE.c) $< ; mv $(notdir $@) $@ 

  %.o: %.f
	gfortran -c  $< 
  %.o: %.c
	$(COMPILE.c) $< 


FIT5DMOLS2 = fit5d-new-mols2.o potenl-AVR-h2oh2.o cubspl-claire.o t.o

	
pesh2oh2:  $(FIT5DMOLS2)
	gfortran $(FIT5DMOLS2) -o $@