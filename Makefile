LIBRARY_DIR=$(ITENSOR)
include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk
INC+=-I${ITENSOR}
INC+=-I$(PWD)/include


LIBSPATH=-L$(ITENSOR)/lib
LIBSPATH+=$(LIBSLINK)



LIBS=-litensor 
LIBSG=-litensor-g 


#########################

CCFLAGS+=-I. $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -Wno-unused-variable -std=c++17 -O2 -std=gnu++1z
CCGFLAGS+=-I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS) 




CPPFLAGS_EXTRA += -O2 -std=c++17



DB=-g
CXX=g++
ND=-DNDEBUG

##################################################



# spinless electron_phonon model

spinless_electron_phonon_model: src/spinless_electron_phonon_model.cpp $(ITENSOR_LIBS) 
	$(CCCOM) $< -o bin/$@ $(CCFLAGS) $(INC) $(LIBFLAGS) 

# spinfulls electron_phonon model

#spinless_electron_phonon_model: src/spinless_electron_phonon_model.cpp $(ITENSOR_LIBS) 
#	$(CCCOM) $< -o bin/$@ $(CCFLAGS) $(INC) $(LIBFLAGS) 


clean:
	rm bin/*

