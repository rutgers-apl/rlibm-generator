objects = 	Log2
		

CC=g++

INCLUDES= -I$(SOFTPOSITPATH)/source/include -I$(SOPLEXPATH)/src
CPP_FLAGS= -O2 -fvisibility=hidden -fvisibility-inlines-hidden --std=c++17 -I../../../helper/
LD_FLAGS= $(SOFTPOSITPATH)/build/Linux-x86_64-GCC/softposit.a $(SOPLEXPATH)/build/lib/libsoplex.a -lmpfr -lgmp  -lz -lm

all: $(objects)
	
%.o: %.cpp
	$(CC) $(CPP_FLAGS) -c $< -o $@ $(INCLUDES)
	
$(objects): %: ../../../helper/PolyFinder.o ../../../helper/IntervalRange.o %.o 
	$(CC) $^ $(LD_FLAGS) -o $@

clean:
	rm -rf *.o ../../../helper/*.o $(objects) *~ *.dSYM dump.lp
