objects = 	Exp \
		Exp2 \
 		Exp10 \
		Log \
		Log10 \
		Log2 \
		Sinpi \
		Cospi \
		Sqrt \
		Cbrt

CC=g++

INCLUDES= -I$(SOFTPOSITPATH)/source/include -I$(SOPLEXPATH)/src
CPP_FLAGS= -fvisibility=hidden -fvisibility-inlines-hidden --std=c++17 -I../../../helper/
LD_FLAGS= $(SOFTPOSITPATH)/build/Linux-x86_64-GCC/softposit.a $(SOPLEXPATH)/build/lib/libsoplex.a -lmpfr -lgmp  -lz -lm

all: $(objects)
	
%.o: %.cpp
	$(CC) $(CPP_FLAGS) -c $< -o $@ $(INCLUDES)
	
$(objects): %: %.o ../../../helper/PolyFinder.o ../../../helper/IntervalRange.o
	$(CC) $^ $(LD_FLAGS) -o $@

clean:
	rm -rf *.o ../../../helper/*.o $(objects) *.lp *~ *.dSYM
