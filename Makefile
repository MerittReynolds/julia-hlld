FFLAGS=-Wall

run: hlld_test
	./hlld_test

hlld_test: hlld_test.f90 hlld.o
	gfortran $(FFLAGS) $^ -o $@

%.o : %.f90
	gfortran $(FFLAGS) -c $< 

clean:
	rm -f *.o *.mod hlld_test
