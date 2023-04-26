FFLAGS  = -I /usr/include
FFCOMP  = gfortran
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lnetcdff

thwt-EOBS_v15.0_Rhine3.o:	thwt-EOBS_v15.0_Rhine3.f
				$(FFCOMP) -c $(FFLAGS) -o $@ $<

resample:			resample.f 
				$(FFCOMP) -o $@ $@.f
fvectors:			fvectors.f 
				$(FFCOMP) -o $@ $@.f
fvectors_select:			fvectors_select.f 
							$(FFCOMP) -o $@ $@.f
backtran:			backtran.f 
				$(FFCOMP) -o $@ $@.f

%.o:	%.f routines.f
	$(FFCOMP) -c $(FFLAGS) -o $@ $<

clean:
	rm -rf *.o
	rm -f thwt-EOBS_v15b
