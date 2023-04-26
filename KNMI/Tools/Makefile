FFLAGS  = -I /usr/include
FFCOMP  = gfortran
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lnetcdff

thwt-EOBS.o:	thwt-EOBS.f
				$(FFCOMP) -c $(FFLAGS) -I $(NETCDFHOME)/include -o $@ $<
thwt-EOBS:		thwt-EOBS.o thwt-EOBS.f
				$(FFCOMP) -o $@ $< $(LDFLAGS)
thwt-HYRAS:			thwt-HYRAS.o thwt-HYRAS.f routines.f
				$(FFCOMP) -o $@ $< $(LDFLAGS)
scanlist:			scanlist.f routines.f
				$(FFCOMP) -o $@ $@.f90
thwt:				thwt.o pinpok.o thwt.f routines.f
				$(FFCOMP) -o $@ $< pinpok.o $(LDFLAGS)
petf:				petf.o petf.f routines.f
				$(FFCOMP) -o $@ $< $(LDFLAGS)
ncGRaDE:			ncGRaDE.o ncGRaDE.f routines.f
				$(FFCOMP) -o $@ $< $(LDFLAGS)
ncGRaDE0:			ncGRaDE0.o ncGRaDE0.f routines.f
				$(FFCOMP) -o $@ $< $(LDFLAGS)
rr-trans2.1_GRADE:		rr-trans2.1_GRADE.f routines.f
				$(FFCOMP) -o $@ $@.f
tm-trans2.1_GRADE:		tm-trans2.1_GRADE.f routines.f
				$(FFCOMP) -o $@ $@.f
seco02:			        seco02.f routines.f
				$(FFCOMP) -o $@ $@.f
seco03:			        seco03.f routines.f
				$(FFCOMP) -o $@ $@.f
xtreme_all:			xtreme_all.f routines.f
				$(FFCOMP) -o $@ $@.f
xtreme_sum:			xtreme_sum.f routines.f
				$(FFCOMP) -o $@ $@.f
xtreme_win:			xtreme_win.f routines.f
				$(FFCOMP) -o $@ $@.f
gumgev:			gumgev.f routines.f
				$(FFCOMP) -o $@ $@.f
areaavg:			areaavg.f routines.f
				$(FFCOMP) -o $@ $@.f
joinasciiseries:	joinasciiseries.f
				$(FFCOMP) -o $@ $@.f
	

pinpok.o:	pinpok.f90
	$(FFCOMP) -c $(FFLAGS) -o pinpok.o pinpok.f90
%.o:	%.f routines.f
	$(FFCOMP) -c $(FFLAGS) -o $@ $<

clean:
	rm -rf *.o
	rm -f thwt-EOBS_v15b
