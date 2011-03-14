FC         = gfortran
 
FFLAGS	   = -ggdb -dH -c -w -O3

LDFLAGS    = 
 
LINKER	   = $(FC) 
 
LIB        = -L/opt/Cern/2003/lib
LIBB       = -L/usr/lib
LIBS       = -lpacklib -lmathlib -lkernlib  
LIBS2      = -lpdflib -lpacklib -lmathlib  
LIBS3      = `cernlib packlib mathlib`
LIBP        = -L/opt/pythia/lib
LIBPS       = -lpythia
 
.f.o:	 
	$(FC) $(FFLAGS) $*.f
 

files	= lepton_tev5.o setlhc.o opensub.o sigandnev.o lepcuts.o decaychain2.o ntuplehisto-dummy.o position_dist.o subroutinesv4.o pythia-6.4.20.o
all:	$(files)  
	$(LINKER) $(LDFLAGS) -o lepton_tev  $(files)
files10	= branch_tau_new.o setlhc.o opensub.o sigandnev2.o lepcuts.o decaychain2.o ntuplehisto-dummy.o position_dist.o subroutinesv4.o pythia-6.4.20.o smear.o
etau:	$(files10)  
	$(LINKER) $(LDFLAGS) -o branch_etau  $(files10)
files2	= ensayo.o setlhc.o opensub.o sigandnev.o lepcuts.o decaychain.o ntuplehisto-dummy.o position_dist.o subroutines.o
ensayo:	$(files2)  
	$(LINKER) $(LDFLAGS) -o ensayo  $(files2) $(LIBP) $(LIBPS)
files3	= decay.o setlhc.o opensub.o lepcuts.o ntuplehisto-dummy.o decaychain.o
decay:	$(files3)  
	$(LINKER) $(LDFLAGS) -o decay  $(files3) $(LIBP) $(LIBPS)
files4	= lepton_tev4.o setlhc.o opensub.o sigandnev.o lepcuts.o decaychain2.o position_dist.o subroutinesv4.o pythia-6.4.20.o inithisto.o savehisto.o
histo:	$(files4)  
	$(LINKER) $(LDFLAGS) -o lepton_histo  $(files4) $(LIBB) $(LIBS3)
files5	= decaycs.o setlhc.o opensub.o lepcuts.o ntuplehisto-dummy.o decaychain.o pythia-6.4.13.o
decayc:	$(files5)  
	$(LINKER) $(LDFLAGS) -o decaycs  $(files5) 
files6 = higgs.o setlhc.o opensub.o lepcuts.o ntuplehisto-dummy.o	\
decaychain.o pythia-6.4.13.o
higgs:	$(files6)  
	$(LINKER) $(LDFLAGS) -o higgs  $(files6) 
files7	= decaycsfull.o setlhc.o opensub.o sigandnev.o lepcuts.o decaychain.o subroutines.o pythia-6.4.13.o
decayf:	$(files7)  
	$(LINKER) $(LDFLAGS) -o decaycsfull  $(files7) 
files8	= decaycshiggs.o setlhc.o opensub.o sigandnev.o lepcuts.o decaychain.o subroutines.o pythia-6.4.14.o
decayh:	$(files8)  
	$(LINKER) $(LDFLAGS) -o decaycshiggs  $(files8) 
files9	= dat-rz.o
dat-rz:	$(files9)  
	$(LINKER) $(LDFLAGS) -o dat-rz  $(files9) $(LIB) $(LIBS)

