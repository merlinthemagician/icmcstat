# Makefile
#
CC      = gcc#
CFLAGS	=-Wall -pedantic#

rm	= rm -f#				delete file

all:	_
clean::		;@ $(MAKE) T='$T' _clean
_clean:	_	;  $(rm) *.o $T a.out core *.tmp *.ps *.bak
run::	_
_:		;@ echo -------------------- $D --------------------


D =	Statistical analysis of mode changes

T = mp_parameter.o mp_mcmc.o likelihood.o bernoulliDist.o mp_proposal.o icmcstat$x cp_proposal.o cp_likelihood.o cp_prior.o generateTestCurrents$x icmcstatGEO$x icmcstatNEGBINK1$x icmcstatNEGBINK2$x icmcstatNJgeoNEGBINK1$x icmcstatNJgeoNEGBINK2$x icmcstatFixed$x icmcstatBin$x dat2SNPindex$x gcmcstat$x

all:	$T

run::	icmcstat$x; 	./icmcstat testdata/IPR2_10uMIP3_5mMATP_10nMCa_C4_T02.dat 10000 42 results/IPR2_10uMIP3_5mMATP_10nMCa_C4_T02_s42

INCLUDE=.
GSLINC = `gsl-config --cflags`
INC2=/opt/local/include/

GSL = `gsl-config --cflags` `gsl-config --libs` -lm

MCMC= mp_parameter.o mp_mcmc.o likelihood.o mp_proposal.o
BER = cp_proposal.o bernoulliDist.o cp_likelihood.o cp_prior.o

mp_parameter.o:	mp_parameter.c;
		$(CC) $(CFLAGS) -c -o $@ mp_parameter.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

mp_proposal.o:	mp_proposal.c;
		$(CC) $(CFLAGS) -c -o $@ mp_proposal.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

mp_mcmc.o:	mp_mcmc.c;
		$(CC) $(CFLAGS) -c -o $@ mp_mcmc.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

likelihood.o:	likelihood.c
		$(CC) $(CFLAGS) -c -o $@ likelihood.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

bernoulliDist.o:	bernoulliDist.c
			$(CC) $(CFLAGS) -c -o $@ bernoulliDist.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

cp_proposal.o:	cp_proposal.c
		$(CC) $(CFLAGS) -c -o $@ cp_proposal.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

icmcstat$x:	$(MCMC) $(BER) icmcstat.c
		$(CC) $(CFLAGS) -DDEFAULT -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatBin$x:	$(MCMC) $(BER) icmcstat.c
		$(CC) $(CFLAGS) -DDEFAULT -DBIN -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatGEO$x:	$(MCMC) $(BER) icmcstat.c
		$(CC) $(CFLAGS) -DGEO -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatNEGBINK1$x:	$(MCMC) $(BER) icmcstat.c
			$(CC) $(CFLAGS) -DNEGATIVEBINK1 -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatNEGBINK2$x:	$(MCMC) $(BER) icmcstat.c
			$(CC) $(CFLAGS) -DNEGATIVEBINK2 -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

cp_likelihood.o:	cp_likelihood.c
			$(CC) $(CFLAGS) -c -o $@ cp_likelihood.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

cp_prior.o:	cp_prior.c
		$(CC) $(CFLAGS) -c -o $@ cp_prior.c -I$(INCLUDE) $(GSLINC) -I$(INC2)

generateTestCurrents$x:	generateTestCurrents.c
			$(CC) $(CFLAGS) -o $@ generateTestCurrents.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatNJgeoNEGBINK1$x:	$(MCMC) $(BER) icmcstat.c
				$(CC) $(CFLAGS) -DNJGEONEGATIVEBINK1 -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatNJgeoNEGBINK2$x:	$(MCMC) $(BER) icmcstat.c
				$(CC) $(CFLAGS) -DNJGEONEGATIVEBINK2 -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

icmcstatFixed$x:	$(MCMC) $(BER) icmcstat.c
		$(CC) $(CFLAGS) -DFIXED -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)

gcmcstat$x:	$(MCMC) $(BER) icmcstat.c
		$(CC) $(CFLAGS) -DBERNOULLIDATA -DDEFAULT -o $@ $(MCMC) $(BER) icmcstat.c -I$(INCLUDE) $(GSLINC) $(GSL)
