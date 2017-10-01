DBG = -Wall
FC = ftn
CC = cc
FLAGS = -O3
LIBS = -lm -lgfortran
SRCDIR = src
BINDIR = bin
OBJDIR = obj

all:
	make mpio
	make spline
	make xapiir
	make utils
	make params
	make stf
	make krg 

mpio:
	$(CC) -c $(DBG) $(FLAGS) $(SRCDIR)/sord_mpio.c -o $(OBJDIR)/$@.o

spline:
	$(CC) -c $(DBG) $(FLAGS) $(SRCDIR)/spline.c -o $(OBJDIR)/$@.o

xapiir:
	$(FC) -c $(DBG) $(FLAGS) $(SRCDIR)/xapiir.f -o $(OBJDIR)/$@.o

utils:
	$(CC) -c $(DBG) $(FLAGS) $(SRCDIR)/utils.c -o $(OBJDIR)/$@.o

params:
	$(CC) -c $(DBG) $(FLAGS) $(SRCDIR)/params.c -o $(OBJDIR)/$@.o
    
stf:
	$(CC) -c $(DBG) $(FLAGS) $(SRCDIR)/stf.c -o $(OBJDIR)/$@.o
	
krg:
	$(CC) $(DBG) $(FLAGS) -o $(BINDIR)/krg $(SRCDIR)/krg.c $(OBJDIR)/*.o $(LIBS)

clean:
	-@rm $(OBJDIR)/*.o $(BINDIR)/krg 2>/dev/null || true

install: 
	cp $(BINDIR)/krg ~/executables
