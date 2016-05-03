CC = gcc
CXX = g++
LDFLAGS =  -lmpfr -lgmp -Wl,--rpath,/usr/local/lib64

INCLIST = 

OBJS = alg_remez.o poly.o

CFLAGS=  -O3
CXXFLAGS=  -O3
DFLAGS +=  -DGMP

# Targets

poly:	poly.o alg_remez.o
	${CXX} alg_remez.o poly.o ${LDFLAGS} -o poly

.SUFFIXES:
.SUFFIXES: .o .C .c

$(BIN):  $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(BIN)

${OBJS}: Makefile

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<

clean:
	rm -f *.dat *.o  *~ *\#* *.stackdump $(BIN)
