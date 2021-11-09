CC = gcc
CFLAGS = -O3 -Wall -I/home/liam/repos/SuiteSparse-5.10.1/include
LDFLAGS = -L/home/liam/repos/SuiteSparse-5.10.1/lib -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis -llapack -lblas -lklu -lm

kludemo: kludemo.c
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS)

clean:
	rm kludemo

