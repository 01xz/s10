CC = gcc
CFLAGS = -O3 -Wall -I/home/liam/repos/SuiteSparse-5.10.1/include
LDFLAGS = -L/home/liam/repos/SuiteSparse-5.10.1/lib -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis -llapack -lblas -lklu -lm

.PHONY: all
all: kludemo

kludemo: kludemo.c
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm kludemo

