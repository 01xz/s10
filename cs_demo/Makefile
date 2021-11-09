CC = gcc
CFLAGS = -O3 -Wall -I/home/liam/repos/SuiteSparse-5.10.1/include
LDFLAGS = -L/home/liam/repos/SuiteSparse-5.10.1/lib -lcxsparse -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis -llapack -lblas

demo: cs_demo2.c cs_demo.c
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm demo

