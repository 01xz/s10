CC = gcc
CFLAGS = -std=c99 -O3 -Wall -I/home/liam/repos/SuiteSparse-5.10.1/include
LDFLAGS = -L/home/liam/repos/SuiteSparse-5.10.1/lib -lcxsparse -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis -llapack -lblas

TARGET = fieldsolver2d
SRC = main.c mesh.c solver.c
OBJ = main.o mesh.o solver.o

$(TARGET): $(OBJ)
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $<

clean:
	rm $(TARGET)
	rm *.o

