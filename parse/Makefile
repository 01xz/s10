CC = gcc
CFLAGS = -std=c99 -O3 -Wall

TARGET = main
SRC = parser_demo.c parser.c
OBJ = parser_demo.o parser.o

$(TARGET): $(SRC)
	$(CC) -o $@ $< $(CFLAGS)

clean:
	rm $(TARGET)
	rm *.o

