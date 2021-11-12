CC = gcc
CFLAGS = -std=c99 -O3

TARGET = main
SRC = parser_demo.c parser.c

$(TARGET): $(SRC) parser.h
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm $(TARGET)

