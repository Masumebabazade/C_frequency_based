CC=gcc
CFLAGS=-O2 -std=c11 -fopenmp
TARGET=freqprog

all: $(TARGET)

$(TARGET): main.c
	$(CC) $(CFLAGS) main.c -o $(TARGET) -lm

clean:
	rm -f $(TARGET)
