## OS specific stuff - don't change ##
ifeq ($(OS),Windows_NT)
TARGET = himeno.exe
RM = del
RUN = $(TARGET)
else
TARGET = himeno
RM = rm
RUN = ./$(TARGET)
endif


## Makefile code ##
CC = gcc

CFLAGS = -lm -Wall -Wextra -fmax-errors=10 -pthread
CSTD = --std=gnu11
JACOBI_FLAGS = -O3

SRCS = timer.c himeno_contest.c main.c
OBJS = $(SRCS:.c=.o)

default: $(TARGET)

run: $(TARGET)
	$(RUN) l

$(TARGET): jacobi.o $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) jacobi.o $(OBJS)

jacobi.o: jacobi.c
	$(CC) -c jacobi.c $(JACOBI_FLAGS) $(CFLAGS) $(CSTD) 

%.o: %.c
	$(CC) -c $<  -o $@ $(JACOBI_FLAGS) $(CFLAGS) $(CSTD) -O0

clean:
	$(RM) $(TARGET) jacobi.o $(OBJS)
