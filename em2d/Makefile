# GCC options
#CC = gcc
#CFLAGS = -Ofast -std=c99 -pedantic
#LDFLAGS = -lm

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic
#LDFLAGS =

# Clang options
CC = clang
CFLAGS = -Ofast -std=c99 -Wall -pedantic
#CFLAGS = -Og -g -std=c99 -Wall -pedantic
#LDFLAGS = -lm

SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
