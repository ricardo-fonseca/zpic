# GCC options
CC = gcc
CFLAGS = -Ofast -std=c99 -pedantic

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic -Wall

# Clang options
#CC = clang
#CFLAGS = -Ofast -std=c99 -pedantic -Wall

LDFLAGS =


SOURCE = charge.c current.c emf.c particles.c grid.c fft.c random.c timer.c main.c simulation.c \
         zdf.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC)    $(CFLAGS) $(LDFLAGS) $(OBJ) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
