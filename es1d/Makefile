# GCC options
CC = gcc
CFLAGS = -Ofast -std=c99 -pedantic

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic

# Clang options
# CC = clang
# CFLAGS = -Ofast -std=c99 -pedantic

LDFLAGS =


SOURCE = charge.c field.c particles.c grid.c fft.c random.c timer.c main.c simulation.c \
         zdf.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
