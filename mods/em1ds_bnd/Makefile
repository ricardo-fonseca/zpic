# GCC options
#CC = gcc-6.3.0
#CFLAGS = -Ofast -std=c99 -pedantic 
#LDFLAGS = -lm

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic
#LDFLAGS =

# Clang options
CC = clang
CFLAGS = -Ofast -std=c99 -pedantic
LDFLAGS = -lm


SOURCE = charge.c current.c emf.c particles.c grid.c fft.c random.c timer.c main.c simulation.c

# HDF5 output
#SOURCE += hdf5.c
#CFLAGS += -DHDF5 -I/opt/hdf5/1.8.16/include/ 
#LDFLAGS += -L/opt/hdf5/1.8.16/lib -lhdf5

# ZDF output
# No external libraries required
SOURCE += zdf.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -I. -c $(CFLAGS) $< -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
