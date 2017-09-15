# GCC options
#CC = gcc-6.3.0
#CFLAGS = -Ofast 

# Intel icc compiler
CC = icc
CFLAGS = -restrict -Ofast -std=c99 -pedantic

# Clang options
#CC = clang
#CFLAGS = -Ofast 

LDFLAGS = 


SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJ) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
