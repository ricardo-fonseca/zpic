# GCC options
CC = gcc
CFLAGS = -Ofast -std=c99 -pedantic -Wall
LDFLAGS = -lm

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic -Wall
#LDFLAGS =

# Clang options
#CC = clang
#CFLAGS = -Ofast -std=c99 -pedantic -Wall
#LDFLAGS = -lm


SOURCE = filter.c charge.c current.c emf.c particles.c grid.c \
         fft.c random.c timer.c main.c simulation.c zdf.c

TARGET = zpic

DOCSBASE = docs

DOCS = $(DOCSBASE)/html/index.html

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

docs : $(DOCS)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

$(DOCS) : $(SOURCE)
	@doxygen ./Doxyfile

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
	rm -rf $(DOCSBASE)
