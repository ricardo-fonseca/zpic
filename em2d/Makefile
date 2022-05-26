# GCC options
CC = gcc
CFLAGS = -Ofast -std=c99 -Wall -pedantic
LDFLAGS = -lm

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic
#LDFLAGS =

# Clang options
#CC = clang
#CFLAGS = -Ofast -std=c99 -Wall -pedantic
#CFLAGS = -Og -g -std=c99 -Wall -pedantic
#LDFLAGS = -lm

SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c

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
