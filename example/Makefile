INC_DIR = $(HOME)/univ-linalg/flint-extras/include

#SOURCES = $(wildcard *.c)
#HEADERS = $(wildcard ../*.h)
#INSTALLED_HEADERS = $(patsubst ../%.h, $(INC_DIR)/%.h, $(HEADERS))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I/opt/homebrew/include -I/usr/local/include 

# ! With setenv DYLD_LIBRARY_PATH=/usr/local/lib for pmlflint.dylib

LIBS = -lflint -lmpfr -lgmp -lm -lpthread -lpmlflint
LDFLAGS = -L/opt/homebrew/lib -L/usr/local/lib 

.SUFFIXES: .o .c

.PHONY: clean


.c:
	$(CC) -o $* $*.c  $(CFLAGS) $(LDFLAGS) $(LIBS)

