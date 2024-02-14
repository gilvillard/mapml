INC_DIR = ../../univ-linalg/flint-extras/include
LIB_DIR = ../../univ-linalg/flint-extras/lib
OBJS = $(wildcard $(LIB_DIR)/*)

#$(info $$OBJS is [${OBJS}])

#SOURCES = $(wildcard *.c)
#HEADERS = $(wildcard ../*.h)
#OBJS = $(patsubst %.c, $(LIB_DIR)/%.o, $(SOURCES))
#INSTALLED_HEADERS = $(patsubst ../%.h, $(INC_DIR)/%.h, $(HEADERS))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I/opt/homebrew/include -I/usr/local/include 
LIBS = -lflint -lmpfr -lgmp -lm -lpthread
LDFLAGS = -L/opt/homebrew/lib -L/usr/local/lib

.SUFFIXES: .o .c

.PHONY: clean


.c:
	$(CC) -o $* $*.c $(OBJS) $(CFLAGS) $(LDFLAGS) $(LIBS)

