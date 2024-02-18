MAPLE_INC = /Library/Frameworks/Maple.framework/Versions/2023/extern/include
MAPLE_LIBDIR = /Library/Frameworks/Maple.framework/Versions/2023/bin.APPLE_ARM64_MACOS

INC_DIR = $(HOME)/univ-linalg/flint-extras/include

#SOURCES = $(wildcard *.c)
#HEADERS = $(wildcard ../*.h)
#INSTALLED_HEADERS = $(patsubst ../%.h, $(INC_DIR)/%.h, $(HEADERS))
OBJS =  $(wildcard  /Users/gvillard/univ-linalg/flint-extras/lib/*.o)

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I$(MAPLE_INC) -I/opt/homebrew/include -I/usr/local/include 

INSTALL_DIR = /usr/local/lib

# ! With setenv DYLD_LIBRARY_PATH=/usr/local/lib for pmlflint.dylib

#LIBS = -lflint -lmpfr -lgmp -lm -lpthread -lpmlflint -lmaplec
LIBS = -lflint -lmpfr -lgmp -lm -lpthread  -lmaplec
LDFLAGS = -L/opt/homebrew/lib -L/usr/local/lib -L$(MAPLE_LIBDIR) 

.SUFFIXES: .o .c

.PHONY: clean

.c:
	$(CC) -c $*.c  $(CFLAGS) 


%.o: %.c 
	$(CC) -c $*.c  $(CFLAGS)

shared: dum.o
	$(CC) -shared -fPIC -o dum.dylib dum.o $(OBJS) $(CFLAGS)  $(LDFLAGS) $(LIBS)
	#$(CC) -shared -fPIC -o dum.dylib dum.o $(CFLAGS)  $(LDFLAGS) $(LIBS)

install:
	cp dum.dylib $(INSTALL_DIR)/



