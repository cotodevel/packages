TARGET = TGDSPKGBuilder
LIBS = -lm	-lz
CC = gcc
CFLAGS = -g -Wall

.PHONY: default all clean

default: installzlib $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))	$(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

installzlib:
	-@echo 'setup zlib'
	-@7z x $(CURDIR)/zlib-1.2.11.zip	-o$(CURDIR)
	-@cd	$(CURDIR)/zlib-1.2.11 &&	./configure --prefix=/usr/local/zlib
	-@sudo $(MAKE)	-C	$(CURDIR)/zlib-1.2.11
	-@sudo mv	$(CURDIR)/zlib-1.2.11/libz.a	$(CURDIR)

%.o: %.cpp $(HEADERS)
	 g++  $(CPPFLAGS) -c $< -o $@	-static -static-libstdc++

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@	-static -static-libgcc

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	-@g++	$(OBJECTS)	-Wall	$(LIBS)	-o	$@
	-@sudo mv	$(CURDIR)/$(TARGET)	$(GCC_BUILD_ENV)$(GCC_BIN_PATH)$(TARGET)

clean:
	-rm -f *.o libz.a	$(TARGET)
	-@sudo rm	-fr	$(CURDIR)/zlib-1.2.11