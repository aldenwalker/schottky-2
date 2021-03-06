UNAME = $(shell uname -s)
ifeq ($(UNAME),Darwin)
  CC=clang++
else
  CC=g++
endif
CFLAGS=-g -Wall -Wextra -pedantic
IFLAGS=-I/usr/X11R6/include -I/opt/X11/include
LFLAGS=-L/usr/X11R6/lib -L/opt/X11/lib -lX11
all: schottky

graphics.o: graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

schottky.o: schottky.cc ifs.cc ifs_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c schottky.cc

trap_grid.o: trap_grid.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trap_grid.cc
	
movie.o: movie.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c movie.cc

ifs_gui.o: ifs_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs_gui.cc

ifs.o: ifs.cc ifs_draw.cc ifs_trap.cc ifs_interface.cc ifs_connected.cc ifs_trap_like.cc ifs_set_A.cc ifs_set_B.cc ifs_nifs.cc ifs_gifs.cc ifs_2d.cc ifs_picture.cc ifs_boundary.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs.cc

schottky: schottky.o graphics.o ifs.o trap_grid.o movie.o ifs_gui.o
	$(CC) $(CFLAGS) -o schottky schottky.o graphics.o trap_grid.o movie.o ifs.o ifs_gui.o $(LFLAGS) -lm

clean:
	rm *.o
	rm schottky
