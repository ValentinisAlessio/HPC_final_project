CC = mpicc
INCDIR = headers
SRCDIR = source
OBJDIR = objects
DEPS = $(INCDIR)/par_quicksort.h

CFLAGS = -c -I$(INCDIR)
OFLAGS = -Wall -Wextra -O3 -march=native -fopenmp

OBJECTS = $(OBJDIR)/par_quicksort.o $(OBJDIR)/main.o

main: $(OBJECTS)
	$(CC) $(OFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c $(OFLAGS) $< -o $@

clean:
	rm -f $(OBJDIR)/*.o main