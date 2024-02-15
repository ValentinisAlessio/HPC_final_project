CC = mpicc
INCDIR = headers
SRCDIR = source
OBJDIR = objects
DEPS = $(INCDIR)/par_quicksort.h

CFLAGS = -I$(INCDIR) -c -fopenmp
OFLAGS = -O3 -march=native -fopenmp

OBJECTS = $(OBJDIR)/par_quicksort.o $(OBJDIR)/main.o

main: $(OBJECTS)
	$(CC) $(OFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJDIR)/*.o main