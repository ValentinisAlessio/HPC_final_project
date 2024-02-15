CC = mpicc
INCDIR = headers
SRCDIR = source
OBJDIR = objects
DEPS = $(INCDIR)/par_quicksort.h

CFLAGS = -I$(INCDIR)
OFLAGS = -O3 -march=native -fopenmp

OBJECTS = $(OBJDIR)/par_quicksort.o $(OBJDIR)/main.o

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@

clean:
	rm -f $(OBJDIR)/*.o main