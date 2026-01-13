CC      = mpicc
CFLAGS  = -O3
TARGET  = parallel_Thomas
SRC     = parallel_Thomas.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)