PROJECT := aggregation
CC := "clang++"
CC_OPTS := -g -O2
INCLUDE_DIRS := external/inc
LINKER_OPTS := -lX11 -lc -lm

$(PROJECT): $(PROJECT).cpp
	$(CC) $(CC_OPTS) -I$(INCLUDE_DIRS) -o $(PROJECT) $(PROJECT).cpp $(LINKER_OPTS)

clean:
	rm $(PROJECT)
