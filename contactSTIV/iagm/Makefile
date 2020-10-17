
# Things for users to define
#
include make.inc


IAGM_LIB := libiagm.$(LIB_EXT)
all: $(IAGM_LIB)

CXXFLAGS := -std=c++0x -Wall -Wextra -pipe -fPIC -O3
CXX := g++
INCS := -Iinc -Isrc -I$(EIGEN_INC)
LIBS := 

VPATH := ./src
SRCS := $(wildcard ./src/*.cpp)

BUILD_DIR := build/
OBJS := $(addprefix $(BUILD_DIR), $(addsuffix .o, $(basename $(notdir $(SRCS)))))


$(IAGM_LIB): $(OBJS)
	mkdir -p lib
#	$(CXX) -dynamiclib -Wl,-install_name,./lib/$(IAGM_LIB) -o lib/$(IAGM_LIB) $(OBJS) $(LIBS)
	$(CXX) -shared -Wl,-install_name,lib/$(IAGM_LIB) -o lib/$(IAGM_LIB) $(OBJS) $(LIBS)
#	cp -rf lib ../../../examples/

$(BUILD_DIR)%.o: %.cpp
	$(shell mkdir -p $(BUILD_DIR))
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<

depend:
	$(CXX) $(INCS) -MM $(SRCS) | sed 's/.*\.o/build\/&/' > .depend

ifeq (.depend, $(wildcard .depend))
  include .depend
endif

clean: 
	-rm -f lib/$(IAGM_LIB) $(BUILD_DIR)*.o .depend

rebuild: clean depend all

