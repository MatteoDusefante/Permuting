CURRENT_PATH=.
PERM=permutation
#BEN=best_bucket
BEN=benchmark
#BEN=papi_test
SPMXV=spmxv

#BOOST_LIB=./boost_lib/
#BOOST_BIN=./boost_lib/
#BOOST_P=/usr/local/boost_1_60_0/
#BOOST_B=/usr/local/lib/

SOURCES_PERM=${PERM}.cpp
OBJECTS_PERM=$(SOURCES_PERM:%.cpp=%.o)

SOURCES_BEN=${BEN}.cpp
OBJECTS_BEN=$(SOURCES_BEN:%.cpp=%.o)

SOURCES_SPMXV=${SPMXV}.cpp
OBJECTS_SPMXV=$(SOURCES_SPMXV:%.cpp=%.o)

ifeq ($(shell uname),Darwin)
	CXX = g++-5
else
	CXX = g++
endif

CFLAGS = -Wall -Wextra
CFLAGS += -std=c++11
CFLAGS += -fopenmp
CFLAGS += -g #valgrind rows
CFLAGS += -lpapi
CFLAGS += -O3 -pipe -mfpmath=sse -march=native  -funroll-loops #-ffast-math  #-lm O2 #-funroll-all-loops

#TIME = #-lboost_timer
#LIB = #-I $(BOOST_P) #for boost
#BIN = #-L $(BOOST_B) #for boost binaries
#CFLAGS = -liomp5 -Wall -Wextra -O2 -ffast-math -fopenmp
#CCFLAGS = $(DFLAGS) $(TIME) #$(BIN) #$(LIB)


all: compile_ben #compile_ben #compile_spmxv

compile_perm: $(OBJECTS_PERM)
	$(CXX) -o $(PERM) $(CFLAGS) $(OBJECTS_PERM:%=$(CURRENT_PATH)/obj/$(PERM).o)

$(OBJECTS_PERM): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $(CURRENT_PATH)/obj/$(PERM).o

compile_ben: $(OBJECTS_BEN)
	$(CXX) -o $(BEN) $(CFLAGS) $(OBJECTS_BEN:%=$(CURRENT_PATH)/obj/$(BEN).o)

$(OBJECTS_BEN): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $(CURRENT_PATH)/obj/$(BEN).o


compile_spmxv: $(OBJECTS_SPMXV)
	$(CXX) -o $(SPMXV) $(CFLAGS) $(OBJECTS_SPMXV:%=$(CURRENT_PATH)/obj/$(SPMXV).o)

$(OBJECTS_SPMXV): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $(CURRENT_PATH)/obj/$(SPMXV).o

.PHONY: clean

clean:
	rm -f obj/*.o
	rm ${BEN}
