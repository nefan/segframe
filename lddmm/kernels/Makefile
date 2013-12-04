#
# This file is part of networkred.
#
# Copyright (C) 2012, Technical University of Denmark
# https://github.com/nefan/networkred
#
# networkred is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# networkred is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with networkred.  If not, see <http://www.gnu.org/licenses/>.
# 

include ../../common.mk

SOURCES=gaussianTKsC.cpp
CLJMEX=gaussianTKsC.mexa64
HEADERS=gaussianKernels.hpp
OBJECTS=$(SOURCES:.cpp=.mexa64)

# libraries
CLIBS+=$(CXXLIBS)

CXXFLAGS+=$(CFLAGS)

INCLUDEPATH+=../include
INC=$(foreach d, $(INCLUDEPATH), -I'$d')

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $*.c -o $@
%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INC) $*.cpp -o $@

# header dependencies
#

%.cljmex.hpp: %.cljmex
	clojure -cp $(CLJMEXSRC) $*.cljmex > $@
$(CLJMEX): %.mexa64: %.cpp %.cljmex.hpp

$(OBJECTS): $(HEADERS) Makefile

%.mexa64: %.cpp
	$(MEX) $(MEXFLAGS) CFLAGS='$$CFLAGS $(CFLAGS) -DMEX' CXXFLAGS='$$CXXFLAGS $(CXXFLAGS) -DMEX' LDFLAGS='$$LDFLAGS $(LDFLAGS)' CXXLIBS='$$CXXLIBS $(CXXLIBS)' $(INC) $(@:.mexa64=.cpp)

clean:
	rm -f ${OBJECTS} *.cljmex.hpp
