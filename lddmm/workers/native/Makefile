#
#  segframe, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
#  https://github.com/nefan/segframe.git
# 
#  This file is part of segframe.
# 
#  segframe is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  segframe is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with segframe.  If not, see <http://www.gnu.org/licenses/>.
#  

include ../../../common.mk

.SUFFIXES: .mexa64 .c

HEADERS=common.h
SOURCES=fastPointPathOrder0.cpp fastPointPathOrder1.cpp fastPointPathEnergyOrder0.cpp fastPointGradTransportOrder0.cpp fastPointTransportOrder0.cpp fastPointGradTransportOrder1.cpp fastPointTransportOrder1.cpp
OBJECTS=$(SOURCES:.cpp=.mexa64)

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(OBJECTS)

$(OBJECTS): $(HEADERS)

fastPointGradTransportOrder1.mexa64: dkernelDerivativeMatrix.hpp gaussianKernels.hpp fastPointGradTransportOrder1.hpp order1.hpp

.c.mexa64:
	$(MEX) CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' -I'../../../$(INCLUDEPATH)' $(@:.mexa64=.c)
.cpp.mexa64:
	$(MEX) CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' -I'../../../$(INCLUDEPATH)' $(@:.mexa64=.cpp)
	
clean:
	rm -f ${OBJECTS}
