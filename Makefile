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

include common.mk

native := lddmm/workers/native
utils := utils
startup = startup.m

DIRS := native utils
DOCSDIR=docs

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(native) $(startup)

$(native):
	$(MAKE) --directory=$@

$(startup):
	echo "addpath(genpath('registration'))" > startup.m
	echo "addpath(genpath('lddmm'))" >> startup.m
	echo "addpath(genpath('utils'))" >> startup.m
	echo "addpath(genpath('tests'))" >> startup.m
	echo "addpath(genpath('thirdparty'))" >> startup.m
	echo "addpath(genpath('$(MINFUNC)'))" >> startup.m

doc:
	$(MATLAB) -nodesktop -nosplash -r "addpath(genpath('thirdparty')); addpath('utils'); gendoc; exit;"

clean:
	rm -rf $(DOCSDIR)
	rm $(startup)
	$(MAKE) --directory=$(native) clean

.PHONY: all $(native)
