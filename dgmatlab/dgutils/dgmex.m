function dgmex(cfile,verbose)

if nargin<2, verbose=''; end

cmexext=mexext;
ARCH3DG=cmexext(4:end);
iswindows=(ARCH3DG(1)=='w');

% Assume this file is {p0}/dgmatlab/dgutils/dgmex.m
p0=fileparts(fileparts(fileparts(mfilename('fullpath'))));
include=[p0,'/include'];
lib=[p0,'/lib/',ARCH3DG,'/libdgmex.a'];

include_hdf5='/usr/lib/x86_64-linux-gnu/hdf5/serial/lib';
lib_hdf5='/usr/lib/x86_64-linux-gnu/hdf5/serial';

if iswindows
  add1=[p0,'/os/winlibs/libumfpack.a'];
  add2=[p0,'/os/winlibs/libamd.a'];
  add3=fullfile(matlabroot,'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
  add4=fullfile(matlabroot,'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
  mex(['-I',include],cfile,lib,add1,add2,add3,add4);
else
  opt='CXXOPTIMFLAGS=-O3 -std=c++0x -pedantic -mtune=native -funroll-loops -ftree-vectorize -DNDEBUG';
  add={'-lhdf5_cpp','-lhdf5_hl','-lhdf5','-lumfpack','-lamd','-llapack','-lblas'};
  add0 = {['-I',include],['-I',include_hdf5],
          ['-L',lib_hdf5], ...
         '-L/usr/lib/x86_64-linux-gnu'};
  if ~isempty(verbose)
    mex(opt,verbose,add0{:},'-largeArrayDims',cfile,lib,add{:});
  else
    mex(opt,add0{:},'-largeArrayDims',cfile,lib,add{:});
  end
end
