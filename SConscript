# -*- python -*-
#
# $Header$
# Authors: James Chiang <jchiang@slac.stanford.edu>
# Version: celestialSources-01-03-03
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('celestialSourcesLib', depsOnly = 1)
celestialSourcesLib = libEnv.StaticLibrary('celestialSources',
                                           listFiles(['src/*.cxx']))

progEnv.Tool('celestialSourcesLib')
test_celestialSourcesBin = progEnv.Program('test_celestialSources',
                                           listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', package='celestialSources',
             staticLibraryCxts=[[celestialSourcesLib, libEnv]], 
             testAppCxts=[[test_celestialSourcesBin, progEnv]],
             includes=listFiles(['celestialSources/*.h']))




