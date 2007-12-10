import glob,os,platform

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

celestialSourcesLib = libEnv.StaticLibrary('celestialSources', listFiles(['src/*.cxx']))

progEnv.Tool('celestialSourcesLib')
test_celestialSourcesBin = progEnv.Program('test_celestialSources', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'celestialSources', libraries = [celestialSourcesLib], testApps = [test_celestialSourcesBin], includes = listFiles(['celestialSources/*.h']))
