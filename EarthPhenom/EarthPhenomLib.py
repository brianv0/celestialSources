#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['EarthPhenom'])
	if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'EarthPhenom') 
	    env.Tool('findPkgPath', package = 'flux') 
	    env.Tool('findPkgPath', package = 'astro') 
	    env.Tool('findPkgPath', package = 'facilities') 
    env.Tool('genericSourcesLib')
    env.Tool('fluxLib')
    env.Tool('SpectObjLib')
    env.Tool('astroLib')
    env.Tool('facilitiesLib')
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'facilities') 
        env.Tool('findPkgPath', package = 'celestialSources') 
        env.Tool('findPkgPath', package = 'flux') 
        env.Tool('findPkgPath', package = 'astro')
        env.Tool('findPkgPath', package = 'xmlBase')


def exists(env):
    return 1
