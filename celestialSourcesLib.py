def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['celestialSources'], package = 'celestialSources')
	if env['PLATFORM'] == 'win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'celestialSources') 
	    env.Tool('findPkgPath', package = 'facilities') 

    env.Tool('genericSourcesLib')
    env.Tool('SpectObjLib')
    env.Tool('GRBLib')
    env.Tool('GRBobsLib')
    env.Tool('GRBtemplateLib')
    env.Tool('PulsarLib')
    env.Tool('eblAttenLib')
    env.Tool('microQuasarLib')
    env.Tool('EarthPhenomLib')
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'flux') 


def exists(env):
    return 1
