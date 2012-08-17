def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['genericSources'])
	if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'genericSources') 
            env.Tool('findPkgPath', package = 'facilities') 
            env.Tool('findPkgPath', package = 'flux') 
            env.Tool('findPkgPath', package = 'astro')
            env.Tool('findPkgPath', package = 'xmlBase')
    env.Tool('facilitiesLib')
    env.Tool('astroLib')
    env.Tool('fluxLib')
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'facilities') 
        env.Tool('findPkgPath', package = 'celestialSources') 
        env.Tool('findPkgPath', package = 'eblAtten')
        env.Tool('findPkgPath', package = 'flux') 
        env.Tool('findPkgPath', package = 'astro')
        env.Tool('findPkgPath', package = 'xmlBase')
        env.Tool('findPkgPath', package = 'tip')

def exists(env):
    return 1
