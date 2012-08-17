def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['SpectObj'])
	if env['PLATFORM'] == 'win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'SpectObj') 
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'eblAtten') 
        env.Tool('findPkgPath', package = 'facilities') 


def exists(env):
    return 1
