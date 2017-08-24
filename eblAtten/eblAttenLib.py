def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['eblAtten'])
	if env['PLATFORM'] == 'win32' and env.get('CONTAINERNAME','') == 'GlastRelease':
	    env.Tool('findPkgPath', package = 'eblAtten') 
    #  No need for incsOnly section since eblAtten references no other pkgs
def exists(env):
    return 1
