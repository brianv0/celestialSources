def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['eblAtten'])
        if (env['PLATFORM'] == 'win32' and
            env.get('CONTAINERNAME','') == 'GlastRelease'):
	    env.Tool('findPkgPath', package = 'eblAtten') 
            env.Tool('findPkgPath', package = 'facilities') 
            env.Tool('findPkgPath', package = 'st_facilities')
    env.Tool('facilitiesLib')
    env.Tool('st_facilitiesLib')
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'facilities') 
        env.Tool('findPkgPath', package = 'st_facilities') 

def exists(env):
    return 1
