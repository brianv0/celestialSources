#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['EarthPhenom'])
    env.Tool('genericSourcesLib')
    env.Tool('fluxLib')
    env.Tool('SpectObjLib')
    env.Tool('astroLib')
    env.Tool('facilitiesLib')

def exists(env):
    return 1
