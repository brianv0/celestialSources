def generate(env, **kw):
    env.Tool('addLibrary', library = ['microQuasar'], packages = 'celestialSources/microQuasar')
    env.Tool('fluxLib')
    env.Tool('SpectObjLib')
    env.Tool('astroLib')
    env.Tool('facilitiesLib')

def exists(env):
    return 1
