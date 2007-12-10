def generate(env, **kw):
    env.Tool('addLibrary', library = ['Pulsar'], package = 'celestialSources/Pulsar')
    env.Tool('fluxLib')
    env.Tool('SpectObjLib')
    env.Tool('astroLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1
