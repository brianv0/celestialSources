def generate(env, **kw):
    env.Tool('addLibrary', library = ['GRBobs'], package = 'celestialSources/GRBobs')
    env.Tool('fluxLib')
    env.Tool('astroLib')
    env.Tool('SpectObjLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1
