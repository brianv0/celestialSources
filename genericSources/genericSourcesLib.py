def generate(env, **kw):
    env.Tool('addLibrary', library = ['genericSources', 'g2c'], package = 'celestialSources/genericSources')
    env.Tool('facilitiesLib')
    env.Tool('astroLib')
    env.Tool('fluxLib')
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1
