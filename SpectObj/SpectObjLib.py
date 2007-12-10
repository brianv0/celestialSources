def generate(env, **kw):
    env.Tool('addLibrary', library = ['SpectObj'], package = 'celestialSources/SpectObj')
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])

def exists(env):
    return 1
