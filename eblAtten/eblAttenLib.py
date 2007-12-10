def generate(env, **kw):
    env.Tool('addLibrary', library = ['eblAtten'], package = 'celestialSources/eblAtten')

def exists(env):
    return 1
