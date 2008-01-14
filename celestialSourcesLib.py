def generate(env, **kw):
    env.Tool('addLibrary', library = ['celestialSources'], package = 'celestialSources')
    env.Tool('genericSourcesLib')
    env.Tool('SpectObjLib')
    env.Tool('GRBLib')
    env.Tool('GRBobsLib')
    env.Tool('GRBtemplateLib')
    env.Tool('PulsarLib')
    env.Tool('eblAttenLib')
    env.Tool('microQuasarLib')

def exists(env):
    return 1