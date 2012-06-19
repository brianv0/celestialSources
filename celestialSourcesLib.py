def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['celestialSources'], package = 'celestialSources')
    env.Tool('genericSourcesLib')
    env.Tool('SpectObjLib')
    env.Tool('GRBLib')
    env.Tool('GRBobsLib')
    env.Tool('GRBtemplateLib')
    env.Tool('PulsarLib')
    env.Tool('eblAttenLib')
    env.Tool('microQuasarLib')
    env.Tool('EarthPhenomLib')

def exists(env):
    return 1
