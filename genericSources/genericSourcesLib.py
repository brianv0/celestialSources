def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['genericSources'])
    env.Tool('facilitiesLib')
    env.Tool('astroLib')
    env.Tool('fluxLib')
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1
