def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['Pulsar'])
    env.Tool('fluxLib')
    env.Tool('SpectObjLib')
    env.Tool('astroLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1
