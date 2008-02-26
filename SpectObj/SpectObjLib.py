def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['SpectObj'])
    env.Tool('eblAttenLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])

def exists(env):
    return 1
