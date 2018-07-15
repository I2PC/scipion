
from pyworkflow.em import Domain


plugins = Domain.getPlugins()

for k, v in plugins.iteritems():
    print("plugin: ", k)

