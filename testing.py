
from pyworkflow.em import Domain


plugins = Domain.getPlugins()

for k, v in Domain.getPlugins():
    print("plugin: ", k)

