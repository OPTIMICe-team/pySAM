# Recorded script from Mayavi2
from numpy import array
from mayavi import mlab
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
vtk_file_reader = engine.open('/home/dori/develop/pySAM/vtk/_55_1.11595_1.70081e-05.vtk')
from mayavi.modules.surface import Surface
surface = Surface()
engine.add_filter(surface, vtk_file_reader)
mlab.savefig('scene.png')
