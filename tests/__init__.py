"""General parameters for tests."""
import pathlib


rootdir = pathlib.Path(".")
datadir = rootdir / "tests" / "data"

defaultmdl = str(datadir / 'pdb' / 'default.pdb')
multimodel = str(datadir / 'pdb' / 'multimodel.pdb')
