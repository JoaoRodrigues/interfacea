# interfacea - interface energy analysis


## Quick example
```
import interfacea as ia

mol = ia.read('example.pdb')
mol.add_missing_atoms()

ff = ia.get_forcefield('amber14')
mol.forcefield(ff)
mol.minimize()

results = mol.analyze(level='residue')
results.write('energetics.json')
```