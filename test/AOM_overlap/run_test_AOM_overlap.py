#!/usr/bin/env python3

import os

single_molecules={
	'ethylene':{
		'src':'../STO_proj/ethylene.xyz',
		'atoms':6,
		'aom_file':'../STO_proj/AOM_COEFF_single_ethylene_MO_6.include',
		'direction':'y',
	},
	'cyclopropene':{
		'src':'../STO_proj/cyclopropene.xyz',
		'atoms':7,
		'aom_file':'../STO_proj/AOM_COEFF_single_cyclopropene_MO_8.include',
		'direction':'y',
	},
	'cyclobutadiene':{
		'src':'../STO_proj/cyclobutadiene.xyz',
		'atoms':8,
		'aom_file':'../STO_proj/AOM_COEFF_single_cyclobutadiene_MO_10.include',
		'direction':'z',
	},
	'cyclopentadiene':{
		'src':'../STO_proj/cyclopentadiene.xyz',
		'atoms':11,
		'aom_file':'../STO_proj/AOM_COEFF_single_cyclopentadiene_MO_13.include',
		'direction':'y',
	},
	'furane':{
		'src':'../STO_proj/furane.xyz',
		'atoms':9,
		'aom_file':'../STO_proj/AOM_COEFF_single_furane_MO_13.include',
		'direction':'x',
	},
	'pyrrole':{
		'src':'../STO_proj/pyrrole.xyz',
		'atoms':10,
		'aom_file':'../STO_proj/AOM_COEFF_single_pyrrole_MO_13.include',
		'direction':'x',
	},
	'thiophene':{
		'src':'../STO_proj/thiophene.xyz',
		'atoms':9,
		'aom_file':'../STO_proj/AOM_COEFF_single_thiophene_MO_13.include',
		'direction':'x',
	},
	'imidazole':{
		'src':'../STO_proj/imidazole.xyz',
		'atoms':9,
		'aom_file':'../STO_proj/AOM_COEFF_single_imidazole_MO_13.include',
		'direction':'z',
	},
	'phenol':{
		'src':'../STO_proj/phenol.xyz',
		'atoms':13,
		'aom_file':'../STO_proj/AOM_COEFF_single_phenol_MO_18.include',
		'direction':'z',
	},
}

distance=[3.5,4.0,4.5,5.0]

mode='dimer'

verb='none'

STO_exponents=[
'AOM_mu	H	1.0000\n',
'AOM_mu	C	1.6083	1.0300\n',
'AOM_mu	N	1.9237	1.2050\n',
'AOM_mu	O	2.2458	2.2266\n',
'AOM_mu	S	2.1223	1.5850\n'
]

aug={'x':[1,0,0],'y':[0,1,0],'z':[0,0,1]}
for molecule,values in single_molecules.items():
	aom_filename=f'AOM_COEFF_{molecule}.include'
	fp=open(values['aom_file'],mode='r')
	aom=fp.readlines()
	fp.close()
	with open(aom_filename,mode='w') as fp:
		for line in aom:
			print(line,file=fp,end='')
		for line in aom:
			print(line,file=fp,end='')
	for d in distance:
		filename=f'config_{molecule}_{d}.txt'
		dimer_filename=f'dimer_{molecule}_{d}.xyz'
		with open(filename,mode='w') as fp:
			print(f'mode {mode}',file=fp)
			print(f'verb {verb}',file=fp)
			print(f'dimer {dimer_filename}',file=fp)
			print(f'AOM_include {aom_filename}',file=fp)
			print(f'atoms_frag1 {values["atoms"]}',file=fp)
			print('#-----------------------------------------------------------------------------------------',file=fp)
			print(f'{"".join(STO_exponents)}',file=fp)
		fp=open(values['src'],mode='r')
		xyz=fp.readlines()
		fp.close()
		atoms=int(xyz[0])
		species=[i.split()[0] for i in xyz[2:atoms+2]]
		x,y,z=[[float(i.split()[j]) for i in xyz[2:atoms+2]] for j in [1,2,3]]
		with open(dimer_filename,mode='w') as fp:
			print(f'{atoms*2}\n',file=fp)
			for i in range(atoms):
				print(f'{species[i]} {x[i]} {y[i]} {z[i]}',file=fp)
			for i in range(atoms):
				print(f'{species[i]} {x[i]+d*aug[values["direction"]][0]} {y[i]+d*aug[values["direction"]][1]} {z[i]+d*aug[values["direction"]][2]}',file=fp)
		os.system(f'../../src/STO_proj_AOM_overlap {filename}')
	

