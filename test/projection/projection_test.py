import os
import json
import time
import ast
import numpy as np

# executable relative path
path_to_bin = '../../src/STO_proj_AOM_overlap'

# STO parameters for GTO-STO projection step
STO_proj_dict = {
    'H':{'1s':1.0000},
    'C':{'2s':1.6083,'2p':1.442657},
    'N':{'2s':1.9237,'2p':1.646703},
    'O':{'2s':2.2458,'2p':1.858823},
    'F':{'2s':2.5638,'2p':2.136394},
    'S':{'3s':2.1223,'3p':1.651749},
}

# common parameters for input files
basis = 'DZVP-GTH'
basis_file = 'cp2k/GTH_BASIS_SETS'
MO_channel = 1
cube = 'no'
cube_grid = 0.25

# test parameters
# np.allclose params
rtol = 1.0e-3
atol = 1.0e-6
# leave empty for all molecules or specify, e.g. ['84a', 'n11']
target = []

# -------------------------------------------------------------------

# read json files
with open('single_molecules.json') as fp:
    ref_data = json.load(fp)
with open('test_configuration.json') as fp:
    test_configuration = json.load(fp)

# run tests!
total = len(target) if target != [] else len(ref_data)
total_time = 0
counter = 0
passed = 0
for key, value in test_configuration.items():
    if key in target or target == []:
        counter += 1
        config = []
        config.append('mode\tmolecule')
        config.append('verb\tnone')
        config.append(f'name\t{key}')
        config.append(f'molecule\t{value["xyz"]}')
        config.append(f'basis\t{basis}')
        config.append(f'basis_file\t{basis_file}')
        config.append(f'MO_file\t{value["MOLog"]}')
        config.append(f'MO\t{value["MO"]}')
        config.append(f'MO_channel\t{MO_channel}')
        config.append(f'cube\t{cube}')
        config.append(f'cube_grid\t{cube_grid}')
        for species, STO in STO_proj_dict.items():
            config.append(f'proj_mu\t{species}\t' + '\t'.join(list(map(str,STO.values()))))
        with open(f'{key}_projection.txt', mode = 'w') as fp:
            for i in config:
                print(i, file = fp)
        tic = time.perf_counter()
        os.system(f'{path_to_bin} {key}_projection.txt')
        toc = time.perf_counter()
        test_time = toc - tic
        total_time += test_time
        with open(f'{key}_state.dat') as fp:
            raw = fp.read()
        test_data = ast.literal_eval(raw)

        check=[]
        check.append(np.allclose(test_data['pvecs']['px'],ref_data[key]['pvecs']['px'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['pvecs']['py'],ref_data[key]['pvecs']['py'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['pvecs']['pz'],ref_data[key]['pvecs']['pz'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['S_matrix'],ref_data[key]['S_matrix'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['V_array'],ref_data[key]['V_array'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['STO_matrix'],ref_data[key]['STO_matrix'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data['singular_values'],ref_data[key]['singular_values'],rtol=rtol, atol=atol))
        check.append(np.allclose(list(test_data['AOM_dict'].values()),list(ref_data[key]['AOM_dict'].values()),rtol=rtol, atol=atol))
        check.append(np.allclose(list(test_data['compl_dict'].values()),list(ref_data[key]['compl_dict'].values()),rtol=rtol, atol=atol))
        
        print(f'[{counter}/{total}] ', end = '')
        if check == [True for i in check]:
            print(f'PASS\t{key}\t{test_time:.2f}s')
            passed += 1
        else:
            print(f'! FAIL\t{key}\t{test_time:.2f}s')
print(f'Total tests: {total}; passed: {passed}/{total}; failed {total-passed}/{total}')
print(f'Execution time: {total_time:.2f}s')

# clean up
for key, value in test_configuration.items():
    if key in target or target == []:
        os.system(f'rm {key}_state.dat')
        os.system(f'rm {key}_projection.txt')
        os.system(f'rm AOM_*{key}*.include')
        os.system(f'rm log*{basis}*dat')

