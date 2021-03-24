import os
import json

with open('overlap.json') as fp:
    ref_data=json.load(fp)

AOM_dict = {
    'H':{'1s':1.0000},
    'C':{'2s':1.6083,'2p':1.385600},
    'N':{'2s':1.9237,'2p':1.617102},
    'O':{'2s':2.2458,'2p':1.505135},
    'F':{'2s':2.5638,'2p':1.665190},
    'S':{'3s':2.1223,'3p':1.641119},
}

res = {}
target_molecules = sorted([i for i in os.listdir('dimers') if not i.startswith('.')])
for mol in target_molecules:
    prefix = mol.split('_')[0] + '_'
    molecule = mol.split(prefix)[-1]
    res[molecule] = {}
    aom_coeff = [i for i in os.listdir('aom_coeff') if i.find(f'AOM_COEFF_dimer_{molecule}_MO') == 0][0]
    with open(f'{molecule}_input.txt', mode = 'w') as fp:
        print('mode\tdimer_simplex', file = fp)
        print('verb\tnone', file = fp)
        dimer_list=sorted([i for i in os.listdir(f'dimers/{mol}') if not i.startswith('.')])
        print(f'entries\t{len(dimer_list)}', file = fp)
        for dimer in dimer_list:
            xyz = f'{"dimers"}/{mol}/{dimer}'
            with open(xyz) as fp_:
                raw = fp_.readlines()
            atoms=int(raw[0])
            print(f'{molecule}\t{xyz}\t{"aom_coeff"}/{aom_coeff}\t{int(atoms/2)}\tN/A', file = fp)
        for species, STO in AOM_dict.items():
            print(f'AOM_mu\t{species}\t' + '\t'.join(list(map(str,STO.values()))), file = fp)
        print('simplex_steps\t0', file = fp)
    stream = os.popen(f'../../src/STO_proj_AOM_overlap {molecule}_input.txt')
    output = stream.readlines()
    output = output[-1].strip().split(';')[:-1]
    for c,i in enumerate(dimer_list):
        res[molecule][i] = {}
        res[molecule][i]['value'] = float(output[c])

atol = 1.0e-6
total_dimers = 0
total = len(ref_data)
passed = 0
for mol_counter, (molecule, record) in enumerate(ref_data.items()):
    test_status = 'PASS'
    for dimer, value in record.items():
        total_dimers += 1
        if abs(ref_data[molecule][dimer] - res[molecule][dimer]['value']) < atol:
            res[molecule][dimer]['status'] = 'PASS'
        else:
            res[molecule][dimer]['status'] = 'FAIL'
            test_status = 'FAIL'
    if test_status == 'PASS':
        passed += 1
        print(f'[{mol_counter + 1}/{total}] PASS\t{molecule}')
    else:
        print(f'[{mol_counter + 1}/{total}] ! FAIL\t{molecule}')
print(f'Total tests: {total}; passed: {passed}/{total}; failed {total-passed}/{total}')
print(f'Total dimers: {total_dimers}')

# clean up
for mol in target_molecules:
    prefix = mol.split('_')[0] + '_'
    os.system(f'rm {mol.split(prefix)[-1]}_input.txt')

