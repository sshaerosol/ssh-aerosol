import json
import sys

tag_ratio = True # Used only for GENOA output
invalid_line = '!!!\nKINETIC ARR 0 0 0'

all_keys = ['ARR1', 'ARR2', 'ARR3', 'MCM1', 'MCM2', 'MCM3', 'TROE4', 'TROE5', 'TROE7', 'TROE', 'SPEC']

mcm_keys = {'KRO2NO': '92 1',
            'KRO2HO2':'92 2',
            'KAPHO2':'92 3',
            'KAPNO':'92 4',
            'KRO2NO3':'92 5',
            'KNO3AL':'92 6',
            'KDEC':'92 7',
            'KROPRIM':'92 8',
            'KROSEC':'92 9',
            'KCH3O2':'92 10',
            'K298CH3O2':'92 11',
            'K14ISOM1':'92 12',
            'KFPAN': '93 21',
            'KBPAN': '93 22',
            'KBPPN': '93 23',
            'KRO2': '93 30'}
# EXTRA 99 [num] [ratio]
spec_keys = {
            '-10': '1',
            '-11': '2',
            '-12': '3',
            '-13': '4',
            '-14': '5',
            '-15': '6',
            '-16': '7',
            '-17': '8'
            }
                 
def check_inputs():
    if len(sys.argv) >= 3:
        arg1, arg3 = sys.argv[1], sys.argv[2]
        print(f'Read file name: {arg1}, new file name: {arg3}.')
        if len(sys.argv) >= 4:
            arg2 = sys.argv[3]
            if arg2 not in ['reaction', 'mol', 'all']:
                raise ValueError(f'Invalid file type: {arg2}. It shoud be reaction, mol or all.')
        else:
            if '.reactions' in arg1: arg2 = 'reaction'
            elif '.mol' in arg1: arg2 = 'mol'
            else: arg2 = 'all'
    else:
        print('This script requires at lease two arguments for file name and new file name.')
        print('Usage: python convert.py <file_name> <new_file_name> <file_type>')

    return arg1, arg2, arg3

def change_sign(val):
    return val[1:] if val.startswith('-') else '-' + val

def convert_reaction_file(old, new):
    """Convert old ssh-aerosol reaction list to the new one"""
    n = 0 # Count
    with open (old, 'r') as f: lines = f.read().splitlines()
    with open (new, 'w') as fout:
        fout.write(f'% Converted from {old} to {new}\n')
        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith('SET'): continue
            elif line.startswith('%'): fout.write(line + '\n')
            elif line.endswith('//'): fout.write(line.replace('//',''))
            elif line.startswith('KINETIC'):
                # No changes
                if 'PHOTOLYSIS' in line:
                    fout.write(line + '\n')
                    continue

                # RO2 pool
                if 'TB RO2' in line: line = line.replace('TB RO2','RO2 1')
                
                # Check keywords
                parts = [i for i in line.split(' ') if i != '']
                pkeys = [i for i in parts if i in all_keys]
                if len(pkeys) != 1: # Not find keyword or find multiple keywords
                    n += 1
                    print(f'# Can not convert kinetic rate: {line}, Need to convert manually.')
                    parts[0] = '%'+parts[0]
                    parts.append(invalid_line)
                else:
                    # Treat one keyword
                    key = pkeys[0]
                    ind = parts.index(key)
                    iratio = None # default ratio is 1
                    
                    if 'ARR' in key:
                        parts[ind] = 'ARR'
                        if key == 'ARR1': parts.extend(['0','0'])
                        elif key == 'ARR2': parts.insert(ind + 2, '0')
                    elif 'MCM' in key:
                        if key == 'MCM3': parts[ind] = 'EXTRA 91'
                        else:
                            ikey = lines[i-1].replace('%','').strip()
                            if '*' in ikey:
                                ikey, iratio = ikey.split('*',1)
                                print(f'Find ratios and key: {iratio}, {ikey}')
                            if 'KMT' in ikey:
                                val = int(ikey.replace('KMT',''))
                                if iratio is None: parts = parts[:ind] + [f"EXTRA 93 {val}"]
                                else: parts = parts[:ind] + [f"EXTRA 93 {val} {iratio}"]
                            elif ikey in mcm_keys: parts = parts[:ind] + [f'EXTRA {mcm_keys[ikey]}']
                            else: raise ValueError(f'Can not find MCM key {ikey}.')
                    elif 'TROE' in key:
                        parts[ind] = 'FALLOFF'
                        if key in ['TROE4', 'TROE5']:
                            # Change to negative: (T/300)^-C to (300/T)^C
                            parts[ind+2] = change_sign(parts[ind+2])
                            parts[ind+4] = change_sign(parts[ind+4])
                            # Add missing coefficients
                            parts.insert(ind+3, '0')
                            parts.insert(ind+6, '0')
                            if key == 'TROE4': parts.append('0.6')
                        else: # For TROE7 and TROE10 (10 never tested yet)
                            # Change to negative for C3 and C6: (T/300)^-C to (300/T)^C
                            parts[ind+3] = change_sign(parts[ind+3])
                            parts[ind+6] = change_sign(parts[ind+6])
                            # Change orders of C3 and C6
                            parts[ind+2], parts[ind+3] = parts[ind+3], parts[ind+2]
                            parts[ind+5], parts[ind+6] = parts[ind+6], parts[ind+5]
                        # Chage k0 (C1-C3) and kinf locations (C4-C6)
                        parts[ind+1:ind+3], parts[ind+4:ind+6] = parts[ind+4:ind+6], parts[ind+1:ind+3]
                        
                        if key != 'TROE10': parts.append('0 0 0')
                        else: print('Check TROE10 coefficients order!!!')
                    elif key == 'SPEC' and len(parts) == ind + 2 and parts[ind+1] in spec_keys:
                        parts[ind] = 'EXTRA 99'
                        parts[ind+1] = spec_keys[parts[ind+1]]
                    else:
                        n += 1
                        parts.append(invalid_line)
                        print(f'Not processed yet: {line}')

                    # Add ratio if need
                    if tag_ratio and 'ARR' not in key and iratio is None: parts.append('1.0')

                # Write
                fout.write(' '.join(parts) + '\n')
            else: fout.write(line+'\n')

    print(f'In total {n} kinetics need to convert manually. Check the reactions list with !!!.')

def convert_species_file(old, new):
    """Convert an old ssh-aerosol species mol file to the new format."""
    with open (old, 'r') as f: lines = f.read().splitlines()
    with open (new, 'w') as fout:
        for line in lines:
            if any(key in line for key in ['nameUsed', 'lump', 'jump', 'source', 'temp']): continue
            if 'condensed' in line: line = line.replace('condensed', 'condensable')
            elif 'SOAPStructure' in line and 'null' not in line:
                val = json.loads(line.split('\t',1)[1])
                soap = {}
                for i, j in enumerate(val[0]): soap[j] = val[1][i]
                line = f'SOAPStructure\t{json.dumps(soap)}'
            # Write       
            fout.write(line + '\n')
    

if __name__ == '__main__':

    filename, mode, new = check_inputs()
    
    if mode == 'reaction':
        convert_reaction_file(filename, f'{new}.reactions')
    elif mode == 'mol':
        convert_species_file(filename, f'{new}.mol')
    else:
        convert_reaction_file(f'{filename}.reactions', f'{new}.reactions')
        convert_species_file(f'{filename}.mol', f'{new}.mol')
