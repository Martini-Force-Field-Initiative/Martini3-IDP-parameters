#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:12:29 2023
merge folded and IDR topol, and add restraint for linking residues;
have to manually edit the linking residue protonation.
also have to manually remove the Crosslinking SC terms, where Gly has no SC
@author: chris Liguo
"""
import numpy as np
import vermouth
from vermouth.forcefield import ForceField
import os

def interaction_adding(ff, molname, first_BB, first_SC, second_BB, second_SC):

    '''
    join two components with interactions
    '''
    # link the last and first BB of two sections
    ff.blocks[molname].add_interaction('bonds',
                                       atoms = [first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['1', '0.360', '8000'],
                                       meta={'comment':'BB linking'}
                                       )

    # BBB angles
    ff.blocks[molname].add_interaction('angles',
                                       atoms = [first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['10', '137', '25'],
                                       meta={'comment':'BBB linking'}
                                       )
    ff.blocks[molname].add_interaction('angles',
                                       atoms = [first_BB[-1],
                                                second_BB[0],
                                                second_BB[1]],
                                       parameters = ['10', '137', '25'],
                                       meta={'comment':'BBB linking'}
                                       )
    #SBB and BBS angles
    ff.blocks[molname].add_interaction('angles',
                                       atoms = [first_SC[-1],
                                                first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'SBB linking'}
                                       )
    ff.blocks[molname].add_interaction('angles',
                                       atoms = [first_BB[-1],
                                                second_BB[0],
                                                second_SC[0]],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'BBS linking'}
                                       )
    #BBBB dihedrals
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-3],
                                                first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['9', '60', '2.8', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-3],
                                                first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['9', '150', '-0.60', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-3],
                                                first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0]],
                                       parameters = ['9', '130', '-1.20', '2'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )                                       
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0],
                                                second_BB[1]],
                                       parameters = ['9', '60', '2.8', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0],
                                                second_BB[1]],
                                       parameters = ['9', '150', '-0.60', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-2],
                                                first_BB[-1],
                                                second_BB[0],
                                                second_BB[1]],
                                       parameters = ['9', '130', '-1.20', '2'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )                                       
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-1],
                                                second_BB[0],
                                                second_BB[1],
                                                second_BB[2]],
                                       parameters = ['9', '60', '2.8', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-1],
                                                second_BB[0],
                                                second_BB[1],
                                                second_BB[2]],
                                       parameters = ['9', '150', '-0.60', '1'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_BB[-1],
                                                second_BB[0],
                                                second_BB[1],
                                                second_BB[2]],
                                    parameters = ['9', '130', '-1.20', '2'],
                                       meta={'comment':'BBBB linking using generic definition in IDP'}
                                       )                                      
    #SBBS dihedral
    #used in the case where the 2 linking residues don't contain Gly.
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_SC[-1],
                                                first_BB[-1],
                                                second_BB[0],
                                                second_SC[0]],
                                       parameters = ['9', '-130', '-1.4', '1'],
                                       meta={'comment':'SBBS linking using definition in IDP'}
                                       )
    ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [first_SC[-1],
                                                first_BB[-1],
                                                second_BB[0],
                                                second_SC[0]],
                                       parameters = ['9',  '100', '-1.4', '2'],
                                       meta={'comment':'SBBS linking using definition in IDP'}
                                       )
    #give up BB-BB(-1)-BB(+1)-SC improper dihedral in Martini-IDP

    return ff
def entry_edit(mols, pre_mol, adding_mol):
    '''
    based on the last index of pre_mol, give directives in adding_mol new indices
    '''
    nres_prev = mols[pre_mol]['atoms'][-2][2]
    natom_prev= mols[pre_mol]['atoms'][-2][0]
    
    for i in mols[adding_mol]['atoms']:
        if len(i)>1:
            i[0] = i[5] = str(int(i[0]) + int(natom_prev))     
            i[2] = str(int(i[2]) + int(nres_prev))
        
    for i in mols[adding_mol]['position_restraints']:
        if len(i)>1:
            try:
                i[0] = str(int(i[0]) + int(natom_prev))     
            except ValueError:
                pass
    for i in mols[adding_mol]['bonds']:
        if len(i)>1:
            try:
                i[0] = str(int(i[0]) + int(natom_prev))     
                i[1] = str(int(i[1]) + int(natom_prev))     
            except ValueError:
                pass
    for i in mols[adding_mol]['constraints']:
        if len(i)>1:
            try:
                i[0] = str(int(i[0]) + int(natom_prev))     
                i[1] = str(int(i[1]) + int(natom_prev))     
            except ValueError:
                pass
    for i in mols[adding_mol]['angles']:
        if len(i)>1:
            i[0] = str(int(i[0]) + int(natom_prev))     
            i[1] = str(int(i[1]) + int(natom_prev))     
            i[2] = str(int(i[2]) + int(natom_prev))     
 
    for i in mols[adding_mol]['dihedrals']:
        if len(i)>1:
            i[0] = str(int(i[0]) + int(natom_prev))     
            i[1] = str(int(i[1]) + int(natom_prev))     
            i[2] = str(int(i[2]) + int(natom_prev))     
            i[3] = str(int(i[3]) + int(natom_prev)) 

#    if mols[adding_mol].get('dihedrals1'):
#        for i in mols[adding_mol]['dihedrals1']:
#            if len(i)>1:
#                i[0] = str(int(i[0]) + int(natom_prev))     
#                i[1] = str(int(i[1]) + int(natom_prev))     
#                i[2] = str(int(i[2]) + int(natom_prev))     
#                i[3] = str(int(i[3]) + int(natom_prev))     
            
    if mols[adding_mol].get('virtual_sitesn'):
        l = []
        for i in mols[adding_mol]['virtual_sitesn']:
            if len(i)>1:
                new = []        
                new.append(str(int(i[0]) + int(natom_prev)))
                new.append(i[1])
                for j,k in enumerate(i[2:]):
                    new.append(str(int(k) + int(natom_prev)))
                l.append(new)
        mols[adding_mol]['virtual_sitesn'] = l

    if mols[adding_mol].get('exclusions'):
        l = []
        for i in mols[adding_mol]['exclusions']:
            if len(i)>1:
                new = []
                for j,k in enumerate(i):
                    new.append(str(int(k) + int(natom_prev)))
                l.append(new)
        mols[adding_mol]['exclusions'] = l
    
    pre_BB = np.array([int(i[0])-1 for i in [i for i in mols[pre_mol]['atoms'] if 'BB' in i]], dtype = int)
    pre_SC = np.array([int(i[0])-1 for i in [i for i in mols[pre_mol]['atoms'] if 'SC1' in i]], dtype = int)
    added_BB = np.array([int(i[0])-1 for i in [i for i in mols[adding_mol]['atoms'] if 'BB' in i]], dtype = int)
    added_SC = np.array([int(i[0])-1 for i in [i for i in mols[adding_mol]['atoms'] if 'SC1' in i]], dtype = int)
    return pre_BB, pre_SC, added_BB, added_SC

if __name__ == '__main__':

    basedir = os.getcwd()+'/'
    #list of the itps that you want to merge
    fs = ['mol0.itp', 'mol1.itp']
    #names of the molecules to merge
    na = ['mol0', 'mol1']
    
    mols = {}
    for i,j in zip(fs,na):
        
        with open(basedir+i) as f:
            lines = f.readlines()
        
        mol_data = {}
        d = []
        first = True
        for line in lines:
            if '[' in line:
                if len(d) > 1 and first == False:
                    try:
                        mol_data[directive] = d
                    except NameError:
                        pass
                d = []
                directive = line.split('[')[1].split(']')[0].strip()
                #this needs to be done because martinize will split up the proper and improper dihedrals
                if directive in list(mol_data.keys()):
                    directive = directive+'1'
                first = False
                # print(directive)
        
            elif line[0] != ';':# and len(line.split()) > 1:
                d.append(line.split())
    
        #this makes sure the final directive read in is added to the current molecule
        mol_data[directive] = d        
        #and then we add the molecule to the overall dictionary
        mols[j] = mol_data
    
    #this is very hacky but assume mol0 is folded and has 2 [ dihedrals ] in the original itp. So merge the dict entry.
    #actually mol2 (if assumed as IDR), it also has 2 [dihedrals] in the original itp, so merge the dict entry, too.
    mols['mol0']['dihedrals'] = mols['mol0']['dihedrals'] + mols['mol0']['dihedrals1']
    mols['mol1']['dihedrals'] = mols['mol1']['dihedrals'] + mols['mol1']['dihedrals1']

    #this edits the entries that you have to ensure correct indexing.
    #write them in the order that the function says, the bookkeeping is up to you!
    mol0_BB, mol0_SC, mol1_BB, mol1_SC = entry_edit(mols, 'mol0', 'mol1')
    #mol1_BB, mol1_SC, mol2_BB, mol2_SC = entry_edit(mols, 'mol1', 'mol2')
    #etc
    
            
    keys = ['atoms', 'position_restraints',
            'bonds', 'constraints', 
            'angles', 'dihedrals',
            'virtual_sitesn', 'exclusions']
    final_mol = {}
    
    
    for i in keys:
        entry = []
        for j in mols.keys():
            if mols[j].get(i):
                entry+= mols[j][i]
        final_mol[i] = entry
               
            
    test = []
    test.append('[ moleculetype ]\n')
    test.append('molecule 1\n\n')
    
    for key in final_mol.keys():
    
        test.append(f'[ {key} ]\n')
        for line in final_mol[key]:
            test.append('\t'.join(line)+'\n')
            
    
        
    ff = ForceField('martini3001')
    
    mol = vermouth.gmx.read_itp(test, ff)
    
    molnames = list(ff.blocks.keys())
    molname = molnames[0]
    ff.blocks[molname].meta['moltype'] = molname
    
    '''
    now add the linking interactions between the different components
    can do this for however many components you want as long as the bookkeeping is right.
    '''
    
    interaction_adding(ff, molname, mol0_BB, mol0_SC, mol1_BB, mol1_SC)

    
    mol_out = ff.blocks[molname].to_molecule()    
    mol_out.meta['moltype'] = molname
    
    header_lines = ['itp stitched by Chris Liguo\n;Manually change the protonation state in crosslinking']
    print('Remember manually change the protonation state in crosslinking !!!!!!!!')
    print('Remember manually remove the Crosslinking SC terms, where Gly has no SC!!!!!!!!!!!!')
    
    with open(basedir + 'mol.itp', 'w') as outfile:
        vermouth.gmx.write_molecule_itp(mol_out, outfile=outfile, header=header_lines)

        
        
        
        
        
