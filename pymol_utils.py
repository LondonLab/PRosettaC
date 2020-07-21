import os
import subprocess
import sys
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored, time
import random
harmless_hetatm = ['PO4', 'SO4', 'K', 'NA', 'CA', 'OXM', 'OXL', 'MN', 'MAE', 'CO', 'FE']

def publication_figure(session, pic_file_name, w = 5000, h = 5000):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(session)
    pnghack(pic_file_name, w, h)

def align_chain(ref, query, chain, out):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(ref)
    cmd.load(query)
    ref_name = ref.split('/')[-1].split('.')[0]
    query_name = query.split('/')[-1].split('.')[0]
    cmd.align(query_name + ' and chain ' + chain, ref_name + ' and chain ' + chain)
    cmd.save(out, query_name)

def get_distance(file_name, res1, atom1, res2, atom2):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('atom1', 'resn ' + res1 + ' and name ' + atom1)
    cmd.select('atom2', 'resn ' + res2 + ' and name ' + atom2)
    print(cmd.get_distance('atom1', 'atom2'))

def center_coords(file_name):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    xyz_limits = cmd.get_extent()
    xyz_mean = [(xyz_limits[0][i] + xyz_limits[1][i]) / 2 for i in range(3)]
    cmd.translate([-1 * xyz_mean[i] for i in range(3)], 'all')
    cmd.save(file_name, 'all')

def center_coords_rec(rec):
    xyz_limits = cmd.get_extent(rec)
    xyz_mean = [(xyz_limits[0][i] + xyz_limits[1][i]) / 2 for i in range(3)]
    cmd.translate([-1 * xyz_mean[i] for i in range(3)], 'all')

def delete_hetero(lig):
    cmd.select('env_het', lig + ' around 5 and hetatm')
    cmd.remove('env_het')

def get_chain(lig):
    cmd.select('env', lig + ' around 5 and (not hetatm)')
    chains = []
    cmd.iterate('env', 'chains.append(chain)', space = locals())
    return chains[0]

def get_seq(file_name, chain):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('seq_chain', 'chain ' + chain)
    return cmd.get_fastastr('seq_chain')

def num_interacting_chains(file_name, mol, chain, res):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain + ' and resi ' + res)
    cmd.select('env', 'lig around 5 and (not hetatm)')
    chains = []
    cmd.iterate('env', 'chains.append(chain)', space = locals())
    return ''.join(list(set(chains)))

def is_pure_env(file_name, mol, chain, res):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain + ' and resi ' + res + ' and alt A+ and not hydrogen')
    if cmd.count_atoms('lig') < 10 or cmd.count_atoms('lig') > 33:
        return False
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain + ' and resi ' + res)
    line = ''
    for hetatm in harmless_hetatm:
        line += ' and (not resn ' + hetatm + ')'
    cmd.select('env', 'lig around 5 and hetatm and (not resn HOH)' + line)
    if cmd.count_atoms('env') == 0:
        return True
    else:
        return False

def env_cysteine(file_name, allowed_het):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('het', 'org')
    stored.list=[]
    cmd.iterate('org', "stored.list.append((resn))")
    het = list(set(stored.list))
    return_table = []
    #with open(het_list, 'r') as f:
    #    het = [l.split()[0] for l in f]
    with open(allowed_het, 'r') as f:
        allowed = [l.split()[0] for l in f]
    for h in het:
        if not h in allowed:
            continue
        cmd.select('lig', 'resn ' + h + ' and not elem H')
        cmd.select('env', 'lig around 6 and resn CYS and elem S')
        stored.list=[]
        cmd.iterate('env', "stored.list.append((resi, chain))")
        res_chain = set(stored.list)
        new_res_chain = []
        for c in res_chain:
            cmd.select('sul', 'resi ' + c[0] + ' and chain ' + c[1] + ' and elem S')
            cmd.select('neighbors', 'not name CB and neighbor sul')
            if cmd.count_atoms('neighbors') == 0:
                new_res_chain.append(c)
        res_chain = new_res_chain
        if len(res_chain) > 0:
            resi = set([c[0] for c in res_chain])
            for r in resi:
                entry = [c for c in res_chain if c[0] == r][0]
                if is_lig_single(h, entry[1]):
                    print(h + '\t' + entry[0] + '\t' + entry[1])
                    return_table.append([h, entry[0], entry[1]])
    return return_table

def get_rec_plus_lig(pdb_id, lig, rec_file, lig_file, new_chain):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.fetch(pdb_id)
    center_coords_rec(pdb_id)
    cmd.select('lig', 'resn lig')
    stored.list=[]
    cmd.iterate('org', "stored.list.append((chain))")
    chain = list(set(stored.list))[0][0]
    cmd.select('lig', 'resn ' + lig + ' and chain ' + chain)
    cmd.select('rec', 'poly and chain ' + chain)
    cmd.alter('rec', 'chain=\'' + new_chain + '\'')
    cmd.save(lig_file, 'lig')
    cmd.save(rec_file, 'rec')
    os.system('rm *.cif')

def save_residue(file_name, resi, save_file):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('res', 'resi ' + str(resi))
    cmd.save(save_file, 'res')

def seperate_rec_res(file_name, resi, resn):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + resn + ' and resi ' + resi)
    if cmd.count_atoms('lig') == 0:
        return 0
    cmd.select('rec', 'not hetatm')
    center_coords_rec('rec')
    cmd.save('rec.pdb', 'rec')
    cmd.save('xtal-lig.pdb', 'lig')

def seperate_rec_lig(file_name, mol, chain, res):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain + ' and resi ' + res)
    if cmd.count_atoms('lig') == 0:
        return 0
    chain = get_chain('lig')
    cmd.select('rec', 'not lig and not hetatm and chain ' + chain)
    center_coords_rec('rec')
    cmd.save('rec.pdb', 'rec')
    cmd.save('xtal-lig.pdb', 'lig')
    cmd.remove(file_name[:-4])

def seperate_rec_lig(file_name, mol, chain):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain)
    if cmd.count_atoms('lig') == 0:
        return 0
    stored.list=[]
    cmd.iterate('lig', "stored.list.append((resi))")
    num_of_lig = len(set(stored.list))
    if not num_of_lig == 1:
        return 0
    chain = get_chain('lig')
    cmd.select('rec', 'not lig and not hetatm and chain ' + chain)
    center_coords_rec('rec')
    cmd.save('rec.pdb', 'rec')
    cmd.save('xtal-lig.pdb', 'lig')

def pnghack(filepath, width=1024, height=768):
    #Workaround if cmd.png() doesn't work
    cmd.set('ray_trace_frames', 1)  # Frames are raytraced before saving an image.
    cmd.set('ray_shadows', 0)
    cmd.viewport(width, height)  # Set resolution
    cmd.mpng(filepath, 1, 1)  # Use batch png mode with 1 frame only
    cmd.mplay()  # cmd.mpng needs the animation to 'run'

def scene_photo(rec, xtal, covalentized, photo_name):
    pymol.finish_launching()
    cmd.delete('all')
    for file_name in [rec, xtal, covalentized]:
        cmd.load(file_name)
    cmd.select('lig', 'org')
    if cmd.count_atoms('lig') == 0:
        return 0
    cmd.set('valence', 0)
    cmd.color('cyan', xtal.split('.')[0])
    cmd.select('cov', covalentized.split('.')[0])
    cmd.select('cysteine', 'br. ' + covalentized.split('.')[0] + ' around 2 and resn CYS')
    cmd.select('cys_cov', 'cysteine or cov')
    cmd.save('tmp.pdb', 'cys_cov')
    cmd.load('tmp.pdb')
    cmd.hide('(hydro)')
    cmd.delete('cov')
    cmd.select('cov', 'tmp and org')
    cmd.select('cysteine', 'tmp and resn cys')
    cmd.color('white', 'cysteine')
    cmd.color('magenta', 'cov')
    cmd.color('green', rec.split('.')[0])
    cmd.util.cnc('all')
    cmd.select('all_mols', 'tmp or ' + xtal.split('.')[0])
    cmd.orient('all_mols')
    pnghack('./tmp.png')
    os.rename('tmp0001.png', photo_name)
    os.remove('tmp.pdb')

def is_lig_single(file_name, mol, chain):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain)
    if cmd.count_atoms('lig') == 0:
        return False
    stored.list=[]
    cmd.iterate('lig', "stored.list.append((resi))")
    num_of_lig = len(set(stored.list))
    if not num_of_lig == 1:
        return False
    return True

def is_lig_single(mol, chain):
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain)
    if cmd.count_atoms('lig') == 0:
        return False
    stored.list=[]
    cmd.iterate('lig', "stored.list.append((resi))")
    num_of_lig = len(set(stored.list))
    if not num_of_lig == 1:
        return False
    return True

def seperate_lig(file_name, mol, chain, output_name):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.select('lig', 'resn ' + mol + ' and chain ' + chain)
    cmd.save(output_name, 'lig')

def get_group_env(file_name, group, output_name):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)
    cmd.load(group)
    cmd.select('group_env', 'br. ' + group.split('.')[0] + ' around 4 and poly')
    cmd.save(output_name, 'group_env')

def refresh():
    cmd.sync()
    cmd.refresh()
    time.sleep(0.1)

def get_image(ind, pdb_name, mol, chain, res, original = True):
    data_folder = 'Data/'
    father_folder = data_folder + 'PNG_' + str(original) + '/'
    #if not os.path.exists(father_folder):
    #    os.mkdir(father_folder)
    folder = father_folder + ind + '_' + pdb_name + '/'
    if not os.path.exists(folder):
        os.mkdir(folder)
    ugly_id = 'cand4'
    true_id = 'cand1'
    docking_folder = 'RosettaDock/'
    file_name = docking_folder + ind + '_' + pdb_name + '/rec1.pdb'
    ugly_rec = docking_folder + ind + '_' + pdb_name + '/rec4.pdb'
    ugly_name = docking_folder + ind + '_' + pdb_name + '/' + ugly_id + '.mol2'
    xtal = docking_folder + ind + '_' + pdb_name + '/' + true_id + '.mol2'

    steps = 3
    step_angle = 360.0 / steps
    pixels = 100
    num_rand = 9

    #This is done for negative of random rotation in the pocket
    '''if not original:
        fake_rotation = []
        for _ in range(3):
            fake_rotation.append(random.random() * 360)'''

    pymol.finish_launching()
    cmd.delete('all')
    if original:
        cmd.load(file_name, 'protein')
    else:
        if not os.path.exists(ugly_rec):
            os.rmdir(folder)
            return 0
        cmd.load(ugly_rec, 'protein')
    cmd.load(xtal)
    #cmd.select('lig', 'resn ' + mol + ' and chain ' + chain + ' and resi ' + res)
    cmd.select('lig', true_id)
    obj_name = true_id
    if not original:
        if not os.path.exists(ugly_name):
            os.rmdir(folder)
            return 0
        cmd.load(ugly_name)
        cmd.select('lig', ugly_id)
        obj_name = ugly_id
    delete_hetero('lig')
    cmd.remove('hydrogen')
    if cmd.count_atoms('lig') == 0:
        os.rmdir(folder)
        return 0
    #cmd.select('rec', 'not lig')
    cmd.select('rec', 'br. lig around 4 and poly')
    cmd.hide("all")
    cmd.color('white', 'rec')
    cmd.color('green', 'lig')
    cmd.util.cnc("all")
    cmd.show("sticks", 'lig')
    #Show H bonds
    cmd.set('h_bond_cutoff_center', 3.5)
    cmd.distance('hbonds', 'lig', 'rec', mode=2)
    cmd.hide('labels')
    #Show double bonds
    #cmd.set('valence', 1, obj_name)
    cmd.orient('lig')
    cmd.zoom('lig', 3)
    #Slab move -5
    cmd.clip('move', -5)
    cmd.set('ray_opaque_background', 0)                                                                         
    cmd.set('ray_shadows', 0)

    #This is done for negative of random rotation in the pocket
    '''if not original:
        cmd.rotate('x', fake_rotation[0], 'lig')
        cmd.rotate('y', fake_rotation[1], 'lig')
        cmd.rotate('z', fake_rotation[2], 'lig')'''
    
    #Transperent surface
    '''cmd.set('surface_carve_cutoff', 4.5)
    cmd.set('surface_carve_selection', 'lig')
    cmd.set('surface_carve_normal_cutoff', -0.1)
    cmd.set('surface_color', 'white')
    cmd.set('surface_type', 3)
    cmd.set('transparency', 0.5)'''

    for show in ['sticks']:#['surface', 'sticks']:
        #Set transparency
#        tra = 0
#        if show == 'surface':
#            tra = 0.5
#        cmd.set('transparency', tra)
        cmd.show(show, 'rec')
        refresh()
        for x in range(steps):
            for y in range(steps):
                refresh()
                cmd.png(folder + show + '_' + str(x) + '_' + str(y) + '.png', pixels, pixels)
                refresh()
                cmd.turn('y', step_angle)
                refresh()
            cmd.turn('x', step_angle)
            refresh()
        '''for x in range(num_rand):
            refresh()
            for ax in ['x', 'y', 'z']:
                cmd.turn(ax, random.random() * 360)
            refresh()
            cmd.png(folder + show + '_' + str(x) + '.png', pixels, pixels)
            refresh()'''
        cmd.hide(show, 'rec')
        refresh()
    return 1

def get_dude_image(entry, mol, original = True):
    data_folder = 'Data/'
    father_folder = data_folder + 'PNG_' + str(original) + '/'
    if not os.path.exists(father_folder):
        os.mkdir(father_folder)
    folder = father_folder + entry + '/'
    if not os.path.exists(folder):
        os.mkdir(folder)
    mol_name = mol[:-5]
    folder += mol_name + '/'
    if not os.path.exists(folder):
        os.mkdir(folder)
    docking_folder = entry
    if original:
        docking_folder += '/actives_Dock/poses/'
    else:
        docking_folder += '/decoys_Dock2/poses/'
    file_name = entry + '/rec.pdb'
    xtal = docking_folder + mol

    steps = 3
    step_angle = 360.0 / steps
    pixels = 100
    num_rand = 9

    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name, 'protein')
    cmd.load(xtal, mol_name)
    cmd.select('lig', mol_name)
    obj_name = mol_name
    delete_hetero('lig')
    cmd.remove('hydrogen')
    if cmd.count_atoms('lig') == 0:
        os.rmdir(folder)
        return 0
    cmd.select('rec', 'br. lig around 4 and poly')
    cmd.hide("all")
    cmd.color('white', 'rec')
    cmd.color('green', 'lig')
    cmd.util.cnc("all")
    cmd.show("sticks", 'lig')
    #Show H bonds
    cmd.set('h_bond_cutoff_center', 3.5)
    cmd.distance('hbonds', 'lig', 'rec', mode=2)
    cmd.hide('labels')
    #Show double bonds
    #cmd.set('valence', 1, obj_name)
    cmd.orient('lig')
    cmd.zoom('lig', 3)
    #Slab move -5
    cmd.clip('move', -5)
    cmd.set('ray_opaque_background', 0)                                                                         
    cmd.set('ray_shadows', 0)

    for show in ['sticks']:#['surface', 'sticks']:
        #Set transparency
#        tra = 0
#        if show == 'surface':
#            tra = 0.5
#        cmd.set('transparency', tra)
        cmd.show(show, 'rec')
        refresh()
        for x in range(steps):
            for y in range(steps):
                refresh()
                cmd.png(folder + show + '_' + str(x) + '_' + str(y) + '.png', pixels, pixels)
                refresh()
                cmd.turn('y', step_angle)
                refresh()
            cmd.turn('x', step_angle)
            refresh()
        cmd.hide(show, 'rec')
        refresh()
    return 1

def get_surface_area(file_name):    
    pymol.finish_launching()
    cmd.load(file_name)
    area = cmd.get_area('resi 1')
    cmd.remove(file_name[:-5])
    return area

def pymol_mutate(file_name, chain, res_index, number, mutant):
    pymol.finish_launching()
    cmd.delete(file_name[:-4])
    selection = chain + '/' + res_index + '/'
    cmd.wizard("mutagenesis")
    pdb = file_name[:-4]
    cmd.load(file_name)
    cmd.refresh_wizard()
    cmd.get_wizard().set_mode(mutant)
    cmd.get_wizard().do_select(selection)
    nStates = cmd.count_states("mutation")
    for i in range(1, nStates + 1):
        cmd.get_wizard().do_select(selection)
        cmd.frame(i)
        cmd.get_wizard().apply()
        cmd.save("rec_" + str(res_index) + "_" + str(i) + ".pdb")

    cmd.set_wizard()
    cmd.remove(file_name[:-4])

def pymol_mutate(file_name, chain, res_index):
    pymol.finish_launching()
    cmd.delete('all')
    selection = chain + '/' + res_index + '/'
    mutant = 'CYS'
    cmd.wizard("mutagenesis")
    pdb = file_name[:-4]
    cmd.load(file_name)
    cmd.remove('not (alt ''+A)')

    cmd.select('mut', 'resi ' + res_index + ' and chain ' + chain)
    if cmd.count_atoms('mut') == 0:
        return False

    cmd.refresh_wizard()
    cmd.get_wizard().set_mode(mutant)
    cmd.get_wizard().do_select(selection)
    nStates = cmd.count_states("mutation")
    for i in range(1, nStates + 1):
        cmd.get_wizard().do_select(selection)
        cmd.frame(i)
        cmd.get_wizard().apply()
        cmd.save("rec_" + str(res_index) + "_" + str(i) + ".pdb")

    cmd.set_wizard()
    cmd.remove(file_name[:-4])
    return True

def show_bumps(selection):
    print(selection)
    name = 'bump_check'
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        print('VDW Strain in state %d: %f' % (state, strain))

def findSurfaceResidues(file_name, objSel="(all)", cutoff=2.5, doShow=False, verbose=False, only_cysteine=False):
    """
    findSurfaceResidues
finds those residues on the surface of a protein
that have at least 'cutoff' exposed A**2 surface area.

    PARAMS
objSel (string)
    the object or selection in which to find
    exposed residues
    DEFAULT: (all)

cutoff (float)
    your cutoff of what is exposed or not. 
    DEFAULT: 2.5 Ang**2

asSel (boolean)
    make a selection out of the residues found

    RETURNS
(list: (chain, resv ) )
    A Python list of residue numbers corresponding
    to those residues w/more exposure than the cutoff.

    """
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)

    tmpObj="__tmp"
    #if only_cysteine:
    #    cmd.create( tmpObj, objSel + " and polymer and resn CYS");
    #else:
    cmd.create( tmpObj, objSel + " and polymer");
    if verbose!=False:
        print("WARNING: I'm setting dot_solvent.  You may not care for this.")
    cmd.set("dot_solvent");
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    if only_cysteine:
        cmd.remove( tmpObj + " and not resn CYS")
        cmd.remove( tmpObj + " and not elem S")
        cmd.remove( tmpObj + " and CYS/SG and bound_to CYS/SG")
    cmd.remove( tmpObj + " and b < " + str(cutoff) )
    cmd.iterate(tmpObj, "b")

    stored.tmp_dict = {}
    if only_cysteine:
        cmd.iterate(tmpObj + " and resn CYS", "stored.tmp_dict[(chain,resv)]=1")
    else:
        cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
    exposed = stored.tmp_dict.keys()
    exposed.sort()

    randstr = str(random.randint(0,10000))
    selName = "exposed_atm_" + randstr
    if verbose!=False:
        print("Exposed residues are selected in: " + selName)
    cmd.select(selName, objSel + " in " + tmpObj ) 
    selNameRes = "exposed_res_" + randstr
    cmd.select(selNameRes, "byres " + selName )
    cmd.delete(tmpObj)

    exposed = [i[1] for i in exposed]
    return exposed

def solventExposure(file_name, resi):
    pymol.finish_launching()
    cmd.delete('all')
    cmd.load(file_name)

    tmpObj="__tmp"
    cmd.create( tmpObj, "(all) and polymer");
    cmd.set("dot_solvent");
    cmd.get_area(selection=tmpObj, load_b=1)
    stored.list=[]
    cmd.remove( tmpObj + " and not resn CYS")
    cmd.remove( tmpObj + " and not elem S")
    cmd.remove( tmpObj + " and CYS/SG and bound_to CYS/SG")
    cmd.remove( tmpObj + " and not resi " + resi)
    cmd.iterate(tmpObj, "stored.list.append((b))")
    #return sum(stored.list) / len(stored.list)
    return max(stored.list)
