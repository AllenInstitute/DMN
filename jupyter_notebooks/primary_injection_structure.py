# -*- coding: utf-8 -*-

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.ontologies_api import OntologiesApi

mcc = MouseConnectivityCache(manifest_file='../connectivity/mouse_connectivity_manifest.json', resolution=25)
structure_tree = mcc.get_structure_tree()

def get_primary_injection_structure(experiment_id):
    '''
    Parameters
    ----------
    experiment ids : integer or list of integer experiment ids

    Returns
    -------
    structure : string or array of strings
        Acronym for injection structure with highest projection volume in each experiment
    '''

    oapi = OntologiesApi()
    structures = oapi.get_structures(structure_set_names="'Mouse Connectivity - Summary'")
    ss_ids = [item['id'] for item in structures]
    structure_dat = {item['id']: item['acronym'] for item in structures}

    unionizes = mcc.get_structure_unionizes(experiment_ids=[experiment_id], is_injection=True)
    valid_structures = []
    for structure in unionizes['structure_id']:
        if structure in ss_ids:
            valid_structure = structure
        else:
            valid_structure = find_ss_parent(structure, ss_ids)
        valid_structures.append(valid_structure)
    if len(unionizes['structure_id']) == len(valid_structures):
        unionizes['structure_id'] = valid_structures
    else:
        print('length mismatch')
    inj_site_by_ss = unionizes.groupby('structure_id')['projection_volume'].mean()
    inj_site_max = inj_site_by_ss.keys()[inj_site_by_ss.values == max(inj_site_by_ss)].values[0]
    inj_site_second = inj_site_by_ss.delete(inj_site_max)
    primary_injection_acronym = structure_dat[int(inj_site_max)]
    secondary_injection_acronym = structure_dat[int(inj_site_second)]
    return primary_injection_acronym, secondary_injection_acronym


def find_ss_parent(structure_id, summary_structure_ids):
    '''
    Parameters
    ----------
    structure_id : integer structure id to start the search
    summary_structure_ids: list of structure ids to search within for valid results

    Returns
    -------
    parent_id: integer code for next ancestor in summary structures
    '''
    parent_id = structure_tree.parent_id([structure_id])[0]
    if parent_id == None:
        return parent_id
    else:
        while parent_id not in summary_structure_ids:
            parent_id = structure_tree.parent_id([parent_id])[0]
            if parent_id == None:
                return parent_id
        return parent_id
