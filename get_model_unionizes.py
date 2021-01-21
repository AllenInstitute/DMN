# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 10:52:29 2019

@author: jenniferwh
"""

import os
from mcmodels.core import VoxelModelCache
from anatomy.anatomy_api import AnatomyApi
import pandas as pd
import numpy as np
import scipy.stats as st
import platform
import json

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi

def main():
    cache = VoxelModelCache(manifest_file='connectivity/voxel_model_manifest.json', resolution = 100)
    mcc = MouseConnectivityCache(manifest_file = 'connectivity/mouse_connectivity_manifest.json', resolution = 100)
    mca = MouseConnectivityApi()
    aapi = AnatomyApi()
    
    ss = aapi.get_summary_structure_data('id')
    
    td_experiments = pd.DataFrame(mcc.get_experiments(cre=['Ai75(RCL-nt)']))
    print(len(td_experiments))
    
    voxel_array, source_mask, target_mask = cache.get_voxel_connectivity_array()
    
    isids = [500836840,584194481,523718823,478376197,500836105,
             294481346,512315551,657046319,501483880,166153483,598605738,
             606929366,593018150,496554237,522409371,521600943,501484658,
             526783054,479981981,647806688,502966396,518742338,182803137,
             479673174,512314723,266249483,500837552,166055636,125833030,
             272928602,480994108,249396394,606100558,559878074,120916102,
             574941472,294482052,287769286,501006221,478491090,530574594,
             530555292,482578964,579203888,576332845,540145406,572588941,
             648253235,272735030,266250195,510834706,286313491,495562600,
             583747816,511550172,503809372,584895127,182793477,482581199,
             648252402,657041814,496113850,126117554,646525997,576341623,
             536298726,479755622,518606617,183330908,571100135,182090318,
             510417993,479700629,590988082,479701339,479984127,518605900,
             583748537,575782182,272698650,524667618,657172447,596253012,
             167794131,592540591,298179622,478376911,166083557,479982715,
             554334685,657334568,168164230,502074651,286482701,272929308,
             120875816,166054929,566992832,166082842,535692871,563179067,
             286314623,496151571,501115762,603479758,156741826,536299435,
             657042668,479983421,501004474,562061175,159753308,562521707,
             159832064,546103149,576340889,656959718,286417464,562671482,
             297652799,168163498,517072832,121510421,562520963,593016678,
             183171679,478582494,120814821,303785454,646525156,479673887,
             503323656,298178204,286312782,161458737,168162771,524874308,
             648050335,597007143,606785720,272930013,482580380,496576666,
             661987195,584511827,479980810,510417255,593277684,263106036,
             485237081,480069939,591535205,303784745,249327301,562674923,
             167902586,524666904,523718075,166082128,263242463,601476074,
             183173527,656632388,580150240,598604288,642809043,554421791,
             479756361,296048512,503069254] #Emx1 and Rbp4
    '''
    isids = [475829896, 477435412, 528741104, 567301515, 571401645, 607059419,
       607321130, 609475867, 637855050, 649362978, 475733382, 528511254,
       569932566, 607052300, 623286726, 475616836, 592698087, 592698832,
       601900484, 606278526, 605092364, 614435699, 478995566, 
       528511967, 483094671, 522274859, 523704892, 525413115, 528328648,
       502005076, 527390805, 527393013, 528512680, 530000865, 531233132,
       478877019, 521954271, 527713163, 527713883, 529435133, 524266253,
       578332611, 609475139, 569993539, 571652998, 572753855, 573977678,
       657424207, 571816813, 572390577, 601804603, 636799953, 475617622,
       524267323, 479115470, 520615681, 521955016, 523177830, 526783792,
       532005897, 617898760, 617900105, 523481517, 575683857, 606260719,
       609157409, 614738437, 496965687, 496964969, 528732005, 607289053,
       607316031, 572595932, 573331541, 595259180, 604100536, 617901499,
       518605181, 529129011, 536920234, 475828414, 475830603, 521255975,
       524063357, 525412369, 531443949, 555012592, 561307215, 561511939,
       561512675, 561910766, 569904687, 570071403, 571647261, 592522663,
       592724077, 605496542, 616674416, 623838656, 664716091, 666090944,
       570461017, 601904029, 521617657, 569623967, 591394174,
       544454674, 478678606, 560045081, 561986735, 591168591, 613898292,
       651703553, 501785691, 501883865, 502590301, 502955689, 502956560,
       504176074, 518012479, 518013943, 523180728, 531397136, 552543088,
       553080579, 560965104, 571653937, 572388976, 560029094,
       485553574, 495346667, 501711996, 501786400, 501787135,
       501837158, 502592260, 515920693, 526784559, 539323512,
       539511058, 539514325, 543297033, 546389260, 561918178, 561918904,
       563205064, 563352720, 565146821, 574637452, 568502684, 571410278,
       561506791, 605093074, 601885751, 572771153] #target-defined experiments
    
    isids = [484612961, 529428776] #failed
    '''
    
    for isid in isids:
        print(isid)
        flip = False
        data_mask = mcc.get_data_mask(experiment_id=isid)[0]
        density = mcc.get_injection_density(isid)[0]
        injection_density = np.multiply(density, data_mask)
        injection_fraction = mcc.get_injection_fraction(isid)[0]
        inj_unionize = mcc.get_structure_unionizes(experiment_ids = [isid],
                                                    is_injection = True,
                                                    structure_ids = ss)
        r_hem = inj_unionize[inj_unionize['hemisphere_id'] == 2]['projection_volume'].sum()
        l_hem = inj_unionize[inj_unionize['hemisphere_id'] == 1]['projection_volume'].sum()
        if l_hem > r_hem:
            flip = True
            injection_density = np.flip(injection_density, 2) 
            injection_fraction = np.flip(injection_fraction, 2)
        centroid = mca.calculate_injection_centroid(
                injection_density, injection_fraction, resolution=1)
        centroid = tuple(map(int, centroid))
        valid_voxels = np.where(injection_density > 0)
        valid_voxels = np.swapaxes(valid_voxels, 0, 1)
        
        try:
            row = source_mask.get_flattened_voxel_index(centroid)
            volume = np.zeros_like(voxel_array[row])
            for voxel in valid_voxels:
                try:
                    row = source_mask.get_flattened_voxel_index(voxel)
                    new_volume = voxel_array[row] * injection_density[tuple(voxel)]
                    volume = np.add(volume, new_volume)
                except:
                    print(voxel)
            
            model = target_mask.map_masked_to_annotation(volume)
            
            experimental_unionize = mcc.get_experiment_structure_unionizes(
                    experiment_id = isid, structure_ids = ss, hemisphere_ids = [1,2])
            exp_unionize = experimental_unionize[['structure_id', 'hemisphere_id', 'normalized_projection_volume']]
            exp_unionize.loc[(np.log10(exp_unionize['normalized_projection_volume']) < -2.5),
                                      'normalized_projection_volume'] = 0
            exp_unionize['data'] = 'experimental'
            
            if flip:
                exp_unionize.loc[exp_unionize['hemisphere_id'] == 1, 'hemisphere_id'] = 0
                exp_unionize.loc[exp_unionize['hemisphere_id'] == 2, 'hemisphere_id'] = 1
                exp_unionize.loc[exp_unionize['hemisphere_id'] == 0, 'hemisphere_id'] = 2
                
            structures = []
            hemis = []
            weight = []
            for structure_id in exp_unionize['structure_id'].unique():
                mask = mcc.get_structure_mask(structure_id)[0]
                halfdim = int(mask.shape[2]/2)
                l = mask[:,:,:halfdim]
                r = mask[:,:,halfdim:]
                blank = np.zeros_like(l)
                maskl = np.concatenate((l, blank), axis = 2)
                maskr = np.concatenate((blank, r), axis = 2)
                hemispheres = exp_unionize[exp_unionize['structure_id'] == structure_id]['hemisphere_id']
                for hemisphere in hemispheres:
                    if hemisphere == 1:
                        structures.append(structure_id)
                        hemis.append(hemisphere)
                        weight.append(np.sum(model*maskl))
                    if hemisphere == 2:
                        structures.append(structure_id)
                        hemis.append(hemisphere)
                        weight.append(np.sum(model*maskr))
            
            model_unionize = pd.DataFrame({'structure_id': structures,
                                           'hemisphere_id': hemis,
                                           'normalized_projection_volume': weight,
                                           'data': 'model'})
            unionizes = pd.concat([exp_unionize, model_unionize])
            unionizes.sort_values(by=['hemisphere_id', 'structure_id'], inplace = True)
            spearmanr, _ = st.spearmanr(unionizes[unionizes['data'] == 'experimental']['normalized_projection_volume'], 
                                        unionizes[unionizes['data'] == 'model']['normalized_projection_volume'])
            pearsonr, _ = st.pearsonr(unionizes[unionizes['data'] == 'experimental']['normalized_projection_volume'], 
                                      unionizes[unionizes['data'] == 'model']['normalized_projection_volume'])
            dat = {'image_series_id': isid, 
                   'spearman_correlation': spearmanr, 
                   'pearson_correlation': pearsonr,
                   'source': td_experiments[td_experiments['id'] == isid]['structure_abbrev'].values[0]}
            savepath = 'cluster_code/voxel_model/output'
            if platform.system() == 'Windows':
                savepath = r'C:\Users\jenniferwh\Dropbox (Personal)\DMN data\correlations\model_unionize_correlations'
            with open(os.path.join(savepath, '{}.json'.format(isid)), 'w') as outfile:
                json.dump(dat, outfile, sort_keys = False, indent = 4)
        except:
            continue
        

if __name__ == '__main__':

    main()
            

    
