# -*- coding: utf-8 -*-
"""
Created on Sun Nov 05 14:18:41 2017

@author: jenniferwh
"""
import os
import nrrd
import numpy as np

from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

print('first test')
output_directory = r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/DMN paper/remove_CAV_signal/grid'
datpaths = {561511939: u'/projects/ctyconn/vol1/prod0/image_series_561511939/', 614738437: u'/allen/programs/celltypes/production/ctyconn/prod20/image_series_614738437/', 543297033: u'/allen/programs/celltypes/production/0378/prod447/image_series_543297033/', 605092364: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605092364/', 523481517: u'/allen/programs/celltypes/production/0378/prod463/image_series_523481517/', 551737360: u'/allen/programs/celltypes/production/0378/prod478/image_series_551737360/', 525412369: u'/allen/programs/celltypes/production/0378/prod465/image_series_525412369/', 544454674: u'/allen/programs/celltypes/production/0378/prod449/image_series_544454674/', 477435412: u'/allen/programs/celltypes/production/0378/prod454/image_series_477435412/', 560045081: u'/allen/programs/celltypes/production/0378/prod455/image_series_560045081/', 501785691: u'/allen/programs/celltypes/production/0378/prod423/image_series_501785691/', 595459108: u'/allen/programs/celltypes/production/ctyconn/prod9/image_series_595459108/', 521255975: u'/allen/programs/celltypes/production/0378/prod461/image_series_521255975/', 607321130: u'/allen/programs/celltypes/production/ctyconn/prod16/image_series_607321130/', 522274859: u'/allen/programs/celltypes/production/0378/prod462/image_series_522274859/', 569904687: u'/projects/ctyconn/vol1/prod1/image_series_569904687/', 529129011: u'/allen/programs/celltypes/production/0378/prod470/image_series_529129011/', 613898292: u'/allen/programs/celltypes/production/ctyconn/prod19/image_series_613898292/', 496965687: u'/allen/programs/celltypes/production/0378/prod458/image_series_496965687/', 501923898: u'/allen/programs/celltypes/production/0378/prod424/image_series_501923898/', 523704892: u'/allen/programs/celltypes/production/0378/prod464/image_series_523704892/', 607316031: u'/allen/programs/celltypes/production/ctyconn/prod16/image_series_607316031/', 609475139: u'/allen/programs/celltypes/production/ctyconn/prod17/image_series_609475139/', 475616836: u'/allen/programs/celltypes/production/0378/prod453/image_series_475616836/', 495877413: u'/allen/programs/celltypes/production/0378/prod418/image_series_495877413/', 595859014: u'/allen/programs/celltypes/production/ctyconn/prod9/image_series_595859014/', 585761377: u'/allen/programs/celltypes/production/ctyconn/prod23/image_series_585761377/', 617900105: u'/allen/programs/celltypes/production/ctyconn/prod20/image_series_617900105/', 573977678: u'/projects/ctyconn/vol1/prod5/image_series_573977678/', 561307215: u'/projects/ctyconn/vol1/prod1/image_series_561307215/', 561918904: u'/projects/ctyconn/vol1/prod1/image_series_561918904/', 527390805: u'/allen/programs/celltypes/production/0378/prod466/image_series_527390805/', 527713883: u'/allen/programs/celltypes/production/0378/prod467/image_series_527713883/', 502956560: u'/allen/programs/celltypes/production/0378/prod425/image_series_502956560/', 528964707: u'/allen/programs/celltypes/production/0378/prod442/image_series_528964707/', 523177830: u'/allen/programs/celltypes/production/0378/prod463/image_series_523177830/', 539511058: u'/allen/programs/celltypes/production/0378/prod446/image_series_539511058/', 592724077: u'/allen/programs/celltypes/production/ctyconn/prod8/image_series_592724077/', 478995566: u'/allen/programs/celltypes/production/0378/prod456/image_series_478995566/', 616674416: u'/allen/programs/celltypes/production/ctyconn/prod20/image_series_616674416/', 591394174: u'/allen/programs/celltypes/production/ctyconn/prod27/image_series_591394174/', 589726838: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589726838/', 526515319: u'/allen/programs/celltypes/production/0378/prod440/image_series_526515319/', 539323512: u'/allen/programs/celltypes/production/0378/prod445/image_series_539323512/', 501711996: u'/allen/programs/celltypes/production/0378/prod422/image_series_501711996/', 524063357: u'/allen/programs/celltypes/production/0378/prod464/image_series_524063357/', 567723369: u'/projects/ctyconn/vol1/prod5/image_series_567723369/', 485929089: u'/allen/programs/celltypes/production/0378/prod459/image_series_485929089/', 528511254: u'/allen/programs/celltypes/production/0378/prod469/image_series_528511254/', 596281479: u'/allen/programs/celltypes/production/ctyconn/prod10/image_series_596281479/', 475829896: u'/allen/programs/celltypes/production/0378/prod454/image_series_475829896/', 483094671: u'/allen/programs/celltypes/production/0378/prod458/image_series_483094671/', 563352720: u'/projects/ctyconn/vol1/prod2/image_series_563352720/', 572388976: u'/projects/ctyconn/vol1/prod14/image_series_572388976/', 560965104: u'/projects/ctyconn/vol1/prod0/image_series_560965104/', 545799841: u'/allen/programs/celltypes/production/0378/prod449/image_series_545799841/', 609475867: u'/allen/programs/celltypes/production/ctyconn/prod17/image_series_609475867/', 605094567: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605094567/', 528512680: u'/allen/programs/celltypes/production/0378/prod469/image_series_528512680/', 502955689: u'/allen/programs/celltypes/production/0378/prod425/image_series_502955689/', 536920234: u'/allen/programs/celltypes/production/0378/prod445/image_series_536920234/', 561986735: u'/projects/ctyconn/vol1/prod1/image_series_561986735/', 572390577: u'/projects/ctyconn/vol1/prod3/image_series_572390577/', 523180728: u'/allen/programs/celltypes/production/0378/prod437/image_series_523180728/', 475828414: u'/allen/programs/celltypes/production/0378/prod454/image_series_475828414/', 585931968: u'/allen/programs/celltypes/production/ctyconn/prod23/image_series_585931968/', 565146821: u'/projects/ctyconn/vol1/prod3/image_series_565146821/', 609157409: u'/allen/programs/celltypes/production/ctyconn/prod17/image_series_609157409/', 521955016: u'/allen/programs/celltypes/production/0378/prod462/image_series_521955016/', 595898061: u'/allen/programs/celltypes/production/ctyconn/prod9/image_series_595898061/', 479115470: u'/allen/programs/celltypes/production/0378/prod456/image_series_479115470/', 585022673: u'/projects/ctyconn/vol1/prod23/image_series_585022673/', 605093074: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605093074/', 589322070: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589322070/', 605685976: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605685976/', 601885751: u'/allen/programs/celltypes/production/ctyconn/prod11/image_series_601885751/', 572595932: u'/projects/ctyconn/vol1/prod3/image_series_572595932/', 607289053: u'/allen/programs/celltypes/production/ctyconn/prod16/image_series_607289053/', 605496542: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605496542/', 596575967: u'/allen/programs/celltypes/production/ctyconn/prod10/image_series_596575967/', 561918178: u'/projects/ctyconn/vol1/prod1/image_series_561918178/', 561512675: u'/projects/ctyconn/vol1/prod0/image_series_561512675/', 589702885: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589702885/', 592698087: u'/allen/programs/celltypes/production/ctyconn/prod8/image_series_592698087/', 531443949: u'/allen/programs/celltypes/production/0378/prod443/image_series_531443949/', 528741104: u'/allen/programs/celltypes/production/0378/prod469/image_series_528741104/', 589064435: u'/allen/programs/celltypes/production/ctyconn/prod24/image_series_589064435/', 526784559: u'/allen/programs/celltypes/production/0378/prod440/image_series_526784559/', 527393013: u'/allen/programs/celltypes/production/0378/prod466/image_series_527393013/', 521617657: u'/allen/programs/celltypes/production/0378/prod461/image_series_521617657/', 525413115: u'/allen/programs/celltypes/production/0378/prod465/image_series_525413115/', 567301515: u'/projects/ctyconn/vol1/prod0/image_series_567301515/', 539514325: u'/allen/programs/celltypes/production/0378/prod446/image_series_539514325/', 520615681: u'/allen/programs/celltypes/production/0378/prod460/image_series_520615681/', 553080579: u'/allen/programs/celltypes/production/0378/prod478/image_series_553080579/', 502592260: u'/allen/programs/celltypes/production/0378/prod425/image_series_502592260/', 601359621: u'/allen/programs/celltypes/production/ctyconn/prod10/image_series_601359621/', 584651014: u'/projects/ctyconn/vol1/prod22/image_series_584651014/', 617898760: u'/allen/programs/celltypes/production/ctyconn/prod20/image_series_617898760/', 546389260: u'/allen/programs/celltypes/production/0378/prod450/image_series_546389260/', 524266253: u'/allen/programs/celltypes/production/0378/prod465/image_series_524266253/', 589397775: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589397775/', 575683857: u'/projects/ctyconn/vol1/prod6/image_series_575683857/', 595261714: u'/allen/programs/celltypes/production/ctyconn/prod8/image_series_595261714/', 569932566: u'/projects/ctyconn/vol1/prod1/image_series_569932566/', 523457583: u'/allen/programs/celltypes/production/0378/prod463/image_series_523457583/', 571647261: u'/projects/ctyconn/vol1/prod3/image_series_571647261/', 501786400: u'/allen/programs/celltypes/production/0378/prod423/image_series_501786400/', 526783792: u'/allen/programs/celltypes/production/0378/prod466/image_series_526783792/', 571652998: u'/projects/ctyconn/vol1/prod3/image_series_571652998/', 529428776: u'/allen/programs/celltypes/production/0378/prod470/image_series_529428776/', 595259180: u'/allen/programs/celltypes/production/ctyconn/prod8/image_series_595259180/', 595458351: u'/allen/programs/celltypes/production/ctyconn/prod9/image_series_595458351/', 527711024: u'/allen/programs/celltypes/production/0378/prod467/image_series_527711024/', 571653937: u'/projects/ctyconn/vol1/prod13/image_series_571653937/', 515920693: u'/allen/programs/celltypes/production/0378/prod431/image_series_515920693/', 524267323: u'/allen/programs/celltypes/production/0378/prod465/image_series_524267323/', 563205064: u'/projects/ctyconn/vol1/prod2/image_series_563205064/', 569993539: u'/projects/ctyconn/vol1/prod1/image_series_569993539/', 604101445: u'/allen/programs/celltypes/production/ctyconn/prod13/image_series_604101445/', 521954271: u'/allen/programs/celltypes/production/0378/prod462/image_series_521954271/', 637855050: u'/allen/programs/celltypes/production/ctyconn/prod22/image_series_637855050/', 475830603: u'/allen/programs/celltypes/production/0378/prod454/image_series_475830603/', 572771153: u'/projects/ctyconn/vol1/prod3/image_series_572771153/', 580594514: u'/projects/ctyconn/vol1/prod20/image_series_580594514/', 502005076: u'/allen/programs/celltypes/production/0378/prod459/image_series_502005076/', 501837158: u'/allen/programs/celltypes/production/0378/prod423/image_series_501837158/', 475617622: u'/allen/programs/celltypes/production/0378/prod453/image_series_475617622/', 570461017: u'/projects/ctyconn/vol1/prod2/image_series_570461017/', 478877019: u'/allen/programs/celltypes/production/0378/prod455/image_series_478877019/', 502590301: u'/allen/programs/celltypes/production/0378/prod425/image_series_502590301/', 589399902: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589399902/', 484612961: u'/allen/programs/celltypes/production/0378/prod459/image_series_484612961/', 601804603: u'/allen/programs/celltypes/production/ctyconn/prod11/image_series_601804603/', 528732005: u'/allen/programs/celltypes/production/0378/prod469/image_series_528732005/', 571410278: u'/projects/ctyconn/vol1/prod2/image_series_571410278/', 496964969: u'/allen/programs/celltypes/production/0378/prod458/image_series_496964969/', 591168591: u'/allen/programs/celltypes/production/ctyconn/prod26/image_series_591168591/', 570071403: u'/projects/ctyconn/vol1/prod1/image_series_570071403/', 531233132: u'/allen/programs/celltypes/production/0378/prod472/image_series_531233132/', 571816813: u'/projects/ctyconn/vol1/prod3/image_series_571816813/', 552543088: u'/allen/programs/celltypes/production/0378/prod478/image_series_552543088/', 614435699: u'/allen/programs/celltypes/production/ctyconn/prod19/image_series_614435699/', 518012479: u'/allen/programs/celltypes/production/0378/prod433/image_series_518012479/', 518605181: u'/allen/programs/celltypes/production/0378/prod460/image_series_518605181/', 606278526: u'/allen/programs/celltypes/production/ctyconn/prod15/image_series_606278526/', 475616128: u'/allen/programs/celltypes/production/0378/prod453/image_series_475616128/', 475733382: u'/allen/programs/celltypes/production/0378/prod453/image_series_475733382/', 527713163: u'/allen/programs/celltypes/production/0378/prod467/image_series_527713163/', 574637452: u'/projects/ctyconn/vol1/prod16/image_series_574637452/', 568502684: u'/projects/ctyconn/vol1/prod1/image_series_568502684/', 601904029: u'/allen/programs/celltypes/production/ctyconn/prod11/image_series_601904029/', 605854110: u'/allen/programs/celltypes/production/ctyconn/prod15/image_series_605854110/', 569623967: u'/projects/ctyconn/vol1/prod1/image_series_569623967/', 595890081: u'/allen/programs/celltypes/production/ctyconn/prod10/image_series_595890081/', 485553574: u'/allen/programs/celltypes/production/0378/prod415/image_series_485553574/', 592522663: u'/allen/programs/celltypes/production/ctyconn/prod7/image_series_592522663/', 571401645: u'/projects/ctyconn/vol1/prod2/image_series_571401645/', 531397136: u'/allen/programs/celltypes/production/0378/prod443/image_series_531397136/', 604100536: u'/allen/programs/celltypes/production/ctyconn/prod13/image_series_604100536/', 617901499: u'/allen/programs/celltypes/production/ctyconn/prod20/image_series_617901499/', 572753855: u'/projects/ctyconn/vol1/prod3/image_series_572753855/', 623838656: u'/allen/programs/celltypes/production/ctyconn/prod21/image_series_623838656/', 578332611: u'/projects/ctyconn/vol1/prod7/image_series_578332611/', 601900484: u'/allen/programs/celltypes/production/ctyconn/prod11/image_series_601900484/', 623286726: u'/allen/programs/celltypes/production/ctyconn/prod21/image_series_623286726/', 528328648: u'/allen/programs/celltypes/production/0378/prod468/image_series_528328648/', 504176074: u'/allen/programs/celltypes/production/0378/prod427/image_series_504176074/', 592698832: u'/allen/programs/celltypes/production/ctyconn/prod8/image_series_592698832/', 478678606: u'/allen/programs/celltypes/production/0378/prod455/image_series_478678606/', 589398486: u'/allen/programs/celltypes/production/ctyconn/prod25/image_series_589398486/', 501883865: u'/allen/programs/celltypes/production/0378/prod424/image_series_501883865/', 589065144: u'/allen/programs/celltypes/production/ctyconn/prod24/image_series_589065144/', 607059419: u'/allen/programs/celltypes/production/ctyconn/prod16/image_series_607059419/', 528511967: u'/allen/programs/celltypes/production/0378/prod469/image_series_528511967/', 530000865: u'/allen/programs/celltypes/production/0378/prod471/image_series_530000865/', 560029094: u'/allen/programs/celltypes/production/0378/prod455/image_series_560029094/', 561506791: u'/projects/ctyconn/vol1/prod1/image_series_561506791/', 495346667: u'/allen/programs/celltypes/production/0378/prod416/image_series_495346667/', 563731437: u'/projects/ctyconn/vol1/prod2/image_series_563731437/', 561910766: u'/projects/ctyconn/vol1/prod0/image_series_561910766/', 606260719: u'/allen/programs/celltypes/production/ctyconn/prod15/image_series_606260719/', 555012592: u'/allen/programs/celltypes/production/0378/prod478/image_series_555012592/', 518013943: u'/allen/programs/celltypes/production/0378/prod433/image_series_518013943/', 501787135: u'/allen/programs/celltypes/production/0378/prod423/image_series_501787135/', 529435133: u'/allen/programs/celltypes/production/0378/prod470/image_series_529435133/', 605112318: u'/allen/programs/celltypes/production/ctyconn/prod14/image_series_605112318/', 573331541: u'/projects/ctyconn/vol1/prod4/image_series_573331541/'}

mcc = MouseConnectivityCache(resolution = 10)
print('mcc')
structure_tree = mcc.get_structure_tree()
iso = structure_tree.get_structures_by_acronym(['Isocortex'])[0]
iso_mask = mcc.get_structure_mask(iso['id'])[0]
print('iso_mask')

def remove_CAV_projections():
    
    for imser in [475616128]:
        print(imser)
        grid_path = os.path.join(datpaths[imser], 'grid')
        proj_density, _ = nrrd.read(os.path.join(grid_path, 'projection_density_10.nrrd'))
        CAV_density, _ = nrrd.read(os.path.join(grid_path, 'CAV_density_10.nrrd'))  
        data_mask, _ = nrrd.read(os.path.join(grid_path, 'data_mask_10.nrrd'))
        masked_density = np.multiply(proj_density, data_mask)
        masked_density[np.where(CAV_density > 0)] = np.nan
        iso_masked = np.copy(iso_mask)
        iso_masked[np.where(CAV_density > 0)] = np.nan
        out_dir = os.path.join(
            output_directory, '{0}'.format(imser)
            )
        print(out_dir)
        output_mask_path = os.path.join(out_dir, 'iso_masked_10.nrrd')
        output_nrrd_path = os.path.join(out_dir, 'masked_density_10.nrrd')
        nrrd.write(output_mask_path, iso_masked)
        nrrd.write(output_nrrd_path, masked_density)

remove_CAV_projections()