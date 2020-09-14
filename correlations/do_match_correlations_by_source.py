import sys
import platform
if platform.system() == 'Linux':
    sys.path.append(r'/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/pbstools')

from pbstools import PythonJob

python_file = r"/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/correlations/correlate_matches_by_source.py"


jobdir = '/allen/programs/celltypes/workgroups/mousecelltypes/T503_Connectivity_in_Alzheimer_Mice/Jennifer/cluster_code/cluster_out'

job_settings = {'queue': 'celltypes',
                'jobdir': jobdir,
                'mem':'20g',
                'walltime':'8:00:00',
                }

sources = ['SUB', 'AV', 'VISpl', 'AIp', 'GU', 'VISal', 'AUDpo', 'VISli', 
           'VISam', 'VISC', 'SSs', 'VISpm', 'PL', 'VISa', 'VISpor',
           'RSPd', 'SSp-m', 'AUDd', 'AUDp', 'SSp-ll', 'ORBvl', 'SSp-n', 'VISrl', 
           'AM', 'CA1', 'RE', 'AIv', 'CA3', 'AId', 'ORBm', 'MOp', 'SSp-un', 
           'SSp-ul', 'ProS', 'SSp-tr', 'ORBl', 'ILA', 'POST', 'APr', 'PAR', 
           'PRE', 'ECT', 'TEa', 'FRP', 'AUDv', 'CLA', 'PERI', 'CA2', 'FC', 'IG', 
           'HATA', 'SSp-bfd', 'RSPv', 'MOs', 'ACAd', 'VISl', 'ACAv', 'VISp',  
           'ENTl', 'ENTm', 'DG', 'RSPagl', 'MD']


'''
sources = ['VAL', 'VM', 'VPL', 'VPLpc', 'VPM', 'VPMpc', 'PoT', 'SPFm', 'SPFp', 
            'SPA', 'PP', 'MG', 'LGd', 'LP', 'PO', 'POL', 'SGN', 
            'AD', 'IAM', 'IAD', 'LD', 'IMD', 'SMT', 'PR', 'PVT', 'PT', 
            'Xi', 'RH', 'CM', 'PCN', 'CL', 'PF', 'PIL', 'RT', 'IGL', 
            'IntG', 'LGv', 'SubG', 'MH', 'LH']
'''
for source in sources:
    PythonJob(
        python_file,
        python_args= source, 
        conda_env='drjigsaw',
        jobname='match_corr_{}'.format(source),
        **job_settings
    ).run(dryrun=False)