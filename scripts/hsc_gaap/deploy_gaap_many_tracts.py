import os, sys
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from hsc_gaap.deploy_gaap_array import deploy_training_job

# tracts = [8282, 8283, 8284, 8285, 8286, 9326, 9327, 9328, 9329, 8522, 8523, 8526, 8527, 8768, 8769, 9010, 9011, 9071, 9072, 9073,
#           9074, 9075, 9076, 9077, 9078, 9079, 9080, 9081, 9082, 9083, 9084, 9086, 9086, 9087, 9088, 9089, 9090, 9091, 9092, 9096,
#           9097,  9098,  9099,  9100,  9101,  9102,  9103,  9104,  9131,  9132,  9133,  9134,  9135,  9136]

# tracts = [9226,  9227,  9313,  9314,
#         9315,  9316,  9317,  9318,  9319,  9320,  9321,  9322,  9323,
#         9324,  9325,  9330,  9331,  9332,
#         9333,  9334,  9335,  9336,  9339,  9340,  9341,  9342,  9343,
#         9344]
# tracts = [9345,  9346,  9347,  9373,  9374,  9375,  9376,  9377]

tracts = [9378,  9379,  9455,  9456,  9457,  9458,  9459,  9465,  9466,
          9467,  9468,  9469,  9470,  9471,  9556,  9557,  9558,  9559,
          9560,  9561,  9562,  9563,  9564,  9565,  9568,  9586,  9587,
          9594,  9813,  9836, 10059, 10284]

for tract in tracts:
    deploy_training_job(tract, filter_jobs=5,
                        python_file='lambo/scripts/hsc_gaap/run_gaap.py',
                        name='gaap', email="am2907@princeton.edu", outname = None, 
                        repo='/scratch/gpfs/am2907/Merian/gaap', submit=True, fixpatches=False)