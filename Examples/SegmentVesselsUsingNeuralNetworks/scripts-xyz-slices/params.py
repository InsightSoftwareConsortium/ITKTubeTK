data_id = 'SVUNN'
filters = 32
design = 'xyz'

if design == 'xyz':
    train_bs = 24 * 1024 / filters
elif design == 'full3d':
    train_bs = 1024 / filters
else:
    raise ValueError("Bad value for NETWORK_DESIGN")

params = {
    "OUTPUT_DATA_ROOT": "/invalid/path/to/data",
    "NETWORK_DESIGN": design,
    "PATCH_RADIUS": 31,
    "SMOOTHING_RADIUS": 2,
    "POSITIVE_PATCHES_PER_INPUT_FILE": 5000,
    "NEGATIVE_TO_POSITIVE_RATIO": 1.0,
    "NUM_TRAIN_EPOCHS": 12,
    "NUM_TEST_EPOCHS": 1,
    "FILTERS_PER_LAYER": filters,
    # These parameters are set here in this file as a reminder that
    # they are guesses and are machine-specific.  These guesses will
    # not work automatically for all configurations, even on the
    # machine they were created for.
    "TRAIN_BATCH_SIZE": train_bs,
    "TEST_BATCH_SIZE": train_bs,
    "DEPLOY_BATCH_SIZE": 4 * train_bs / 3,
    "DEPLOY_TOP_WINDOW": 5,
    "SOLVER_TYPE": "SGD",
    "SOLVER_PARAMS": {
        "MOMENTUM": 0.9,
    },
    "BASE_LR": 0.2,
    "WEIGHT_DECAY": 0.0001,
    "GAMMA": 0.02,
    "VESSEL_SEED_PROBABILITY": 0.95,
}

if data_id == 'SVUNN':
    params.update({
        "INPUT_DATA_ROOT": "/home/nealsiekierski/caffe/data/SegmentVesselsUsingNeuralNetworks",
        "TYPE_SUBDIR_STRUCTURE": "*",
        "TYPES": {
            "Controls": "controls",
            "LargeTumor": "tumors",
        },
        "RESAMPLE_SPACING": None,
        "RESHUFFLE": False,
        "VESSEL_SCALE": 0.1,
    })
elif data_id == 'B-MRA':
    params.update({
        "INPUT_DATA_ROOT": "/home/nealsiekierski/data/Bullitt-MRA-converted/reflattened",
        "TYPE_SUBDIR_STRUCTURE": "",
        "TYPES": {
            "Male": "male",
            "Female": "female",
        },
        "RESAMPLE_SPACING": 0.513393,
        "RESHUFFLE": False,
        "VESSEL_SCALE": 0.9,
    })
elif data_id == 'stroke':
    params.update({
        "INPUT_DATA_ROOT": "/home/nealsiekierski/data/Lee-StrokeCollaterals",
        "TYPE_SUBDIR_STRUCTURE": "",
        "TYPES": {
            "restructured-with-skull": "data",
        },
        "RESAMPLE_SPACING": 0.5,
        "RESHUFFLE": True,
        "VESSEL_SCALE": 0.9,
    })
else:
    raise ValueError("Invalid value for data_id!")
