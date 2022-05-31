##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

"""Downnload experiments data from MIDAS.
"""


__license__ = "Apache License, Version 2.0"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import pydas
import os
import sys
import json
import logging
from optparse import OptionParser


LOGGING_LEVELS = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG}


def item_search(name, folder_id):
    """ Search item in folder

    :param name: Name of the item to search
    :param folder_id: ID of the folder where to search
    :returns: ID of the item (if found), or -1 (not found)
    """
    (dummy, found_item_id) = pydas.api._search_folder_for_item_or_folder(
        name,
        folder_id)
    return found_item_id


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = OptionParser()
    parser.add_option("", "--outDir", help="Directory where to store downloaded files", default=".")
    parser.add_option("", "--config", help="MIDAS configuration file")
    parser.add_option("", "--logto", help="Logging file")
    parser.add_option("", "--logat", help="Logging level (see code)", default="warning")
    (options, args) = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=LOGGING_LEVELS.get(options.logat, logging.NOTSET),
        filename=options.logto,
        format='%(asctime)s [%(funcName)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger()

    # Read MIDAS loging credentials
    if (options.config is None):
        print "Config file missing!"
        return -1
    config_fid = open(options.config).read()
    config = json.loads(config_fid)

    # Folder ID for the top-level data folder in MIDAS
    # Here: Designed Database of MR Brain Images of Healthy Individuals
    base_folder_id = '8051'

    pydas.login(
        email=config['email'],
        api_key=config['api_key'],
        url=config['url'])

    # Get attributes for top-level folder
    cur_folder = pydas.session.communicator.folder_get(
        pydas.session.token,
        base_folder_id)

    logger.debug("Output directory = %s" % options.outDir)
    download_dir = options.outDir;
    if not os.path.exists(download_dir):
        print "Directory %s not existent!" % download_dir

    # Currently, max(patientID) = 109
    patient_id_range = range(1,110)
    for patient_id in patient_id_range:
        folder_search_str = "Normal-%.3d" % patient_id
        (dummy, found_patient_folder_id) = pydas.api._search_folder_for_item_or_folder(
                folder_search_str,
                base_folder_id)

        if (found_patient_folder_id > 0):
            data_dir = os.path.join(download_dir, "Normal-%.3d" % patient_id)

            cur_children = pydas.session.communicator.folder_children(
                pydas.session.token,
                found_patient_folder_id)

            # Check if we have auxillary data available in the
            # current patient folder
            has_aux = False
            for sub_folder in cur_children['folders']:
                if sub_folder['name'] == 'AuxillaryData':
                    has_aux = True

            # In case we found auxillary data, download ...
            if has_aux:
                if (not os.path.exists(data_dir)):
                    os.makedirs(data_dir)

                logger.debug("Downloading patient data to %s" % data_dir)

                for sub_folder in cur_children['folders']:
                    item_id = -1

                    # We need the MRA data
                    if sub_folder['name'] == 'MRA':
                        target_file = "Normal%.3d-MRA.mha" % patient_id
                        item_id = item_search(target_file, sub_folder['folder_id'])
                        if (item_id > 0):
                            pydas.api._download_item(item_id, data_dir)
                        else:
                            logger.warning("%s not found!" % target_file)

                    # ... ,the T1-Flash data
                    elif sub_folder['name'] == 'T1-Flash':
                        target_file = "Normal%.3d-T1-Flash.mha" % patient_id
                        item_id = item_search(target_file, sub_folder['folder_id'])
                        if (item_id > 0):
                            pydas.api._download_item(item_id, data_dir)
                        else:
                            logger.warning("%s not found!" % target_file)

                    # and the 1) extracted vascular network data (.tre) as well
                    # as the skull-stripped T1-Flash MR images
                    elif sub_folder['name'] == 'AuxillaryData':
                        item_id = item_search("VascularNetwork.tre",
                            sub_folder['folder_id'])
                        if (item_id > 0):
                            pydas.api._download_item(item_id, data_dir)
                        else:
                            logger.warning("Could not find VascularNetwork.tre")

                        item_id = item_search("SkullStripped-T1-Flash.mha",
                            sub_folder['folder_id'])
                        if (item_id > 0):
                            pydas.api._download_item(item_id, data_dir)
                        else:
                            logger.warning("Could not find SkullStripped-T1-Flash.mha")
            else:
                logger.warning("No AuxillaryData directory, skipping patient %d!" % patient_id)


if __name__ == "__main__":
    sys.exit( main() )
