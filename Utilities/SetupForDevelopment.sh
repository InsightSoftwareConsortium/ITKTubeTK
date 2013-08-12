#!/usr/bin/env bash

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
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

src_dir="${BASH_SOURCE%/*}/.."
cd "${src_dir}" &&
Utilities/GitSetup/setup-user && echo &&
Utilities/GitSetup/setup-hooks && echo &&
Utilities/GitSetup/tips

# Rebase master by default
git config rebase.stat true
git config branch.master.rebase true
git config hooks.KWStyle.conf "${src_dir}/CMake/KWStyle/KWStyle.kws.xml.in"
git config hooks.KWStyle.overwriteRulesConf "${src_dir}/CMake/KWStyle/KWStyle.Overwrite.txt.in"
