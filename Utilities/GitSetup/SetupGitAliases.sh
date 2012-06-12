#!/usr/bin/env bash
#==========================================================================
#
#   Copyright Insight Software Consortium
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          http://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/

echo "Setting up useful Git aliases..." &&

# Pull all updates - first a general pull and then submodules.
git config alias.pullall "!bash -c \"git pull && git submodule update --init\"" &&
# Useful alias to see what commits are on the current branch with respect
# to origin/master.
git config alias.prepush 'log --graph --stat origin/master..' &&

# Staging aliases - can help streamline staging topic branches.
stage_cmd='ssh git@igstk.org stage IGSTK' &&
git_branch="\$(git symbolic-ref HEAD | sed -e 's|^refs/heads/||')" &&
# General alias to run the SSH command, e.g. git stage-cmd print.
git config alias.stage-cmd "!${stage_cmd}" &&
# Push the current topic branch to the stage.
git config alias.stage-push "!bash -c \"git fetch stage --prune && git push stage HEAD\"" &&
# List all staged topic branches.
git config alias.stage-branch "!bash -c \"${stage_cmd} print\"" &&
# Merge the current topic branch (if staged) into the next branch.
git config alias.stage-merge "!bash -c \"${stage_cmd} merge ${git_branch}\"" &&

true
