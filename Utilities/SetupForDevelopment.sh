#!/usr/bin/env bash

cd "${BASH_SOURCE%/*}/.." &&
CMake/GitSetup/setup-user && echo &&
CMake/GitSetup/setup-hooks && echo &&
CMake/GitSetup/setup-aliases && echo &&
CMake/GitSetup/setup-stage && echo &&
(CMake/GitSetup/setup-ssh ||
 echo 'Failed to setup SSH.  Run this again to retry.') && echo &&
(CMake/GitSetup/setup-gerrit ||
 echo 'Failed to setup Gerrit.  Run this again to retry.') && echo &&
CMake/GitSetup/tips

# Rebase master by default
git config rebase.stat true
git config branch.master.rebase true
