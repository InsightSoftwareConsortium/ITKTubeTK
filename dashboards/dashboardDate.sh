#!/bin/sh
rm -rf ITKTubeTK
git clone git@github.com:/KitwareMedical/ITKTubeTK
cd ITKTubeTK/dashboards
git pull --all
git checkout DashboardDate
echo `date` > dashboardDate.txt
git commit dashboardDate.txt -m "COMP: Update DashboardDate to initiate nightly CI builds"
git push
