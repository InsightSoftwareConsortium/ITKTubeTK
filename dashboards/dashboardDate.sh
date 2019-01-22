#!/bin/sh
cd /Users/aylward/src/dashboards/
rm -rf ITKTubeTK
git clone git@github.com:/KitwareMedical/ITKTubeTK
cd ITKTubeTK/dashboards
echo `date` > dashboardDate.txt
git commit dashboardDate.txt -m "NIGHTLY: Update DashboardDate to initiate nightly CI builds"
git push
