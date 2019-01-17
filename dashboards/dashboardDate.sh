#!/bin/sh
cd /Users/aylward/src/dashboards/ITKTubeTK/dashboards
git checkout master --force
git pull --all
git checkout DashboardDate
echo `date` > dashboardDate.txt
git commit dashboardDate.txt -m "COMP: Update DashboardDate to initiate nightly CI builds"
git push
