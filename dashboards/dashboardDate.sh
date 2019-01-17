#!/bin/sh
cd /Users/aylward/src/ITKTubeTK
git checkout master --force
git pull --all
git checkout DashboardDate
echo `date` > /Users/aylward/src/dashboards/ITKTubeTK/dashboards/dashboardDate.txt
git commit dashboardDate.txt -m "Update DashboardDate to initiate night CI builds"
git push origin dashboardDate
