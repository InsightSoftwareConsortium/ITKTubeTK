#!/bin/sh

echo `date` > .\dashboardDate.txt
git checkout -b DashboardDate
git commit dashboardDate.txt -m "Update DashboardDate to initiate night CI builds"
git push origin dashboardDate
