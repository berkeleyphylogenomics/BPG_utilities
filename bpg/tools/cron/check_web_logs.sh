#!/bin/bash
#
# Create email content for a message to be sent on a daily basis.
#
# The weblogs need to be reviewed on a daily basis. However, that is
# often forgotten. This script will take a snapshot of the log for
# the day in question.
#
# Since, for security reasons, Makana does not send email, email will
# be sent from Ohana at a later time. This is easy, as they share a 
# filesytem.
#
# This script was built to be run as a cron job. It makes the most
# sense to run it at the last minute of every day (since we are 
# grep'ing the snippets of log by timestamp. This has the potential,
# of course, to miss entries created between 11:59 p.m. (after the
# script has been run) and midnight. We accept this loss of alerting
# and race condition.
# Crontab entry:
# 59 23 * * * /clusterfs/ohana/software/prod/lib/python2.4/site-packages/bpg/tools/cron/check_web_logs.sh

source $HOME/bin/cron/common.sh

today=$(date +"\[%a %b %d")
temp_log=$HOME/$CHECK_WEB_LOG

cat > $temp_log << EOH
    This is the daily check of the webserver logs for $(date +"%d-%b").
The first section will list the production logs, and the second will
list the test logs. 

    We strive for zero entries in the production logs, but know
that the test environment will, by its very nature, have errors.
These should be resolved before rolling into production.


PRODUCTION ERROR LOGS
=====================

$(grep "$today" /var/log/httpd/production_error.log*)

STAGING ERROR LOGS
==================

$(grep "$today" /var/log/httpd/staging_error.log*)

EOH
