#!/bin/bash
#
# Send daily webserver logs
#
# The weblogs need to be reviewed on a daily basis. However, that is
# often forgotten. This script will take a snapshot of the log for
# the day in question.
#
# The files were generated on Makana, but need to be sent from Ohana.

source $HOME/bin/cron/common.sh

temp_log=$HOME/$CHECK_WEB_LOG

mutt -s "Makana Webserver Logs for today"  bpg.shailen@gmail.com < $temp_log
