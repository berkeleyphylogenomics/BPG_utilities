#!/usr/bin/env python

import logging
import pickle
import sys
import os

if __name__ == '__main__':
    for path in sys.argv[2:]:
        sys.path.append(path)

from django.core.management import setup_environ
from django.core.mail import send_mail

from pfacts003 import settings

# Constants
START_LEVEL = 23
STOP_LEVEL = 24
STEP_START_LEVEL = 26
STEP_COMPLETE_LEVEL = 27

# Module configuration (including django's environ)
setup_environ(settings)
logging.addLevelName(START_LEVEL, "START")
logging.addLevelName(STEP_START_LEVEL, "STEP_START")
logging.addLevelName(STEP_COMPLETE_LEVEL, "STEP_COMPLETE")
logging.addLevelName(STOP_LEVEL, "STOP")

from pfacts003.queued.utils import arg_name_re, _q_import, \
    q_decorated_function, Q_data

def main(app_name):
    # import the app
    app = _q_import(app_name)

    # parse pickled input from the form
    handle = open('./cleaned_data.pkl')
    cleaned_data = pickle.load(handle)
    handle.close()
    q_data = Q_data(cleaned_data, app_name)

    # make generic handler to be used in both loggers
    handler = logging.FileHandler('./%s.log' % app_name)
    handler.setLevel(logging.DEBUG)

    #handler.setFormatter(logging.Formatter(' - '.join('%%(%s)s'% x \
    #                                for x in log_fields), '%Y-%m-%d %H:%M:%S'))

    handler.setFormatter(logging.Formatter(' - '.join( \
        # put special brackets around the 'message' field
        (x != 'message' and '%%(%s)s' or '<<< %%(%s)s >>>') % x \
        for x in ('asctime', 'name', 'levelname', 'message')), '%Y-%m-%d %H:%M:%S'))

    # make queued_logger for '%s.out' % program_name
    queued_log = logging.getLogger('application')
    queued_log.addHandler(handler)
    queued_log.setLevel(logging.DEBUG)

    # Create step_log master logger
    step_logger = logging.getLogger('step')
    step_logger.setLevel(logging.DEBUG)
    step_logger.addHandler(handler)

    steps = sorted([
        getattr(app, item) \
        for item in dir(app) \
        if isinstance(getattr(app, item), q_decorated_function)
    ], lambda x,y: x.count - y.count)

    # prepare email parameters/messages
    mail_params = {'app_title': app.title, 'app_name': app_name,
        'job_id': os.path.basename(os.getcwd()).replace(app_name,'',1)}
    mail_error = '''Your %(app_title)s submission has encountered an error.

Information regarding this issue can be found at the following link:
http://makana.berkeley.edu/q/%(app_name)s/%(job_id)s/
''' % mail_params
    mail_success = '''Your %(app_title)s submission has finished.

Your results can be found at the following link:
http://makana.berkeley.edu/q/%(app_name)s/%(job_id)s/
''' % mail_params


    # Log start on application level
    queued_log.log(START_LEVEL, "'%s' has started." % app.title)

    for step in steps:
        # Create specific logger for this step (inherit from 'step' logger)
        step_log = logging.getLogger('step.%s' % step.label)
        step_log.log(STEP_START_LEVEL, "Step '%s' has started." % step.title)

        # Pass logger to function for local function use
        q_data._set_logger(step_log)
        q_data._set_step(step.label)

        try:
            step_return = step(q_data)
        except Exception, e:
            step_log.error("Error occurred at step '%s'" % step.title)
            # step_log.exception(e)
            # Log completion on application level
            queued_log.error(app.title)
            if cleaned_data.has_key('email_to') and \
                bool(cleaned_data['email_to']):
                send_mail(cleaned_data['email_subject'],
                    mail_error,
                    'phylo@phylogenomics.berkeley.edu',
                    [cleaned_data['email_to']], fail_silently=False,
                )
            raise e

        if step_return is False:
            step_log.error("Error occurred at step '%s'" % step.title)
            # Log completion on application level
            # queued_log.error(app.title)
            if cleaned_data.has_key('email_to') and \
                bool(cleaned_data['email_to']):
                send_mail(cleaned_data['email_subject'],
                    mail_error,
                    'phylo@phylogenomics.berkeley.edu',
                    [cleaned_data['email_to']], fail_silently=False,
                )
            break
        else:
            step_log.log(STEP_COMPLETE_LEVEL,
                 "Step '%s' has completed successfully." % step.title)
    else:
        # Log completion on application level
        queued_log.log(STOP_LEVEL, app.title)
        if cleaned_data.has_key('email_to') and \
            bool(cleaned_data['email_to']):
            send_mail(cleaned_data['email_subject'],
                mail_success,
                'phylo@phylogenomics.berkeley.edu',
                [cleaned_data['email_to']], fail_silently=False,
            )

if __name__ == '__main__':
    if len(sys.argv) > 1:
        app_name = sys.argv[1]
        if app_name not in settings.QUEUED_APPS:
            sys.exit(1)
        else:
            main(app_name)
