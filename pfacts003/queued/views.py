import os
import tempfile
import re
import pickle
import stat

from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseForbidden, \
    HttpResponseRedirect, Http404
from django.template import Context, loader
from django.forms import Form

from pfacts003.queued.forms import EmailForm
from django.conf import settings
#from pfacts003.queued.config import work_prefix, job_prefix, \
#    recaptcha_private, recaptcha_public, queued_run_prefix, site_packages_prefix
from pfacts003.queued.utils import _q_import, q_job_is_finished, q_read_log
from recaptcha.client import captcha

def return_forbidden(message):
    """Return 403.html template with 403 header"""

    t = loader.get_template('403.html')
    return HttpResponseForbidden(t.render(Context({'message': message})))

def q_index(request, app_name):

    # import the app
    try:
        app = _q_import(app_name)
    except ImportError, e:
        return return_forbidden("Unable to import %s." % app_name)
        
    response_dict = {
        'js': {'jquery': True, 'jquery_ui': True, 'expandable': True},
        'main_viewing': app_name, 'sub_viewing': 'submission_form',
        'recaptcha_public': settings.QUEUED_RECAPTCHA_PUBLIC,
    }

    advanced_form = getattr(app, 'AdvancedForm', Form)(request.POST,
        request.FILES)

    # check for submission
    if not request.POST:

        response_dict.update(
            form=app.AppForm(),
            # replace the email subject, if indicated
            email_form= hasattr(app, 'email_subject') and \
                EmailForm({'email_subject': app.email_subject}) or EmailForm(),
            # provide the advanced form if it exists, otherwise None
            advanced_form=hasattr(app, 'AdvancedForm') and app.AdvancedForm() \
            or None)
        return render_to_response(
            '%s/index.html' % app_name,
            response_dict,
        )

    # load the form page with corrections
    form = app.AppForm(request.POST, request.FILES)
    email_form = EmailForm(request.POST)
    # assign advanced_form, or give it a blank form (for convenience)
    advanced_form = getattr(app, 'AdvancedForm', Form)(request.POST,
        request.FILES)
    if not form.is_valid() or not email_form.is_valid() or \
        not advanced_form.is_valid():
        response_dict.update(
            form=form,
            email_form=email_form,
            # pass False if invalid, otherwise load the error code
            # pass None if there is no AdvancedForm for the app
            advanced_form=hasattr(app,'AdvancedForm') and advanced_form or None,
        )
        return render_to_response(
            '%s/index.html' % app_name,
            response_dict,
        )

    # make temp folder for running the job and set permissions
    work_path = tempfile.mkdtemp(dir=settings.QUEUED_INCOMING_DIR, prefix=app_name)
    work_dir = os.path.basename(work_path)
    # TODO !! FOR TESTING; MAKE A POINT TO FIX !! #
    os.chmod(work_path, stat.S_IRWXU | stat.S_IRWXG)

    # pickle the cleaned_data
    handle = open(os.path.join(work_path,'cleaned_data.pkl'), 'w')
    # pickle the combined dictionaries
    pickle.dump(dict(form.cleaned_data, **dict(email_form.cleaned_data,
        **advanced_form.cleaned_data)), handle)
    handle.close()

    # write the scripts
    #qsub_exec = os.path.join(os.path.split(os.path.abspath(__file__))[0],
    #                                                           'queued_run.py')
    qsub_exec = os.path.join(settings.QUEUED_DEPLOYMENT_ROOT, 'pfacts003', 'queued', 'queued_run.py')

    handle = open(os.path.join(settings.QUEUED_QSUB_DIR, work_dir) + '.sh', 'w')
    handle.write('%s\n%s %s %s %s' % (work_path, qsub_exec, app_name, settings.QUEUED_DEPLOYMENT_ROOT, settings.QUEUED_SITE_PACKAGES_ROOT))
    handle.close()

    # redirect to the execution/results page
    return HttpResponseRedirect('/q/%s/%s/' % (
        app_name, work_dir.replace(app_name, '', 1)))

def q_general_results(request, app_name, rand):

    # construct paths for log files
    work_dir = app_name + rand
    work_path = os.path.join(settings.QUEUED_INCOMING_DIR, work_dir)
    log_path = os.path.join(work_path, '%s.log' % app_name)

    # make sure the job ID is valid
    if not os.path.exists(work_path):
        raise Http404

    # import the app
    try:
        app = _q_import(app_name)
    except ImportError, e:
        return return_forbidden('Unable to import %s.' % app_name)

    # read the log
    job_info = q_read_log(app, work_path)

    response_dict = {'main_viewing': app_name}

    response_dict.update(job_info=job_info,
        js={'jquery': True, 'jquery_ui': True, 'expandable': True})
    if not job_info.errored and not job_info.is_finished():
        response_dict.update(refresh_time=15)
    
    if not job_info.is_finished():
        return render_to_response('%s/exec.html' % app_name, response_dict)

    return app.general_results(request, work_path, dict(response_dict,
        refresh_time=False))

def q_specific_results(request, app_name, rand, specific):
    '''broken for the time being'''
    
    # Construct paths
    work_dir = app_name + rand
    log_path = os.path.join(settings.QUEUED_INCOMING_DIR, work_dir, '%s.log' % app_name)

    # Check for job completion
    if os.path.exists(log_path):
        if q_job_is_finished(_q_import(app_name), work_dir):
            '''totally broken  vvv'''
            return __import__(app_name + '.queued.specific_results')(
                request, app_name, work_path.replace(app_name, 1), specific)
    
    raise Http404

def text_download(text, filename):
    response = HttpResponse(mimetype='text/plain')
    response['Content-Disposition'] = \
        'attachment; filename=%s' % filename
    print >>response, text,
    return response

def q_download(request, app_name, rand, text_dl, dummy_download):

    # Construct paths
    work_dir = app_name + rand
    work_path = os.path.join(settings.QUEUED_INCOMING_DIR, work_dir)
    log_path = os.path.join(work_path, '%s.log' % app_name)
    file_path = os.path.join(work_path, text_dl)

    # Check that the file is available for download
    try:
        app = _q_import(app_name)
    except ImportError, e:
        return return_forbidden("Unable to import %s." % app_name)

    if text_dl not in app.downloadable_files:
        return return_forbidden("%s not in downloadable files." % text_dl)
            
    if not os.path.exists(log_path) or not os.path.exists(file_path) or \
        not q_job_is_finished(app, work_path):
        raise Http404

    file = open(file_path, 'r')
    content = file.read()
    file.close()

    return text_download(content, str(text_dl))
