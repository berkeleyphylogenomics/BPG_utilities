import datetime
import glob
import os
import shutil
import stat

from django.core.mail import send_mail
from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseRedirect

from pfacts003.common.forms import ContactForm

def index(request):
    return render_to_response('common/index.html')

def contact_us(request):
    '''
    Currently this function is not being used since the url has been redirected to
    a simple template. However, this is a nice function, hence I have left it here for future use.
    '''
    if request.method == 'GET':
        # We pass the subject from some places --
        # like the Satchmo Contact Us link
        subject = request.GET.get('subject', None)
    else:
        subject = None

    if request.method == 'POST':
        contact_form = ContactForm(request.POST)
        if contact_form.is_valid():
            sender = contact_form.cleaned_data['sender']
            subject = contact_form.cleaned_data['subject']
            cc_myself = contact_form.cleaned_data['cc_myself']
            message = contact_form.cleaned_data['message']

            if not sender:
                sender = "phylo@phylogenomics.berkeley.edu"

            recipients = ['phylo@phylogenomics.berkeley.edu']
            if cc_myself:
                recipients.append(sender)

            send_mail(subject, message, sender, recipients, fail_silently=False)
            return HttpResponseRedirect('/phylofacts/')
    else:
        contact_form = ContactForm()

    return render_to_response('common/contact_us.html', {
        'contact_form': contact_form,
        'subject': subject,
    })

def clean_webserver_temp(request):
    """Remove webserver temp files/directories owned by Apache

    The Ohana cluster has Ohana as the master node. However, a second Master
    Node, Makana, runs the webserver (under the userid of apache). Makana
    cannot directly submit to the queue and has other limitations that we have
    not been able to get changed (even with multiple requests).

    Therefore, any file that Apache creates must be removed. And, the only way
    to do that is with this script. There is no harm in having others run this
    script - it only assists us in enforcing that directories older than the
    set number of days are deleted. So, only a psuedo-security method is
    implemented. A seemingly random string is submitted as a token or the
    script will not run. However, that seemingly random string is only a
    static random string with three changing points based upon today's date.

    This can easily be scripted so that this process runs daily - and the
    token will appear to change daily.

    This is invoked by calling URL:
    http://...edu/staff/clean_webserver_temp?token=<token>

    Where token is generated from:
        now_ = datetime.datetime.now()
        token = "%dlwx3%d0K91vmL%di" % (now_.day, now_.month, now_.year)
    """
    DAYS = 7 # The number of days to retain a directory in temp
    TEMP_PATH = '/clusterfs/ohana/software/webserver/temp'

    submitted_token = request.GET.get('token', None)
    
    # Generate seemly random (but not really) token based upont today's date
    # This can be used from cron scripts to ensure there's at least some
    # knowledge of what this routine does.
    now_ = datetime.datetime.now()
    token = "%dlwx3%d0K91vmL%di" % (now_.day, now_.month, now_.year)

    if submitted_token != token:
        return HttpResponse('Invalid token or no token given. '
                            'No cleaning done.')

    os.chdir(TEMP_PATH)
    response = []
    for dir in glob.glob('*'):
        # Get timestamp from directory
        timestamp = os.stat(os.path.abspath(dir))[stat.ST_CTIME]
        file_date = datetime.datetime.fromtimestamp(timestamp)
        if file_date < now_ - datetime.timedelta(days=DAYS):
            response.append("Deleting %s.<br />\n" % dir)
            try:
                shutil.rmtree(os.path.abspath(dir))
            except Exception, e:
                response.append("%s <br />\n" % e)
        else:
            response.append("Skipping %s.<br />\n" % dir)

    return HttpResponse("".join(response))

