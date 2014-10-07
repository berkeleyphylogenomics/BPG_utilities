import pickle
import subprocess
import re
import os
from types import FunctionType
from django.conf import settings

ERROR, STARTED, COMPLETED = range(-1,2)

_argset_keys = ('kw__','kw_','fl__','fl_')
_word = r'[a-zA-Z][a-zA-Z0-10]*'

# matches fields associated with queued app steps
field_name_patt = r'^(?P<step_label>%s)__(?P<field_key>.*)$'

# matches the <field_key> portion of args associated with queued app steps
arg_name_re = re.compile(
    r'^(?P<arg_type>%s)(?P<arg_label>%s)$' \
    % ('|'.join(_argset_keys), _word)
)

# matches the full name of args associated with queued app steps
#field_arg_name_re = re.compile(
#    r'^(?P<step_label>%s)__(?P<arg_type>%s)(?P<arg_label>%s)$' \
#    % (_word, '|'.join(_argset_keys), _word)
#)

def _q_import(app_name):
    """Dynamically import queued from app using pfacts003.app_name.queued structure"""

    # only import allowed queued apps
    if app_name not in settings.QUEUED_APPS:
        raise ImportError

    app = __import__('pfacts003.' + app_name + '.queued')
    app = getattr(app, app_name).queued

    return app

class Q_data(object):
    def __init__(self, cleaned_data, app_name):
        
        step_labels = [step.label for step in _q_get_steps(_q_import(app_name))]

        self.cleaned_data = cleaned_data

        steps_data = dict([(label, {'fields': {}, 'args': []}) \
            for label in step_labels])

        for key, val in cleaned_data.items():

            # parse only if the field name matches our pattern
            field_name_match = re.match(
                field_name_patt % '|'.join(step_labels), key)
            if field_name_match:
                step_label, field_key = field_name_match.groups()
            
                this_step = steps_data[step_label]
            
                # parse as an argument only if the field name matches our
                # pattern
                arg_name_match = arg_name_re.match(field_key)
                if arg_name_match:
                    arg_type, arg_label = arg_name_match.groups()
                    new_key = '-' + arg_label
                    # add second hyphen, if necessary
                    if field_key.endswith('__'):
                        new_key += '-' + new_key
                    this_step['fields'][new_key] = val
                    # add the key, if the field isn't a false flag
                    if not (arg_type.startswith('fl') and not val):
                        this_step['args'].append(new_key)
                    # add the value if the field is a keyword
                    if arg_type.startswith('kw'):
                        this_step['args'].append(val)
                # otherwise, parse as a field (not necessarily for the command
                # line)
                else:
                    this_step['fields'][field_name_match.group('field_key')] = val
        
        self._parsed_input = steps_data
        
    def _set_step(self, step_label):
        step_data = self._parsed_input[step_label]
        self.fields = step_data['fields']
        self.args = step_data['args']
        if hasattr(self,'returncode'):
            del self.returncode

    def _set_returncode(self, returncode):
        self.returncode = returncode

    def _set_logger(self, logger):
        self.logger = logger

class q_decorated_function(object):
    def __init__(self, f, count, title=None, percent=None, command_str=None, \
        command_out=None, command_args=None, display_outfile=None, \
        write_field=None):
        self.label = f.__name__
        self.title = title
        self.percent = percent
        self.command_str = command_str
        self.command_out = command_out
        self.command_args = command_args
        self.display_outfile = display_outfile
        self.write_field = write_field
        self.f = f
        self.count = count
    def __call__(self, q_data, *args, **kwargs):

        # write a field to a file, if requested
        if (isinstance(self.write_field, tuple) or \
            isinstance(self.write_field, list)) and len(self.write_field) == 2:

            src, dst = self.write_field
            # make sure the field exists
            if not q_data.cleaned_data.has_key(src):
                # field apparently does not exist
                raise Exception


            # write the data to a file
            handle = open(dst, 'w')
            handle.write(q_data.cleaned_data[src])
            handle.close()
            
        # execute a command-line utility, if requested
        if self.command_str:
            # write to an outfile
            handle = open(
                # by default, use <step label>.out...
                self.command_out is None and ('%s.out' % self.label) or \
                # ...otherwise use command_out as the filename
                self.command_out , 'w')
            handle_e = open('%s.err' % self.label, 'w')
            q_data._set_returncode(subprocess.call(
                # combine the arguments
                [self.command_str] \
                    + q_data.args \
                    + self.command_args,
                # direct output
                stdout=handle,
                stderr=handle_e))
            handle.close()
            handle_e.close()
        
        # return f's return value
        return self.f(q_data, *args, **kwargs)

class _q_function(object):
    count = 0

    def __init__(self, title=None, percent=None, command_str=None, \
        command_out=None, command_args=None, display_outfile=None, \
        write_field=None):

        self.title = title
        self.percent = percent
        self.command_str = command_str
        self.command_out = command_out
        if command_args is None:
            self.command_args = []
        else:
            self.command_args = command_args
        self.display_outfile = display_outfile
        self.write_field = write_field

        self.count = _q_function.count
        _q_function.count += 1

    def __call__(self, f):
        if f.__name__.count('__'):
            # q_functions cannot contain the '__' substring
            raise Exception
        if isinstance(self.title,basestring):
            title = self.title
        else:
            title = f.__name__
        return q_decorated_function(f, self.count, title, self.percent, \
            self.command_str, self.command_out, self.command_args, \
            self.display_outfile, self.write_field)

def q_function(title=None, percent=None, command_str=None, command_out=None,
    command_args=None, display_outfile=None, write_field=None):
    if isinstance(title, FunctionType):
        return _q_function()(title)
    else:
        return _q_function(title, percent, command_str, command_out,
        command_args, display_outfile, write_field)

def q_job_is_finished(app, work_path):
    return JobInfo(app, work_path).is_finished()

def q_read_log(app, work_path):
    return JobInfo(app, work_path)

def _q_get_steps(app):
    return sorted([
        getattr(app, item) \
        for item in dir(app) \
        if isinstance(getattr(app, item), q_decorated_function)
    ], key=lambda x: x.count)

class _LogInfoBase(object):

    def is_finished(self):
        return bool(self.stop)

    def get_info(self):
        return filter(lambda x: x['levelname'] == 'INFO', self.entries)

    def get_entries(self):
        bookends = ('START', 'STOP', 'STEP_START', 'STEP_COMPLETE')
        return filter(lambda x: x['levelname'] not in bookends, self.entries)

    def get_error(self):
        return filter(lambda x: x['levelname'] == 'ERROR', self.entries)

class StepInfo(_LogInfoBase):
    
    def _set_start_stop_errored(self):
        self.start = None
        self.stop = None
        self.errored = False
        for entry in self.entries:
            if entry['levelname'] == 'STEP_START':
                self.start = entry['asctime']
            elif entry['levelname'] == 'STEP_COMPLETE':
                self.stop = entry['asctime']
            elif entry['levelname'] == 'ERROR':
                self.errored = True
    
    def __init__(self, step, percent, display_outfile, entries):
        self.title = step.title
        self.label = step.label
        self.percent = percent
        self.entries = entries
        self.info = []
        self.errors = []
        self._set_start_stop_errored()
        if display_outfile is not None and self.start is not None and \
            os.path.exists(display_outfile):
            f = open(display_outfile)
            self.display_outfile = f.read()
            f.close()

class JobInfo(_LogInfoBase):
    
    def _set_start_stop_errored(self):
        self.start = None
        self.stop = None
        self.errored = False
        for entry in self.entries:
            if entry['levelname'] == 'START':
                self.start = entry['asctime']
            elif entry['levelname'] == 'STOP':
                self.stop = entry['asctime']
            elif entry['levelname'] == 'ERROR':
                self.errored = True
    
    def __init__(self, app, work_path):

        # next line is fairly hackish...
        self.label = app.__name__.split('.')[-2]
        self.title = app.title
        log_path = os.path.join(work_path, '%s.log' % self.label)

        if os.path.exists(log_path):
            log_handle = open(log_path)
            log_data = log_handle.read()
            log_handle.close()
            log_field_patts = {
                'levelname': r'(?P<levelname>START|STOP|STEP_START|STEP_COMPLETE|ERROR|INFO)',
                'name': r'(?P<name>application|step\.\w*?)',
                'message': r'<<< (?P<message>.*?) >>>',
            }

            # make a regex describing a single entry from the log
            entry_re = re.compile('^%s$' % (' - '.join([ \
                log_field_patts.get(field, '([^\n\r]*?)') \
                for field in ('asctime', 'name', 'levelname', 'message')])),
            re.M | re.S)

            # make a list of dictionaries, one for each entry in the list
            self._all_entries = [dict(zip(
                ('asctime', 'name', 'levelname', 'message'), match.groups())) \
                for match in entry_re.finditer(log_data)]

            # make a list of dictionaries for this application
            self.entries = filter(lambda x: x['name'] == 'application',
                self._all_entries)
        else:
            self._all_entries = []
            self.entries = []

        
        self._set_start_stop_errored()
        self.errored |= bool(len(filter(lambda x: x['levelname'] == 'ERROR',
            self._all_entries)))

        app_steps = _q_get_steps(app)
        def_val = (100 - sum(filter(lambda x: x.percent is not None,
            app_steps)))/float(len(app_steps))
        self.steps = [StepInfo(
            # the step
            step,
            # the percent of the step (or default percent)
            step.percent is not None and step.percent or def_val,
            # the outfile to display (with path)
            step.display_outfile is not None and \
                os.path.join(work_path, step.display_outfile) or None,
            # the log entries for the step
            filter(lambda x: x['name'] == 'step.%s' % step.label,
                self._all_entries),
        ) for step in app_steps]

    def get_percent_complete(self):
        return int(sum(map(lambda x: x.stop is not None and x.percent or 0,
            self.steps)))
