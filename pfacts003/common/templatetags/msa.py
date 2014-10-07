from django.template.defaultfilters import stringfilter
from django.utils.html import escape
from django.utils.safestring import mark_safe
from django.utils.html import conditional_escape
from django import template

import re

register = template.Library()


# WARNING: presently this doesn't error-check; use at your own risk
@register.filter
@stringfilter
def parallel_msa(value, arg="20", autoescape=None):
    # Check the input
    if autoescape:
        esc = lambda x: x
        splitter = '>'
    else:
        esc = conditional_escape
        splitter = '&gt;'
    
    size = int(arg)

    return mark_safe('\n'.join((
        # Format the string with left-justified header then the seq
        ('<span title="%%s"><span>&gt;%%-%i.*s</span> %%s</span>' % size) % (
            # Provide title text for the full header
            esc(header),
            # Print at most 20 characters of the header
            size, header,
            # Linearize the seq
            seq.replace('\n',''),
        ) for header, seq in (
            # Break the text into header and seq
            entry.split('\n',1) for entry in value.split(splitter)[1:]
        )
    )))
parallel_msa.needs_autoescape = True
