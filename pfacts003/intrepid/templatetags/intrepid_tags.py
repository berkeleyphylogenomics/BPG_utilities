from django import template

register = template.Library()

@register.filter(name='readable_percent')
def readable_percent(value, d):
    """ Multiplies a float by 100. """
    return "%s %%" % (str(round(100.0*float(value), int(d)))) 
