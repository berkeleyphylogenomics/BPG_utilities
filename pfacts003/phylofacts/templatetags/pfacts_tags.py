from django import template

register = template.Library()

@register.filter
def escapenl(value):
    return value.replace('\n', '%0D%0A')

