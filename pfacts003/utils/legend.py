import copy
import colorsys
from django.utils import simplejson as json
from django.utils.html import escape

class Legend(object):
    
    def __init__(self, headers, links=None, constants=None, color=True, entries=None, invisible='style="display: none;"'):
        self.headers = headers
        self.links = links
        self.constants = constants
        self.color = color
        self.invisible = ' ' + invisible

        if entries:
            self.entries = entries
        else:
            self.entries = []

    def __len__(self):
        return len(self.rows)

    def __unicode__(self):

        new_entries = copy.deepcopy(self.entries)
        length = float(len(new_entries))

        # print the headers
        leg = '<tr>'
        if self.color:
            leg += '<th></th>'
        leg += '<th>%s</th></tr>\n'%'</th><th>'.join(self.headers)

        for i,entry in enumerate(new_entries):
            # assign the color (and print)
            leg += '<tr>'
            if self.color:
                if not entry.has_key('_color'):
                    entry['_color'] = '#' + ''.join('%02x'%(ele*255) for ele in colorsys.hsv_to_rgb(i/length, 0.9, 0.75))
                leg += '<td bgcolor="%s">&nbsp;&nbsp;</td>'%entry['_color']
            # assign the link(s)
            if self.links:
                if isinstance(self.links,dict):
                    entry['_links'] = dict((header,self.links[header][0]%(
                        entry[item] for item in self.links[header][1]
                    )) for header in self.links)
                else:
                    entry['_link'] = self.links[0]%(
                        entry[item] for item in self.links[1]
                    )
            # print the entries
            for header in self.headers:
                cell = '<td>%s</td>'
                if entry.has_key('_links') and entry['_links'].has_key(header):
                    cell = cell%'<a href="%s">%%s</a>' % \
                        escape(entry['_links'][header])
                leg += cell%entry[header]
            leg += '</tr>\n'

        # print hidden JSON data
        leg += '<tr%s><td>%s</td></tr>\n'%(self.invisible, escape(json.dumps({
            'constants': self.constants,
            'entries': new_entries,
        })))

        return u'%s' % leg

    def __str__(self):
        return self.__unicode__()

    def append(self, input):
        self.entries.append(input)
