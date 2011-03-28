#!/usr/bin/env python
"""/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
"""
import os, glob, sys, optparse
from Tkinter import *

class MultiListbox(Frame):
    def __init__(self, master, lists):
       Frame.__init__(self, master)
       self.lists = []
       for l, w in lists:
           frame = Frame(self); 
           frame.pack(side='left', expand=YES, fill=BOTH)
           Label(frame, text=l, borderwidth=1, relief=RAISED).pack(fill=X)
           lb = Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
              relief=FLAT, exportselection=FALSE)
           lb.pack(expand=YES, fill=BOTH)
           self.lists.append(lb)
           lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
           lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
           lb.bind('<Leave>', lambda e: 'break')
           lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
           lb.bind('<Button-2>', lambda e, s=self: s._button2(e.x, e.y))
       frame = Frame(self); frame.pack(side='left', fill=Y)
       Label(frame, borderwidth=1, relief=RAISED).pack(fill=X)
       sb = Scrollbar(frame, orient=VERTICAL, command=self._scroll)
       sb.pack(expand=YES, fill=Y)
       self.lists[0]['yscrollcommand']=sb.set

    def _select(self, y):
       row = self.lists[0].nearest(y)
       self.selection_clear(0, END)
       self.selection_set(row)
       return 'break'

    def _button2(self, x, y):
       for l in self.lists: l.scan_mark(x, y)
       return 'break'

    def _b2motion(self, x, y):
       for l in self.lists: l.scan_dragto(x, y)
       return 'break'

    def _scroll(self, *args):
       for l in self.lists:
           apply(l.yview, args)

    def curselection(self):
        return self.lists[0].curselection()

    def delete(self, first, last=None):
       for l in self.lists:
           l.delete(first, last)

    def get(self, first, last=None):
       result = []
       for l in self.lists:
           result.append(l.get(first,last))
       if last: return apply(map, [None] + result)
       return result
      
    def index(self, index):
        self.lists[0].index(index)

    def insert(self, index, *elements):
       for e in elements:
           i = 0
           for l in self.lists:
              l.insert(index, e[i])
              i = i + 1

    def size(self):
        return self.lists[0].size()

    def see(self, index):
        for l in self.lists:
            l.see(index)

    def selection_anchor(self, index):
       for l in self.lists:
           l.selection_anchor(index)

    def selection_clear(self, first, last=None):
       for l in self.lists:
           l.selection_clear(first, last)

    def selection_includes(self, index):
        return self.lists[0].selection_includes(index)

    def selection_set(self, first, last=None):
       for l in self.lists:
           l.selection_set(first, last)
           
   
def command_line_options():
      """ add command line options here"""
      _usage="""usage: %prog [options]
Example:
   xmipp_metadata_show -i all_images.sel"""
      parser = optparse.OptionParser(_usage)
      parser.add_option("-i", "--input",  dest="input", type="string",
            help="The input metadata")        
      (options, args) = parser.parse_args()

      return options.input

scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path

from xmipp import *

if __name__ == '__main__':
    input = command_line_options()
    md = MetaData(input)
    tk = Tk()
    Label(tk, text='Metadata: ' + input).pack()
    labels = md.getActiveLabels()
    columns = []
    # Create labels columns
    for l in labels:
        max = md.getMaxStringLength(l)
        columns.append((label2Str(l), md.getMaxStringLength(l) + 2))
    # Create grid
    mlb = MultiListbox(tk, columns)
    for id in md:
        row = []
        for l in labels:
            row.append(str(md.getValue(l, id)))            
        mlb.insert(END, row)
        
    mlb.pack(expand=YES,fill=BOTH)
    tk.mainloop()