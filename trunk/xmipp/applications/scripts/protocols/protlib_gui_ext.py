from Tkinter import *

class MultiListbox(PanedWindow):
    def __init__(self,master,lists):
        PanedWindow.__init__(self,master,borderwidth=1,showhandle=False,sashwidth=2,sashpad=0,relief=SUNKEN)
        self.lists = []
        self.columns=[]
        for l,w in lists:
            self.columns.append(l)
            frame = Frame(self); frame.pack(side=LEFT, expand=YES, fill=BOTH)
            tl=Label(frame, text=l, borderwidth=2, relief=GROOVE)
            tl.pack(fill=X)
            tl.bind('<Button-1>',self.clickon)
            lb = Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,relief=FLAT, exportselection=FALSE, selectmode=MULTIPLE ,bg='white')
            lb.pack(expand=YES, fill=BOTH)
            self.lists.append(lb)
            lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
            lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
            lb.bind('<Leave>', lambda e: 'break')
            lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
            lb.bind('<Button-2>', lambda e, s=self: s._button2(e.x, e.y))
            lb.bind('<Button-4>', lambda e, s=self: s._scroll(SCROLL, 1, PAGES))
            lb.bind('<Button-5>', lambda e, s=self: s._scroll(SCROLL, -1, PAGES))
            self.add(frame)
        Label(master, borderwidth=1, relief=FLAT).pack(fill=X)
        sb = Scrollbar(master, orient=VERTICAL, command=self._scroll,borderwidth=1)
        sb.pack(fill=Y,side=RIGHT,expand=NO)
        for l in self.lists:
            l['yscrollcommand']=sb.set
        self.add(frame)
        self.pack(expand=YES,fill=BOTH)
        self.sortedBy=-1
        self.previousWheel=0


    def _select(self, y,state=16):
        row = self.lists[0].nearest(y)
        if state==16:self.selection_clear(0, END)
        self.selection_set(row)
##        print self.curselection()
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
        return 'break'


    def clickon(self,e):
        self._sortBy(self.columns.index(e.widget['text']))


    def _sortBy(self, column):
        """ Sort by a given column. """
        if column == self.sortedBy:
            direction = -1 * self.direction
        else:
            direction = 1

        elements = self.get(0, END)
        self.delete(0, END)
        elements.sort(lambda x, y: self._sortAssist(column, direction, x, y))
        self.insert(END, *elements)

        self.sortedBy = column
        self.direction = direction


    def _sortAssist(self, column, direction, x, y):
        c = cmp(x[column], y[column])
        if c:
            return direction * c
        else:
            return direction * cmp(x, y)

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
        #print self.curselection()


if __name__ == '__main__':
    tk = Tk()
    Label(tk, text='MultiListbox').pack(side=TOP)
    mlb = MultiListbox(tk, (('Subject', 40), ('Sender', 20), ('Date', 10)))
    for i in range(5000):
        mlb.insert(END, ('Important Message: %d' % i, 'John Doe %d' % i, '10/10/%04d' % (1900+i)))
    mlb.pack(expand=YES,fill=BOTH,side=TOP)
    tk.mainloop()

