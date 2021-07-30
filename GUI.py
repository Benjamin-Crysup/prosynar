
import tkinter
import subprocess

import prosynargui
import prosynarFpprobreggui
import prosynarFprovergui
import prosynarFrovergui
import prosynarMaligngui
import prosynarMflashgui
import prosynarMpeargui

allArgs = ["prosynar"]

# main gui
root = tkinter.Tk()
root.title("prosynar GUI Helper")
app = prosynargui.ArgumentGUIFrame(master = root)
app.mainloop()
curArgs = app.getArgs()
if curArgs is None:
    system.exit()
allArgs.extend(curArgs)

#figure out which filters and which merger to run
class FilterMergeGUIFrame(tkinter.Frame):
    def __init__(self,master=None):
        super().__init__(master)
        self.wasOkP = False
        self.master = master
        canv = tkinter.Canvas(self)
        self.filtLab = tkinter.Label(canv, text='Filters')
        self.filtLab.grid(column = 0, row = 0)
        self.varFrover = tkinter.IntVar()
        self.chkFrover = tkinter.Checkbutton(canv, text='Naive Reference Overlap', variable=self.varFrover)
        self.chkFrover.grid(column = 0, row = 1)
        self.varFprover = tkinter.IntVar()
        self.chkFprover = tkinter.Checkbutton(canv, text='Probabilistic Reference Overlap', variable=self.varFprover)
        self.chkFprover.grid(column = 0, row = 2)
        self.varFpprobreg = tkinter.IntVar()
        self.chkFpprobreg = tkinter.Checkbutton(canv, text='Probabilistic Reference Repeat', variable=self.varFpprobreg)
        self.chkFpprobreg.grid(column = 0, row = 2)
        mergeSRow = 3
        self.filtLab = tkinter.Label(canv, text='Merger')
        self.filtLab.grid(column = 0, row = mergeSRow + 0)
        self.varMerge = tkinter.StringVar()
        self.varMerge.set('Malign')
        self.radMalign = tkinter.Radiobutton(canv, text='Align', value='Malign', variable = self.varMerge)
        self.radMalign.grid(column = 0, row = mergeSRow + 1)
        self.radMflash = tkinter.Radiobutton(canv, text='FLASH', value='Mflash', variable = self.varMerge)
        self.radMflash.grid(column = 0, row = mergeSRow + 2)
        self.radMpear = tkinter.Radiobutton(canv, text='PEAR', value='Mpear', variable = self.varMerge)
        self.radMpear.grid(column = 0, row = mergeSRow + 3)
        endRow = mergeSRow + 4
        self.canBtn = tkinter.Button(canv, text = 'Cancel', command = lambda: self.goCancel())
        self.canBtn.grid(column = 0, row = endRow)
        self.goBtn = tkinter.Button(canv, text = 'OK', command = lambda: self.goGo())
        self.goBtn.grid(column = 1, row = endRow)
        canv.pack(side=tkinter.LEFT)
        mainscr = tkinter.Scrollbar(self, command=canv.yview)
        mainscr.pack(side=tkinter.LEFT, fill='y')
        canv.configure(yscrollcommand = mainscr.set)
        self.pack()
    def goCancel(self):
        self.master.destroy()
    def goGo(self):
        self.wasOkP = True
        self.master.destroy()
root = tkinter.Tk()
root.title("prosynar Filter/Merge Select")
app = FilterMergeGUIFrame(master = root)
app.mainloop()
if not app.wasOkP:
    system.exit()
hotApps = []
hotAppNames = []
if app.varFrover.get() != 0:
    hotApps.append(prosynarFrovergui.ArgumentGUIFrame)
    hotAppNames.append("prosynar Reference Overlap Filter")
if app.varFprover.get() != 0:
    hotApps.append(prosynarFprovergui.ArgumentGUIFrame)
    hotAppNames.append("prosynar Probabilistic Reference Overlap Filter")
if app.varFpprobreg.get() != 0:
    hotApps.append(prosynarFpprobreggui.ArgumentGUIFrame)
    hotAppNames.append("prosynar Probabilistic Repeat Region Filter")
if app.varMerge.get() == "Malign":
    hotApps.append(prosynarMaligngui.ArgumentGUIFrame)
    hotAppNames.append("prosynar Alignment Merger")
elif app.varMerge.get() == "Mflash":
    hotApps.append(prosynarMflashgui.ArgumentGUIFrame)
    hotAppNames.append("prosynar FLASH Merger")
elif app.varMerge.get() == "Mpear":
    hotApps.append(prosynarMpeargui.ArgumentGUIFrame)
    hotAppNames.append("prosynar PEAR Merger")

# get arguments for the pieces
for i in range(len(hotApps)):
    root = tkinter.Tk()
    root.title(hotAppNames[i])
    app = hotApps[i](master = root)
    app.mainloop()
    appArgs = app.getArgs()
    if appArgs is None:
        system.exit()
    allArgs.append("--")
    allArgs.extend(appArgs)

# run the stupid thing
subprocess.run(allArgs)

