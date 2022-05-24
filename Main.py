from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog


class GUI:

    def __init__(self):
        self.frame = tk.Frame()
        self.frame.grid()

        self.filepath = None
        self.genome = None
        self.trie = None

        self.label_genome = ttk.Label(self.frame, text="Load in a genome:", style='text.TLabel')
        self.label_genome.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W)

        self.open_file_dialog_button = ttk.Button(self.frame, text='Choose a file', command=self.file_dialog,
                                                  style='my.TButton')
        self.open_file_dialog_button.grid(row=1, column=0, padx=10, pady=5, ipadx=10, ipady=5)

        self.file_name = ttk.Label(self.frame, text="", style="text.TLabel")
        self.file_name.grid(row=1, column=1, padx=10, pady=10, sticky=tk.W)

        self.temp_label = ttk.Label(self.frame, text="Enter annealing temperature:")
        self.temp_label.grid(row=2, column=0, padx=10, pady=10, sticky=tk.W)
        self.enter_temp = ttk.Entry(self.frame, width=35)
        self.enter_temp.grid(row=2, column=1, padx=10, pady=10, sticky=tk.W)

        self.deltaT_label = ttk.Label(self.frame, text="Enter deltaT: \n(lower values = longer runtime):")
        self.deltaT_label.grid(row=3, column=0, padx=10, pady=10, sticky=tk.W)
        self.enter_deltaT = ttk.Entry(self.frame, width=35)
        self.enter_deltaT.grid(row=3, column=1, padx=10, pady=10, sticky=tk.W)

        self.filler1 = ttk.Label(self.frame, text="", style="text.TLabel")
        self.filler1.grid(row=4, column=0, padx=10, pady=10, sticky=tk.W)

        self.run = ttk.Button(self.frame, text='Run the program', command=self.file_dialog,
                                                  style='my.TButton')
        self.run.grid(row=5, column=0, padx=10, pady=5, ipadx=10, ipady=5)

    def file_dialog(self):
        self.filepath = filedialog.askopenfilename(filetypes=(('Fasta files', '*.fasta'),))
        self.file_name["text"] = self.filepath

        if self.filepath:
            print(self.filepath)
            try:
                self.genome = Sequence(self.filepath)
            except:
                tk.messagebox.showerror("Something went wrong")









def run_gui():

    root = tk.Tk()
    root.geometry("700x500")
    root.title("PCR-primers GUI")

    GUI()

    root.mainloop()


if __name__ == "__main__":
    run_gui()
