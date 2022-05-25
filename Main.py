from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, Listbox, Scrollbar
import tkinter.scrolledtext as st


class GUI:

    def __init__(self):
        self.frame = tk.Frame()
        self.frame.grid()

        self.filepath = None
        self.genome = None
        self.trie = None
        self.primer_list = None
        self.circular_list = None

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
        self.filler1 = ttk.Label(self.frame, text="", style="text.TLabel")
        self.filler1.grid(row=7, column=0, padx=10, pady=10, sticky=tk.W)

        self.run = ttk.Button(self.frame, text='Run the program', command=self.run_program,
                              style='my.TButton', state="disabled")
        self.run.grid(row=6, column=0, padx=10, pady=5, ipadx=10, ipady=5)

        self.check_button = ttk.Button(self.frame, text='Check input values', command=self.check_input,
                                       style='my.TButton')
        self.check_button.grid(row=5, column=0, padx=10, pady=5, ipadx=10, ipady=5)

        self.window_label = ttk.Label(self.frame, text="Messages will be displayed here:")
        self.window_label.grid(row=8, column=0, padx=10, pady=5, sticky=tk.W)
        self.text_area = st.ScrolledText(self.frame, width=60, height=20)
        self.text_area.grid(row=9, column=0, padx=10, pady=5, columnspan=3, sticky=tk.W)

        self.display_primers_button = ttk.Button(self.frame, text="View primers", command=self.call_display_primers,
                                                 style="my.TButton")
        self.display_primers_button.grid(row=9, column=4, padx=10, pady=5, ipadx=10, ipady=5)

    def file_dialog(self):
        self.filepath = filedialog.askopenfilename(filetypes=(('Fasta files', '*.fasta'),))
        self.file_name["text"] = self.filepath

        if self.filepath:
            try:
                self.genome = Sequence(self.filepath)
                self.text_area.insert(tk.INSERT, "Genome loaded\n")
            except:
                tk.messagebox.showerror("Error")

    def check_input(self):
        delta_T = self.enter_deltaT.get()
        temp = self.enter_temp.get()

        if self.filepath and delta_T and temp:
            self.run.config(state="normal")
            tk.messagebox.showinfo("Good values", "You can now press the run program button.")
        else:
            tk.messagebox.showinfo("Missing values", "one or more values seems to be missing\nor you did not load in a"
                                                     " fasta file")

    def run_program(self):
        try:
            delta_T = self.enter_deltaT.get()
            temp = self.enter_temp.get()
            self.text_area.insert(tk.INSERT, f"Delta_T = {delta_T}\nTemp = {temp}\n")

        except:
            tk.messagebox.showerror("Error", "Something went wrong")

    def call_display_primers(self):
        self.window = tk.Toplevel()
        self.frame_list = tk.Frame(self.window)
        self.frame_list.grid(row=0, column=0)

        DisplayPrimers(parent=self.frame_list, primer_list=self.primer_list, circular_list=self.circular_list)

        self.window.grab_set()


class DisplayPrimers:
    def __init__(self, parent, primer_list, circular_list):
        self.parent = parent
        self.frame = ttk.Frame(self.parent)
        self.frame.pack()
        self.primer_list = primer_list

        self.box_label = ttk.Label(self.frame, text="forward primer:\t  Reverse primer:\t  Start:\t  Stop:\t  Length:")
        self.box_label.grid(row=0, column=0, sticky=tk.W)

        self.list_box = Listbox(self.frame, width=80, height=40)
        self.list_box.grid(row=1, column=0, padx=10, pady=5)

        self.scrollbar = Scrollbar(self.frame, orient="vertical")
        self.scrollbar.grid(row=1, column=1, sticky=tk.NS)

        self.list_box.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.list_box.yview)


def run_gui():

    root = tk.Tk()
    root.geometry("700x700")
    root.title("PCR-primers GUI")

    GUI()

    root.mainloop()


if __name__ == "__main__":
    run_gui()
