from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, Listbox, Scrollbar
import tkinter.scrolledtext as st


class GUI:

    def __init__(self, root):
        self.frame = tk.Frame()
        self.frame.grid()
        self.root = root

        self.filepath = None
        self.genome = None

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
        self.text_area = st.ScrolledText(self.frame, width=65, height=20)
        self.text_area.grid(row=9, column=0, padx=10, pady=5, columnspan=3, sticky=tk.W)

        self.display_primers_button = ttk.Button(self.frame, text="View primers", command=self.call_display_primers,
                                                 style="my.TButton")
        self.display_primers_button.grid(row=9, column=4, padx=20, pady=5, ipadx=10, ipady=5)

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
            delta_T = int(self.enter_deltaT.get())
            temp = int(self.enter_temp.get())
            self.text_area.insert(tk.INSERT, f"Delta_T = {delta_T}\nTemp = {temp}\n")
            self.text_area.insert(tk.INSERT, "Building trie.\n")

            trie = self.genome.build_trie(20)
            self.text_area.insert(tk.INSERT, "Building of trie is done.\n")
            self.text_area.insert(tk.INSERT, "Generating potential primers.\n")

            primers = Primers(self.genome, 20)
            self.text_area.insert(tk.INSERT, "Potential primers generated, removal of bad primers beginning.\n")

            primers.GC_clamp()
            primers.GC_content(60, 40)
            primers.temp_selection(temp)

            self.text_area.insert(tk.INSERT, "Only acceptable primers remain.\n")
            self.text_area.insert(tk.INSERT, "Running the remaining primers through the trie.\n")

            forward_primers = search(trie, primers.get_frw_primers(), delta_T)
            reverse_primers = search(trie, primers.get_frw_primers(), delta_T)
            print(len(forward_primers))
            print(len(reverse_primers))

            self.text_area.insert(tk.INSERT, "Done.\n")
            self.text_area.insert(tk.INSERT, "Generating primer pairs.\n")

            self.primer_list, self.circular_list = sort_primers(forward_primers, reverse_primers, self.genome)
            print(len(self.primer_list))
            print(len(self.circular_list))
            self.text_area.insert(tk.INSERT, "Primer pairs generated, press ""View Primers"".\n")

    def call_display_primers(self):
        self.window = tk.Toplevel()
        self.frame_list = tk.Frame(self.window)
        self.frame_list.grid(row=0, column=0)

        DisplayPrimers(parent=self.frame_list, primer_list=self.primer_list, circular_list=self.circular_list, root=self.root)


class DisplayPrimers:
    def __init__(self, parent, primer_list, circular_list, root):
        self.parent = parent
        self.frame = ttk.Frame(self.parent)
        self.frame.pack()
        self.primer_list = primer_list
        self.circular_list = circular_list
        self.root = root

        self.title_label = ttk.Label(self.frame, text="Non-circular primer pairs:")
        self.title_label.grid(row=0, column=0, sticky=tk.W)

        self.box_label = ttk.Label(self.frame, text="forward primer:\t|  Reverse primer:\t| Start:\t| Stop:\t"
                                                    "| Length:\t  | EcoRI digestion fragments:")
        self.box_label.grid(row=1, column=0, sticky=tk.W)

        self.list_box = Listbox(self.frame, width=90, height=20)
        self.list_box.grid(row=2, column=0, padx=10, pady=5)

        self.scrollbar = Scrollbar(self.frame, orient="vertical")
        self.scrollbar.grid(row=2, column=1, sticky=tk.NS)

        self.list_box.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.list_box.yview)

        self.title_label = ttk.Label(self.frame, text="Circular primer pairs:")
        self.title_label.grid(row=3, column=0, sticky=tk.W)

        self.box_label2 = ttk.Label(self.frame, text="forward primer:\t|  Reverse primer:\t| Start:\t| Stop:\t"
                                                    "| Length:\t  | EcoRI digestion fragments:")
        self.box_label2.grid(row=4, column=0, sticky=tk.W)

        self.list_box2 = Listbox(self.frame, width=90, height=20)
        self.list_box2.grid(row=5, column=0, padx=10, pady=5)

        self.scrollbar2 = Scrollbar(self.frame, orient="vertical")
        self.scrollbar2.grid(row=5, column=1, sticky=tk.NS)

        self.list_box2.config(yscrollcommand=self.scrollbar2.set)
        self.scrollbar2.config(command=self.list_box2.yview)

        self.copy_primer_button = ttk.Button(self.frame, text="Select", command=self.copy_button,
                                             style='my.TButton')
        self.copy_primer_button.grid(row=4, column=1, padx=10, pady=5)

        # insert items into the list boxes
        for item in self.primer_list:
            self.list_box.insert("end", item)
        for item in self.circular_list:
            self.list_box2.insert("end", item)

    def copy_button(self):
        try:
            primer_info = self.list_box.get("anchor")
            print(primer_info)
            primers = primer_info[:2]
            self.root.clipboard_clear()
            self.root.clipboard_append(primers)
        except:
            tk.messagebox.showerror("Error")


def run_gui():

    root = tk.Tk()
    root.geometry("700x700")
    root.title("PCR-primers GUI")

    GUI(root)

    root.mainloop()


if __name__ == "__main__":
    run_gui()
