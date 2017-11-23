
# coding: utf-8

# In[4]:

#This notebook will construct the relative strenghts list of the DOW 30, define below 
import sys
import numpy as np
import datetime
import os


import Bio
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO 

import time

import csv

import requests

def write_results(ids, seqs, blast_PIDS, blast_HSPS, blast_uniprot):
    results = ''
    for i in range(len(blast_PIDS)):
        for j in range(len(blast_PIDS[i])):
            for k in range(len(blast_uniprot[i][j])):
                
                results += ids[i]+'\t'
                results += seqs[i][1:10] + '\t'
                results += blast_PIDS[i][j] + '\t'
                results += blast_HSPS[i][j] + '\t'
                results += blast_uniprot[i][j][k]['GO'] + '\t'
                results += blast_uniprot[i][j][k]['Function'] + '\t'
                results += blast_uniprot[i][j][k]['Database'] + '\t'
                results += '\n'
    return(results)

def results_to_listdict(ids, seqs, blast_PIDS, blast_HSPS, blast_uniprot):
    listdict = []
    for i in range(len(blast_PIDS)):
        for j in range(len(blast_PIDS[i])):
            for k in range(len(blast_uniprot[i][j])):
                listdict.append({'ID':ids[i],
                                 'Sequence':seqs[i],
                                 'Uniprot': blast_PIDS[i][j],
                                 'E-value': blast_HSPS[i][j],
                                 'GO': blast_uniprot[i][j][k]['GO'],
                                 'Function': blast_uniprot[i][j][k]['Function'],
                                 'Database': blast_uniprot[i][j][k]['Database']})
    return(listdict)
                                                
    
def blast_text(text_box):
    #Write the FASTA file out 
    fasta_file = filedialog.askopenfilename()
    #Update the textbox
    text_box.delete(1.0, tk.END)
    text_box.insert(tk.END, 'BLAST is working...Patience is a virtue')
    xml_name, ids, seqs = blast(fasta_file, text_box)
    blast_records = process_blast(xml_name)
    #list of lists of PIDS
    blast_PIDS = []
    blast_HSPS = []
    for record in blast_records:
        temp = parse_PIDS(record)
        blast_PIDS.append(temp[0])
        blast_HSPS.append(temp[1])
    text_box.delete(1.0, tk.END)
    text_box.insert(tk.END, 'Parsed BLAST')
    blast_uniprot = []
    for blast_record in blast_PIDS:
        temp_gos = []
        for PID in blast_record:
            text_box.delete(1.0, tk.END)
            text_box.insert(tk.END, 'Getting uniprot {}'.format(PID))
            temp_gos.append(parse_Uniprot(get_Uniprot(PID)))
        blast_uniprot.append(temp_gos)
    #write out the results
    text_box.delete(1.0, tk.END)
    text_box.insert(tk.END, write_results(ids, seqs, blast_PIDS, blast_HSPS,blast_uniprot))
    #write to csv
    listdict = results_to_listdict(ids, seqs, blast_PIDS, blast_HSPS,blast_uniprot)
    with open(fasta_file.split('.')[0]+'_2go.csv','w', newline='') as csvfile:
        fieldnames = ['ID', 'Sequence', 'Uniprot', 'E-value', 'GO', 'Function', 'Database']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for ld in listdict:
            writer.writerow(ld)    
    
    return(fasta_file)

def blast(file_name, text_box):
    records = SeqIO.parse(file_name, format='fasta')
    xml_name = time.strftime("%Y%m%d-%H%M%S")+'.xml'
    xml_handle = open(xml_name, 'a')
    ids = []
    seqs = []
    for record in records:
        #update the blast text prepend it to text 
        current_text = text_box.get('1.0',tk.END)
        new_text = 'BLAST...{}\n'.format(record.id)+current_text
        text_box.delete(1.0, tk.END)
        text_box.insert(tk.END, new_text)
        #get the info 
        ids.append(record.id)
        seqs.append(str(record.seq))
        result_handle = NCBIWWW.qblast('blastp', 'swissprot', record.format('fasta'))
        xml_handle.write(result_handle.read())
        result_handle.close()
    xml_handle.close()
    return(xml_name, ids, seqs)

def process_blast(xml_name):
    xml_file = open(xml_name)
    blast_records = NCBIXML.parse(xml_file)
    #can only read this once so we save it in a list 
    blast_records = list(blast_records)
    #list of blast records
    return(blast_records)

def parse_PIDS(blast_record):
    #update the textbox 
    PIDS = []
    HSPS = []
    for alignment in blast_record.alignments:
        hsp_temp = 10000
        for hsp in alignment.hsps:
            hsp_temp = min(hsp.expect, hsp_temp)
        HSPS.append(str(hsp_temp))
        title_split = alignment.title.split('|')
        PIDS.append(title_split[title_split.index('sp')+1].split('.')[0])
    return(PIDS, HSPS)

def get_Uniprot(PID):
    return(requests.get('http://www.uniprot.org/uniprot/{}.txt'.format(PID)).text)

#return a list of [{'GO':GO id, 'Database': database, 'Function': function}]
def parse_Uniprot(text):
    lines = text.split('\n')
    go_ids = []
    for line in lines:
        if line.startswith('DR   GO'):
            go_line = line.split(';')[1:] #strip off the junk identifier 
            go_dict = {}
            go_dict['GO'] = go_line[0]
            go_dict['Function'] = go_line[1] if 1<len(go_line) else 'No Function Identified'
            go_dict['Database'] = go_line[2] if 2<len(go_line) else'No Database Identified'
            go_ids.append(go_dict)
    return(go_ids)

def records_to_go(blast_records):
    go_ids = []
    for record in blast_records:
        temp_gos = []
        for pid in parse_PIDS(record):
            temp_gos.append(parse_Uniprot(get_Uniprot(pid)))
        go_ids.append(temp_gos)
    return(go_ids)

def file_to_go(file_name):
    return(records_to_go(process_blast(blast(file_name)[0])))



import tkinter as tk
from tkinter import ttk, filedialog
LARGE_FONT= ("Verdana", 12)

class Swiss2Goapp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        
        tk.Tk.wm_title(self, "SWISS2GO")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}



        frame = StartPage(container, self)

        self.frames[StartPage] = frame

        frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()

    
def clear_text(textbox):
    textbox.delete(1.0,tk.END)



class StartPage(tk.Frame): 
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label_text = tk.StringVar()
        label_text.set('Blastp -> GO')
        label = tk.Label(self, textvariable=label_text, font=LARGE_FONT) 
        label.pack(pady=10,padx=10)

        text = tk.Text(self,height='30', width='120')
        text.insert(tk.END, 'Please choose your file then be PATIENT')
        text.pack()
        
        label_date = tk.Label(self, text='Choose your FASTA file',font=LARGE_FONT)
        label_date.pack()
        

        button_update = ttk.Button(self, text="Pick file",
                            command=lambda: blast_text(text))
        button_update.pack()
        button_export = ttk.Button(self, text='Clear window',command = lambda: clear_text(text))
        button_export.pack()
        
        

app = Swiss2Goapp()
app.mainloop()





