import sys
import time

import os
from barcode import get_barcode
import gspread

req_per_min = 60 # i think
def parse_abi_name(abi_name):
    abi_name = abi_name.split(".ab1")[0]
    abi_tokens = reversed(abi_name.split("_"))
    well,primer,colony,plate,name = abi_tokens
    return name,plate,colony,primer

def parse_dirname(reads_dir):
    name = reads_dir.split("/")[-1]
    email,order_id,seq_date = name.split("_")
    return email,order_id,seq_date

def init_worksheet():
    gc = gspread.oauth()
    wks = gc.open("Preps")
    return wks

def process_folder(folder,wks):
    processed_sheet = wks.worksheet("processed_folders")
    folder_cols = processed_sheet.find("folder_name").col
    folders = processed_sheet.col_values(folder_cols)
    next_row = len(folders) + 1

    processed_sheet.update_cell(next_row,folder_cols,folder)
    reads_dir = f"{seq_dir}/{folder}"
    process_barcodes(reads_dir,wks)


def process_barcodes(reads_dir,wks,verbose=True):
    _,order_id,seq_date = parse_dirname(reads_dir)
    abis = sorted([f"{reads_dir}/{f}" for f in os.listdir(reads_dir) if f.endswith(".ab1") and not f.startswith(".")])
    items = len(abis)

    print(f"Processing {items} items...")
    print("application maybe throttled, idk how to do this ideally with rate limiting")
    print(" ¯\_(ツ)_/¯")

    sheet = wks.worksheet("sequencing")
    name_col = sheet.find("Name")
    plate_col = sheet.find("Plate")
    colony_col = sheet.find("Colony")
    primer_col = sheet.find("Primer")
    barcode_col = sheet.find("Barcode").col
    order_col = sheet.find("order_id").col
    seq_date_col = sheet.find("seq_date").col

    rows = []
    for abi in abis:
        bcode = get_barcode(abi)
        abi_name = abi.split("/")[-1]
        name,plate,colony,primer = parse_abi_name(abi_name)
        row = [name,plate,colony,primer,bcode,order_id,seq_date]
        rows.append(row)
        if verbose:
            for item in row:
                print(item)
            print("---")

    next_row = len(wks.col_values(name_col))+1
    last_row = next_row + items
    first_cell = f"{name_col}{next_row}"
    last_cell = f"{seq_date_col}{last_row}"
    cell_range = f"{first_cell}:{last_cell}"
    wks.update(cell_range,rows)

def get_unprocessed_folders(seq_dir,wks):
    seqs = [f for f in os.listdir(seq_dir) if not f.startswith(".") and not f == "old"]
    processed_sheet = wks.worksheet("processed_folders")
    folder_cols = processed_sheet.find("folder_name").col
    processed_folders = processed_sheet.col_values(folder_cols)
    unprocessed = [f for f in seqs if f not in processed_folders]
    return unprocessed
def process_folders(seq_dir,wks):
    folders = get_unprocessed_folders(seq_dir,wks)
    for f in folders:
        process_folder(f,wks)

if __name__ == "__main__":
    seq_dir = "/Users/adrian/Documents/research/sequencing"
    wks = init_worksheet()
    process_folders(seq_dir,wks)

    #process_barcodes(reads_dir,wks)
