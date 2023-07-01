import os
from barcode import get_barcode
from gspread.utils import rowcol_to_a1
import gspread

def parse_abi_name(abi_name):
    abi_name = abi_name.split(".ab1")[0]
    abi_tokens = reversed(abi_name.split("_"))
    well,primer,colony,*plate, = abi_tokens
    plate = "_".join(reversed(plate))
    return plate,colony,primer

def parse_dirname(reads_dir):
    name = reads_dir.split("/")[-1]
    print(name)
    print("---"*15)
    email,order_id,*tail = name.split("_")
    if len(tail) == 2:
        re,seq_date = tail
        rerun = "TRUE"
        print("--rerun")
    else:
        seq_date = tail[0]
        rerun = ''
    return email,order_id,seq_date,rerun

def init_worksheet():
    gc = gspread.oauth()
    wks = gc.open("Preps")
    return wks

def process_folder(folder,wks):
    reads_dir = f"{seq_dir}/{folder}"
    process_barcodes(reads_dir,wks)

    processed_sheet = wks.worksheet("processed_folders")
    folder_cols = processed_sheet.find("folder_name").col
    folders = processed_sheet.col_values(folder_cols)
    next_row = len(folders) + 1
    processed_sheet.update_cell(next_row,folder_cols,folder)


def process_barcodes(reads_dir,wks,verbose=False):
    _,order_id,seq_date,rerun = parse_dirname(reads_dir)
    abis = sorted([f"{reads_dir}/{f}" for f in os.listdir(reads_dir) if f.endswith(".ab1") and not f.startswith(".")])
    items = len(abis)

    print(f"Processing {items} items...")
    print("---"*15)

    sheet = wks.worksheet("sequencing")
    first_col = sheet.find("Plate").col
    last_col = sheet.find("seq_date").col

    rows = []
    for abi in abis:
        bcode = get_barcode(abi)
        abi_name = abi.split("/")[-1]
        row = [*parse_abi_name(abi_name),rerun,bcode,order_id,seq_date]
        print(row)
        rows.append(row)
        if verbose:
            for item in row:
                print(item)
            print("---")

    next_row = len(sheet.col_values(first_col))+1
    last_row = next_row + items
    first_cell = rowcol_to_a1(next_row,first_col)
    last_cell = rowcol_to_a1(last_row,last_col)
    cell_range = f"{first_cell}:{last_cell}"
    
    sheet.update(cell_range,rows)
    

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
