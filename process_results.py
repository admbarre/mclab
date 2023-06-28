import sys
import time

import os
from mcdown import main
from barcode import get_barcode
import gspread

req_per_min = 60 # i think
def parse_abi_name(abi_name):
    abi_name = abi_name.split(".ab1")[0]
    abi_tokens = reversed(abi_name.split("_"))
    well,primer,*name = abi_tokens
    name = "_".join(reversed(name))
    return name

def parse_dirname(reads_dir):
    name = reads_dir.split("/")[-1]
    email,order_id,date = name.split("_")
    return email,order_id,date

def process_barcodes(reads_dir,verbose=True,throttle=True):
    _,order_id,date = parse_dirname(reads_dir)
    abis = sorted([f"{reads_dir}/{f}" for f in os.listdir(reads_dir) if f.endswith(".ab1") and not f.startswith(".")])

    items = len(abis)
    wait_time = items/req_per_min * 3 # might need some more margin... 

    print(f"Processing {items} items...")
    print("application maybe throttled, idk how to do this ideally with rate limiting")
    print(" ¯\_(ツ)_/¯")
    print(f"wait time in between writes: {wait_time}")

    gc = gspread.oauth()
    wks = gc.open("Preps")
    sheet = wks.sheet1
    barcode_col = sheet.find("Barcode").col
    order_col = sheet.find("order_id").col
    seq_date_col = sheet.find("seq_date").col

    for abi in abis:
        bcode = get_barcode(abi)
        if bcode == None:
            bcode = "not found"

        abi_name = abi.split("/")[-1]
        name = parse_abi_name(abi_name)
        row = sheet.find(name).row

        if verbose:
            print(name)
            print(row)
            print(bcode)
            print("---")
        if throttle:
            time.sleep(wait_time)

        sheet.update_cell(row,barcode_col, bcode)
        # TODO: this can be handled in a single pass without iterating
        # TODO: these all probably count as separate writes for the quota too
        # damn itttt
        # theres gotta be a way to batch it...
        sheet.update_cell(row,order_col,order_id)
        sheet.update_cell(row,seq_date_col,date)


# TODO: might need to make it interactive to select which to download especially
# if they re run samples for us
if __name__ == "__main__":
    if len(sys.argv) == 2:
        _,credentials_path = sys.argv
        reads_dir = main(credentials_path)
        process_barcodes(reads_dir)
