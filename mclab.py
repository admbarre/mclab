from openpyxl import load_workbook
from dataclasses import dataclass, astuple
from collections import Counter

import datetime
import sys
import pandas as pd

# This setup means that our order info needs to be preprocessed before being sent in here...?
# TODO: thinkab
class McOrder:
    # Static fields... I think
    mcform_template = "res/mcform.xlsx"
    cells = { 
        "date": "J3",
        "name": "D5",
        "institution": "D6",
        "address": "D7",
        "city_state_zip": "D8",
        "phone": "D9",
        "fax": "D10",
        "email": "D11",
        "po": "I5"}
    info_cols = {
        "C": "t_name",
        "D": "t_conc",
        "E": "p_name",
        "F": "p_conc",
        "G": "p_melt",
        "H": "t_type",
        "I": "t_size",
        "J": "p_handle"}
    order_sheet = "Sheet1"

    # TODO: this should be calculated in primer sheet
    # McLab said they run everything normally anyway and only do GC rich 
    # if run fails, so this might be redundant
    GC_RICH = ["p278","p85"]

    def __init__(self,rxns_info,contact_info):
        # Instance fields
        self.rxns_info = rxns_info
        self.contact_info = contact_info
        self.workbook = load_workbook(filename=self.mcform_template)
        self.sheet = self.workbook[McOrder.order_sheet]

        self.generate_order()


    def generate_order(self):
        self.today = datetime.date.today().strftime("%m-%d-%Y")
        self.add_contact()
        self.add_rxns(self.rxns_info)

    def add_contact(self):
        for k,v in self.contact_info.items():
            if k in McOrder.cells:
                self.sheet[McOrder.cells[k]] = v
        self.sheet[McOrder.cells["date"]] = self.today

    # TODO: do we just access our instance fields instead of
    # passing in as parameter?
    def add_rxns(self,info):
        rxns = self.parse_rxns(info)

        # this is the row on the spreadsheet where rxns start
        rxn_row_start = 24 
        for i,rxn in enumerate(rxns):
            row = rxn_row_start + i
            rxn_cells = [(f"{col}{row}",value) for col,value in self.info_cols.items()]
            # TODO: i do not like my naming here, rethink
            for cell,value in rxn_cells:
                self.sheet[cell] = getattr(rxn,value) # is there a better way?

    # Provide converted rxns object from JSON
    def parse_rxns(self, rxns_info):
        # TODO: make this read from an actual primers sheet in database
        # TODO: add counts to see how much primer we need to provide
        melttemps = self.load_meltemps("res/meltemps.csv")
        rxns = []
        primers_count = []
        for plasmid,info in rxns_info.items():
            conc,primers = info.values()

            for p in primers:
                # TODO: convert this section into a primer module or something
                primers_count.append(p)
                try:
                    melt = melttemps[p]
                except:
                    print("Melting temp not available, put in manually")
                    melt = -1

                if p in McOrder.GC_RICH:
                    handle = "GC Rich"
                else:
                    handle = ""

                rxn = mcRxn(t_name = plasmid,
                      p_name = p,
                      p_melt = melt,
                      p_handle = handle)

                rxns.append(rxn)
        primers_count = Counter(primers_count)
        for p in primers_count:
            rxn_count = primers_count[p]
            vol_needed = rxn_count * 3 # 1ul per rxn, McLab asked for 3x margin
            print(f"Primer: {p} minimum volume: {vol_needed}ul")
            
        print(f"Reactions: {len(rxns)}")
        return rxns

    def load_meltemps(self,filename):
        primer_temps = pd.read_csv(filename)
        # Sorting out just the columns we want
        primer_temps = primer_temps[["primer_id","mt"]]
        temp_dict_list = primer_temps.to_dict("records")
        temp_dict = {}
        for item in temp_dict_list:
            primer_id = item["primer_id"]
            primer_id = f"p{primer_id}"
            mt = item["mt"]
            temp_dict[primer_id] = mt
        return temp_dict

    def save(self,directory):
        name = f"{directory}/{self.today} order.xlsx"
        self.workbook.save(filename=name)

        # TODO: return link to saved file location
        return name

# NOTE: the order here doesn't match the spreadsheet order
# because default valued args go last
@dataclass
class mcRxn():
    t_name: str
    p_name: str
    p_melt: int
    p_handle: str # empty usually
    t_conc: int = 100 # ng/ul
    p_conc: float = 3.2  # uM
    t_type: str = "plasmid"
    t_size: int  = 10 # kB


    def __iter__(self):
        return iter(astuple(self))

if __name__ == "__main__":
    print("If you wish to place an order use mcorder.py.")
    print("If you wish to download data use mcdown.py.")
    print("If you wish to download and process barcodes use process_results.py")
