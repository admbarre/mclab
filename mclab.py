from openpyxl import load_workbook
from dataclasses import dataclass, astuple
import datetime
import sys
import pandas as pd

# This setup means that our order info needs to be preprocessed before being sent in here...?
# TODO: thinkab
class McOrder:
    # Static fields... I think
    mcform_template = "res/mcform.xlsx"
    cells = { "date": "J3",
              "name": "D5",
              "institution": "D6",
              "address": "D7",
              "city_state_zip": "D8",
              "phone": "D9",
              "fax": "D10",
              "email": "D11",
              "po": "I5"}
    order_sheet = "Sheet1"

    GC_RICH = ["p278"]

    def __init__(self,rxns_info,contact_info):
        # Instance fields
        self.rxns_info = rxns_info
        self.contact_info = contact_info
        self.workbook = load_workbook(filename=mcform_template)
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

    # TODO: do we just access our instance fields instead of
    # passing in as parameter?
    def add_rxns(self,info):
        rxns = self.parse_rxns(info)
        rxn_row_start = 24

        for i,rxn in enumerate(rxns):
            row = rxn_row_start + i

            # RXN Cells C --> J
            # TEMPLATE  | PRIMER
            # Name Conc | Name Conc Tm DNA_type Size Notes
            rxn_cells = [f"{col}{row}" for col in "CDEFGHIJ"]
            for cell,entry in zip(rxn_cells,rxn):
                self.sheet[cell] = entry


    # TODO: provide an easier to read format
    # google sheets currently has this column to simplify
    def mc_dilute(self,num_rxns,conc,target_conc=100):
        # NOTE:
        # primers they want 3.2uM - 5uM
        # 1ul per rxn but 3ul for safety
        mc_min_vol = 3 * num_rxns # 1ul per rxn, mclab wants 3x to be safe
        our_min_vol = 10 # we dilute into 10ul H2O
        min_min = max(mc_min_vol,our_min_vol)

        dilution_factor = (target_conc)/(conc-target_conc)
        v1 = dilution_factor * min_min
        v1 = round(v1,2)
        return f"+ add {v1}ul into {min_min}ul of H2O"

    # Provide converted rxns object from JSON
    def parse_rxns(self, rxns_info):
        # TODO: make this read from an actual primers sheet in database
        melttemps = self.load_meltemps("res/meltemps.csv")
        rxns = []
        unique_primers = set()
        for plasmid,info in rxns_info.items():
            conc,primers = info.values()


            # TODO: turn this into a separate function and print out an instruction sheet
            # TODO: this could be calculated on the JS side and added into the JSON, and we can keep track of this
            if conc != "":
                num_rxns = len(primers)
                dilute_instructions = self.mc_dilute(num_rxns,conc)
                print(f"{plasmid}: ")
                print(dilute_instructions)


            for p in primers:
                # TODO: convert this section into a primer module or something
                unique_primers.add(p) # trying to make a primer shopping list
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
                      t_conc = 100,
                      p_name = p,
                      p_conc = 3.2,
                      p_melt = melt,
                      t_type = "plasmid",
                      t_size = 10,
                      p_handle = handle)

                rxns.append(rxn)
        print("Get these primers out for use:")
        print(unique_primers)
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

@dataclass
class mcRxn():
    t_name: str
    t_conc: int # 100 ng/ul
    p_name: str
    p_conc: float # 3.2  uM
    p_melt: int
    t_type: str # "plasmid"
    t_size: int # 10 kB

    p_handle: str # empty usually

    def __iter__(self):
        return iter(astuple(self))
